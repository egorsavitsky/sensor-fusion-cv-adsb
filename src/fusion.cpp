#include "fusion.hpp"
#include "overlap.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

/** Складывает координаты двух векторов.
 * @param a Первое слагаемое.
 * @param b Второе слагаемое.
 * @return Покомпонентная сумма.
 */
Vec3 operator+(const Vec3& a, const Vec3& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

/** Вычитает координаты второго вектора из первого.
 * @param a Уменьшаемое.
 * @param b Вычитаемое.
 * @return Покомпонентная разность.
 */
Vec3 operator-(const Vec3& a, const Vec3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

/** Умножает вектор на число.
 * @param a Исходный вектор.
 * @param scale Множитель для каждой координаты вектора.
 * @return Масштабированный вектор.
 */
Vec3 operator*(const Vec3& a, double scale) {
    return {a.x * scale, a.y * scale, a.z * scale};
}

/** Вычисляет квадрат евклидова расстояния в трёх измерениях.
 * @param a Координаты первой точки.
 * @param b Координаты второй точки в тех же осях и единицах.
 * @return Квадрат расстояния без извлечения корня.
 */
double squaredDistance(const Vec3& a, const Vec3& b) {
    const Vec3 d = a - b;
    return d.x * d.x + d.y * d.y + d.z * d.z;
}

namespace {
struct Candidate {
    std::size_t cv;
    std::size_t adsb;
    double cost;
};

/** Проверяет координаты, размеры, время и уникальность ID внутри источника.
 * @param tracks Объекты одного источника на одной временной метке.
 */
void checkTracks(const std::vector<Track3D>& tracks) {
    std::vector<int> ids;
    ids.reserve(tracks.size());
    for (const auto& track : tracks) {
        const Vec3& p = track.position;
        if (!std::isfinite(track.timestamp) || !std::isfinite(p.x) ||
            !std::isfinite(p.y) || !std::isfinite(p.z) || !std::isfinite(track.heading) ||
            !std::isfinite(track.size.x) || !std::isfinite(track.size.y) || !std::isfinite(track.size.z) ||
            track.size.x <= 0 || track.size.y <= 0 || track.size.z <= 0) {
            throw std::invalid_argument("Tracks need finite time, position, heading and positive finite dimensions");
        }
        ids.push_back(track.id);
        if (std::abs(track.timestamp - tracks.front().timestamp) > 1e-6) {
            throw std::invalid_argument("Tracks must refer to the same frame time");
        }
    }
    std::sort(ids.begin(), ids.end());
    if (std::adjacent_find(ids.begin(), ids.end()) != ids.end()) {
        throw std::invalid_argument("Duplicate track ID in one frame");
    }
}

/** Проверяет согласование времени двух источников перед сопоставлением.
 * @param cv Объекты СТЗ в общей системе ENU.
 * @param adsb Объекты АЗН-В на то же время.
 */
void checkTimes(const std::vector<Track3D>& cv, const std::vector<Track3D>& adsb) {
    if (!cv.empty() && !adsb.empty() && std::abs(cv.front().timestamp-adsb.front().timestamp) > 1e-6) {
        throw std::invalid_argument("Synchronize sources before association");
    }
}

/** Собирает выбранные пары и списки объектов без пары.
 * @param cv Объекты СТЗ в общей системе ENU.
 * @param adsb Объекты АЗН-В на то же время.
 * @param selected Выбранные индексы пар без повторного назначения объектов.
 * @param limit Порог расстояния в метрах для оценки близости.
 * @return Результат назначения с расстояниями и оценками близости.
 */
Association collect(const std::vector<Track3D>& cv, const std::vector<Track3D>& adsb,
                    const std::vector<Candidate>& selected, double limit) {
    Association result;
    std::vector<bool> used_cv(cv.size(), false), used_adsb(adsb.size(), false);
    for (const auto& pair : selected) {
        const double distance = std::sqrt(squaredDistance(cv[pair.cv].position,
                                                          adsb[pair.adsb].position));
        result.matches.push_back({cv[pair.cv].id, adsb[pair.adsb].id,
                                  distance, 1.0 - distance / limit});
        used_cv[pair.cv] = true;
        used_adsb[pair.adsb] = true;
    }
    for (std::size_t i = 0; i < cv.size(); ++i) {
        if (!used_cv[i]) result.cv_only.push_back(cv[i].id);
    }
    for (std::size_t i = 0; i < adsb.size(); ++i) {
        if (!used_adsb[i]) result.adsb_only.push_back(adsb[i].id);
    }
    return result;
}

/** Выбирает свободные пары в порядке возрастания стоимости.
 * @param pairs Допустимые пары с индексами и стоимостью; порядок может измениться.
 * @param n Число объектов СТЗ.
 * @param m Число объектов АЗН-В.
 * @return Индексы пар без повторного использования объектов.
 */
std::vector<Candidate> greedy(std::vector<Candidate>& pairs, std::size_t n, std::size_t m) {
    std::sort(pairs.begin(), pairs.end(), [](const Candidate& a, const Candidate& b) {
        if (a.cost != b.cost) return a.cost < b.cost;
        if (a.cv != b.cv) return a.cv < b.cv;
        return a.adsb < b.adsb;
    });
    std::vector<bool> used_cv(n, false), used_adsb(m, false);
    std::vector<Candidate> selected;
    selected.reserve(std::min(n, m));
    for (const auto& pair : pairs) {
        if (!used_cv[pair.cv] && !used_adsb[pair.adsb]) {
            selected.push_back(pair);
            used_cv[pair.cv] = true;
            used_adsb[pair.adsb] = true;
        }
    }
    return selected;
}

/** Решает задачу назначения венгерским алгоритмом с возможностью отказа.
 * @param n Число объектов СТЗ.
 * @param real_columns Число реальных объектов АЗН-В без столбцов отказа.
 * @param candidates Допустимые пары. Стоимость равна −IoU либо нормированному квадрату расстояния.
 * @return Назначение с наибольшим числом пар и наименьшей суммарной стоимостью.
 */
std::vector<Candidate> optimal(std::size_t n, std::size_t real_columns,
                               const std::vector<Candidate>& candidates) {
    const std::size_t m = real_columns + n;
    if (n == 0) return {};
    // В каждом из двух режимов разброс цен не больше 1. Штраф n+1 даёт приоритет числу пар.
    const double unmatched = static_cast<double>(n + 1);
    const double forbidden = unmatched * 2.0;
    std::vector<std::vector<double>> cost(n, std::vector<double>(m, unmatched));
    for (std::size_t i = 0; i < n; ++i) {
        std::fill(cost[i].begin(),cost[i].begin()+real_columns,forbidden);
    }
    for (const auto& pair : candidates) cost[pair.cv][pair.adsb] = pair.cost;
    std::vector<double> u(n + 1), v(m + 1);
    std::vector<std::size_t> p(m + 1), way(m + 1);
    for (std::size_t i = 1; i <= n; ++i) {
        // Ищем увеличивающий путь с минимальной приведённой стоимостью.
        p[0] = i;
        std::size_t column = 0;
        std::vector<double> best(m + 1, std::numeric_limits<double>::infinity());
        std::vector<bool> used(m + 1, false);
        do {
            used[column] = true;
            const std::size_t row = p[column];
            double delta = std::numeric_limits<double>::infinity();
            std::size_t next = 0;
            for (std::size_t j = 1; j <= m; ++j) {
                if (used[j]) continue;
                const double value = cost[row - 1][j - 1] - u[row] - v[j];
                if (value < best[j]) {
                    best[j] = value;
                    way[j] = column;
                }
                if (best[j] < delta) {
                    delta = best[j];
                    next = j;
                }
            }
            for (std::size_t j = 0; j <= m; ++j) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    best[j] -= delta;
                }
            }
            column = next;
        } while (p[column] != 0);
        do {
            const std::size_t previous = way[column];
            p[column] = p[previous];
            column = previous;
        } while (column != 0);
    }
    std::vector<Candidate> selected;
    for (std::size_t j = 1; j <= real_columns; ++j) {
        // Служебные столбцы отказа не попадают в список соответствий.
        if (p[j] != 0 && cost[p[j] - 1][j - 1] < unmatched) {
            selected.push_back({p[j] - 1, j - 1, cost[p[j] - 1][j - 1]});
        }
    }
    return selected;
}
} // namespace

/** Интерполирует положение, размеры и курс между двумя измерениями.
 * @param time Момент оценки в секундах от начала записи.
 * @param left Измерение того же объекта слева по времени.
 * @param right Измерение справа по времени с тем же ID и классом.
 * @param max_gap Максимальный интервал между измерениями в секундах.
 * @param result Выходной трек на заданный момент времени.
 * @return true при допустимых входах; false при несовпадении ID, класса, времени или разрыве.
 */
bool interpolateAdsb(double time, const Track3D& left, const Track3D& right,
                     double max_gap, Track3D& result) {
    if (!std::isfinite(time) || !std::isfinite(max_gap) || max_gap <= 0.0 ||
        left.id != right.id || left.object_class != right.object_class || !std::isfinite(left.timestamp) ||
        !std::isfinite(right.timestamp) || right.timestamp < left.timestamp ||
        time < left.timestamp || time > right.timestamp ||
        right.timestamp - left.timestamp > max_gap) return false;
    if (!std::isfinite(squaredDistance(left.position, right.position))) return false;
    if (!std::isfinite(left.heading) || !std::isfinite(right.heading) ||
        !std::isfinite(squaredDistance(left.size,right.size)) ||
        left.size.x <= 0 || left.size.y <= 0 || left.size.z <= 0 ||
        right.size.x <= 0 || right.size.y <= 0 || right.size.z <= 0) return false;
    if (left.timestamp == right.timestamp) {
        if (squaredDistance(left.position, right.position) != 0.0 ||
            squaredDistance(left.size,right.size) != 0 || left.heading != right.heading) return false;
        result = left;
        return true;
    }
    const double ratio = (time - left.timestamp) / (right.timestamp - left.timestamp);
    result = left;
    result.timestamp = time;
    result.position = left.position + (right.position - left.position) * ratio;
    result.size = left.size + (right.size-left.size)*ratio;
    const double angle = std::atan2(std::sin(right.heading-left.heading),
                                   std::cos(right.heading-left.heading));
    result.heading = left.heading + angle*ratio;
    return true;
}

/** Задаёт предельное расстояние между точками в режиме distance.
 * @param max_distance Порог в метрах. Пара на самой границе исключается.
 */
SensorFusion::SensorFusion(double max_distance) : max_distance_(max_distance) {
    if (!std::isfinite(max_distance) || max_distance <= 0.0 ||
        !std::isfinite(max_distance * max_distance) || max_distance * max_distance == 0.0) {
        throw std::invalid_argument("Distance threshold must be finite and positive");
    }
}

/** Сопоставляет объекты по трёхмерному расстоянию.
 * @param cv Объекты СТЗ в общей системе ENU.
 * @param adsb Объекты АЗН-В на то же время.
 * @param method Жадное или глобальное назначение один к одному.
 * @return Назначенные пары и объекты без соответствия.
 */
Association SensorFusion::associateTracks(const std::vector<Track3D>& cv,
                                          const std::vector<Track3D>& adsb,
                                          Assignment method) const {
    checkTracks(cv);
    checkTracks(adsb);
    checkTimes(cv, adsb);
    const double limit_sq = max_distance_ * max_distance_;
    std::vector<Candidate> pairs;
    for (std::size_t i = 0; i < cv.size(); ++i) {
        for (std::size_t j = 0; j < adsb.size(); ++j) {
            const double dx = cv[i].position.x - adsb[j].position.x;
            if (std::abs(dx) >= max_distance_) continue;
            const double dy = cv[i].position.y - adsb[j].position.y;
            if (std::abs(dy) >= max_distance_) continue;
            const double dz = cv[i].position.z - adsb[j].position.z;
            if (std::abs(dz) >= max_distance_) continue;
            const double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 < limit_sq) pairs.push_back({i, j, d2/limit_sq});
        }
    }
    const auto chosen = method == Assignment::Optimal ? optimal(cv.size(),adsb.size(),pairs)
                                                     : greedy(pairs,cv.size(),adsb.size());
    return collect(cv,adsb,chosen,max_distance_);
}

/** Назначает пары по IoU горизонтальных проекций.
 * @param cv Объекты СТЗ в ENU с габаритами и курсом.
 * @param adsb Объекты АЗН-В на тот же момент времени.
 * @param min_iou Включительный порог [0,1]; нулевая площадь всегда исключается.
 * @param method Жадное либо глобальное назначение один к одному.
 * @return Назначенные пары, объекты без пары и все положительные пересечения.
 */
Association SensorFusion::associateByOverlap(const std::vector<Track3D>& cv,
                                             const std::vector<Track3D>& adsb,
                                             double min_iou, Assignment method) const {
    checkTracks(cv);
    checkTracks(adsb);
    checkTimes(cv,adsb);
    if (!std::isfinite(min_iou) || min_iou < 0 || min_iou > 1) {
        throw std::invalid_argument("Overlap threshold must be in [0,1]");
    }
    std::vector<Candidate> candidates;
    std::vector<OverlapPair> overlaps;
    std::vector<std::vector<std::size_t>> indices(cv.size(),std::vector<std::size_t>(adsb.size()));
    for (std::size_t i=0; i<cv.size(); ++i) {
        for (std::size_t j=0; j<adsb.size(); ++j) {
            const Overlap area = footprintOverlap(cv[i],adsb[j]);
            if (area.iou <= 0) continue;
            const bool same_class = cv[i].object_class < 0 || adsb[j].object_class < 0 ||
                                    cv[i].object_class == adsb[j].object_class;
            const bool passes = area.iou >= min_iou;
            indices[i][j] = overlaps.size();
            overlaps.push_back({cv[i].id,adsb[j].id,area.intersection_area,area.iou,
                                passes,area.height_overlap,same_class,false});
            // Площадное пересечение не связывает разные классы и непересекающиеся высоты.
            if (passes && area.height_overlap && same_class) candidates.push_back({i, j, -area.iou});
        }
    }
    const auto chosen = method == Assignment::Optimal ? optimal(cv.size(),adsb.size(),candidates)
                                                     : greedy(candidates,cv.size(),adsb.size());
    Association result = collect(cv,adsb,chosen,max_distance_);
    for (std::size_t k=0; k<chosen.size(); ++k) {
        auto& pair = overlaps[indices[chosen[k].cv][chosen[k].adsb]];
        pair.assigned = true;
        result.matches[k].score = pair.iou;
        result.matches[k].overlap = pair.iou;
        result.matches[k].intersection_area = pair.intersection_area;
    }
    result.overlaps = std::move(overlaps);
    return result;
}

/** Выполняет исходный расчёт с pow и sqrt для контрольного замера.
 * @param cv Объекты СТЗ в общей системе ENU.
 * @param adsb Объекты АЗН-В на то же время.
 * @return Жадное назначение по расстоянию.
 */
Association SensorFusion::associateReference(const std::vector<Track3D>& cv,
                                             const std::vector<Track3D>& adsb) const {
    checkTracks(cv);
    checkTracks(adsb);
    checkTimes(cv, adsb);
    std::vector<Candidate> pairs;
    for (std::size_t i = 0; i < cv.size(); ++i) {
        for (std::size_t j = 0; j < adsb.size(); ++j) {
            const Vec3 d = cv[i].position - adsb[j].position;
            const double distance = std::sqrt(std::pow(d.x, 2) + std::pow(d.y, 2) +
                                              std::pow(d.z, 2));
            if (distance < max_distance_) {
                pairs.push_back({i, j, distance});
            }
        }
    }
    return collect(cv, adsb, greedy(pairs, cv.size(), adsb.size()), max_distance_);
}
