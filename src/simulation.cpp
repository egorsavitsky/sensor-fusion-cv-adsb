#include "simulation.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

namespace {
/** Проверяет допустимость параметров модели АЗН-В.
 * @param config Период, ошибки, задержка, потери, охват объектов и seed модели.
 */
void validate(const SimulationConfig& config) {
    if (!std::isfinite(config.period) || config.period <= 0.0 ||
        !std::isfinite(config.horizontal_sigma) || config.horizontal_sigma < 0.0 ||
        !std::isfinite(config.vertical_sigma) || config.vertical_sigma < 0.0 ||
        !std::isfinite(config.delay) || config.delay < 0.0 ||
        !std::isfinite(config.dropout) || config.dropout < 0.0 || config.dropout > 1.0 ||
        !std::isfinite(config.coverage) || config.coverage < 0.0 || config.coverage > 1.0 ||
        !std::isfinite(config.heading_sigma) || config.heading_sigma < 0.0) {
        throw std::invalid_argument("Invalid ADS-B simulation parameters");
    }
}
} // namespace

/** Создаёт опорное описание объекта для моделирования второго источника.
 * @param obstacle Габариты, класс и опорное положение из лога.
 * @param time Время кадра в секундах от начала записи.
 * @param frame Общая система ENU.
 * @return Трек с координатами, габаритами и направлением в ENU.
 */
Track3D referenceTrack(const Obstacle& obstacle, double time, const LocalFrame& frame) {
    Track3D result = {time,obstacle.id,frame.fromGeodetic(obstacle.reference)};
    result.size = obstacle.size;
    result.object_class = obstacle.object_class;
    const Vec3 direction = frame.vectorFromEcef(enuToEcef({std::sin(obstacle.reference_heading),
                                                         std::cos(obstacle.reference_heading),0},obstacle.reference));
    result.heading = std::atan2(direction.x,direction.y);
    return result;
}

/** Переводит измеренный камерой объект в общую систему координат.
 * @param obstacle Положение центра, направление и габариты из СТЗ.
 * @param time Время кадра в секундах.
 * @param transform Подготовленное преобразование камеры для данного кадра.
 * @return Трек с координатой основания и курсом в ENU.
 */
Track3D cameraTrack(const Obstacle& obstacle, double time, const CameraTransform& transform) {
    Track3D result = {time,obstacle.id,transform.apply(obstacle.position,obstacle.size.y)};
    result.size = obstacle.size;
    result.object_class = obstacle.object_class;
    const Vec3 direction = transform.direction(obstacle.direction);
    if (std::hypot(direction.x,direction.y) < 1e-12) {
        throw std::invalid_argument("Object direction has no horizontal component");
    }
    result.heading = std::atan2(direction.x,direction.y);
    return result;
}

/** Создаёт небольшой сценарий с самолётом и машиной для проверки пересечения.
 * @param camera Известное смещение и углы камеры для прямого и обратного преобразования.
 * @return 12 кадров по 0,5 с с двумя движущимися объектами разных габаритов и курсов.
 */
LogData twoObjectScenario(const Camera& camera) {
    const Geodetic origin = {radians(43.4),radians(39.9),27};
    const LocalFrame local(origin);
    LogData result;
    for (int index=0; index<12; ++index) {
        Frame frame;
        frame.timestamp = index*0.5;
        const double time = frame.timestamp;
        frame.own = local.toGeodetic({0,time*10,0});
        frame.attitude = {radians(28.64+time*4),0,0};
        for (int object=0; object<2; ++object) {
            Obstacle obstacle;
            obstacle.id = object+1;
            obstacle.object_class = object;
            obstacle.size = object == 0 ? Vec3{12,4,22} : Vec3{2.2,1.8,4.8};
            const Vec3 position = object == 0 ? Vec3{time*5,100+time*4,1+time*.2}
                                              : Vec3{35-time,110+time,0};
            const double heading = object == 0 ? radians(20+time*5) : radians(-45);
            obstacle.reference = local.toGeodetic(position);

            // Сначала задаются опорные объекты, затем вычисляется то, что видит камера.
            Vec3 relative = ecefToEnu(toEcef(obstacle.reference)-toEcef(frame.own),frame.own);
            if (camera.position_is_center) relative.z += obstacle.size.y/2;
            obstacle.position = bodyToCamera(enuToBody(relative,frame.attitude)-camera.offset,camera.mounting);
            obstacle.position.x += .3*std::sin(time+object);
            obstacle.position.y += .02*std::sin(time);
            obstacle.position.z += .3*std::cos(time+object);
            const Vec3 world_direction = enuToEcef({std::sin(heading),std::cos(heading),0},origin);
            const Vec3 ref_direction = ecefToEnu(world_direction,obstacle.reference);
            obstacle.reference_heading = std::atan2(ref_direction.x,ref_direction.y);
            obstacle.direction = bodyToCamera(enuToBody(ecefToEnu(world_direction,frame.own),
                                                       frame.attitude),camera.mounting);
            frame.obstacles.push_back(obstacle);
            ++result.records;
        }
        result.frames.push_back(frame);
    }
    return result;
}

/** Формирует независимые ID и сообщения АЗН-В с ошибками и потерями.
 * @param frames Кадры в порядке возрастания времени.
 * @param frame Общая локальная система ENU.
 * @param config Период, ошибки, задержка, потери, охват объектов и seed модели.
 * @return Модельные сообщения и отдельная таблица опорных ID для оценки результата.
 */
ModelData simulateAdsb(const std::vector<Frame>& frames, const LocalFrame& frame,
                       const SimulationConfig& config) {
    validate(config);
    std::map<int, std::vector<Track3D>> references;
    for (const auto& item : frames) {
        for (const auto& obstacle : item.obstacles) {
            references[obstacle.id].push_back(referenceTrack(obstacle,item.timestamp,frame));
        }
    }
    std::mt19937 coverage_random(config.seed);
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    ModelData result;
    double scheduled = 0;
    int model_id = 100000;
    for (const auto& entry : references) {
        const int adsb_id = model_id++;
        if (uniform(coverage_random) >= config.coverage) continue;
        result.reference_ids[adsb_id] = entry.first;
        std::seed_seq loss_seed{config.seed,static_cast<std::uint32_t>(adsb_id),1u};
        // Разные генераторы позволяют менять потери, сохраняя ошибки общих измерений.
        std::seed_seq noise_seed{config.seed,static_cast<std::uint32_t>(adsb_id),2u};
        std::mt19937 loss_random(loss_seed), noise_random(noise_seed);
        std::normal_distribution<double> normal(0.0, 1.0);
        const auto& truth = entry.second;
        scheduled += std::floor((truth.back().timestamp-truth.front().timestamp)/config.period) + 1.0;
        if (!std::isfinite(scheduled) || scheduled > 10000000) {
            throw std::invalid_argument("More than 10 million synthetic samples; increase --period");
        }
        auto& messages = result.streams[adsb_id];
        std::size_t next = 0;
        for (std::size_t sample = 0; ; ++sample) {
            const double time = truth.front().timestamp + sample * config.period;
            if (time > truth.back().timestamp) break;
            while (next + 1 < truth.size() && truth[next + 1].timestamp <= time) ++next;
            Track3D measurement = truth[next];
            if (next + 1 < truth.size() && truth[next].timestamp != time) {
                const double gap = truth[next+1].timestamp-truth[next].timestamp;
                if (!interpolateAdsb(time,truth[next],truth[next+1],gap+1e-9,measurement)) {
                    throw std::runtime_error("Cannot interpolate reference object");
                }
            }
            Vec3& position = measurement.position;
            position.x += normal(noise_random) * config.horizontal_sigma;
            position.y += normal(noise_random) * config.horizontal_sigma;
            position.z += normal(noise_random) * config.vertical_sigma;
            measurement.heading += normal(noise_random) * config.heading_sigma;
            if (uniform(loss_random) < config.dropout) continue;
            measurement.timestamp = time;
            measurement.id = adsb_id;
            messages.push_back({measurement, time + config.delay});
            ++result.messages;
        }
    }
    return result;
}

/** Задаёт режим согласования времени и ограничения устаревания.
 * @param model Модельные потоки сообщений, упорядоченные по времени.
 * @param max_age Максимальный возраст последнего измерения в секундах.
 * @param max_gap Максимальный интервал между измерениями в секундах.
 * @param offline Разрешение использовать будущие измерения для интерполяции.
 */
AdsbSynchronizer::AdsbSynchronizer(const ModelData& model, double max_age,
                                   double max_gap, bool offline)
    : model_(model), max_age_(max_age), max_gap_(max_gap), offline_(offline) {
    if (!std::isfinite(max_age) || max_age < 0.0 || !std::isfinite(max_gap) || max_gap <= 0.0) {
        throw std::invalid_argument("Invalid synchronization limits");
    }
}

/** Получает доступные положения АЗН-В на очередной кадр.
 * @param time Момент оценки в секундах от начала записи.
 * @return Треки на время кадра; устаревшие измерения исключены.
 */
std::vector<Track3D> AdsbSynchronizer::at(double time) {
    if (!std::isfinite(time) || time < last_time_) {
        throw std::invalid_argument("Synchronization requires ordered finite times");
    }
    last_time_ = time;
    std::vector<Track3D> tracks;
    tracks.reserve(model_.streams.size());
    for (const auto& entry : model_.streams) {
        const auto& messages = entry.second;
        std::size_t& next = cursors_[entry.first];
        // Причинный режим ограничен временем доставки; офлайн — временем измерения.
        while (next < messages.size() &&
               (offline_ ? messages[next].track.timestamp : messages[next].received_at) <= time) ++next;
        if (next == 0) continue;
        const Track3D& previous = messages[next-1].track;
        Track3D current;
        if (offline_ && next < messages.size() && previous.timestamp != time &&
            interpolateAdsb(time, previous, messages[next].track, max_gap_, current)) {
            tracks.push_back(current);
        } else if (time - previous.timestamp <= max_age_) {
            // Для причинного режима используем последнее принятое положение без прогноза.
            current = previous;
            current.timestamp = time;
            tracks.push_back(current);
        }
    }
    return tracks;
}
