#include "overlap.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {
struct Point { double x, y; };
using Polygon = std::vector<Point>;

/** Находит знаковую площадь параллелограмма ABP.
 * @param a Начало ориентированного ребра.
 * @param b Конец ребра.
 * @param p Проверяемая точка.
 * @return Положительное значение для точки слева от AB.
 */
double side(const Point& a, const Point& b, const Point& p) {
    return (b.x-a.x)*(p.y-a.y) - (b.y-a.y)*(p.x-a.x);
}

/** Строит контур основания против часовой стрелки.
 * @param size Ширина, высота и длина в метрах.
 * @param center Центр основания на плоскости сравнения.
 * @param heading Поворот от продольной оси этой плоскости по часовой стрелке, рад.
 * @return Четыре вершины горизонтальной проекции.
 */
Polygon rectangle(const Vec3& size, Point center, double heading) {
    const double c = std::cos(heading), s = std::sin(heading);
    const double half_width = size.x / 2, half_length = size.z / 2;
    const std::array<Point,4> local = {{{-half_width,-half_length}, {half_width,-half_length},
                                       {half_width,half_length}, {-half_width,half_length}}};
    Polygon result;
    result.reserve(local.size());
    for (const auto& point : local) {
        result.push_back({center.x + c*point.x + s*point.y,
                          center.y - s*point.x + c*point.y});
    }
    return result;
}

/** Отсекает выпуклый многоугольник одним ребром контура.
 * @param polygon Текущая область пересечения.
 * @param a Начало отсекающего ребра.
 * @param b Конец ребра. Сохраняется левая полуплоскость.
 * @param roundoff Погрешность координат после поворота, м.
 * @return Оставшаяся часть многоугольника.
 */
Polygon clip(const Polygon& polygon, const Point& a, const Point& b, double roundoff) {
    Polygon result;
    if (polygon.empty()) return result;
    Point previous = polygon.back();
    double before = side(a,b,previous);
    const double edge_error = roundoff * (std::abs(b.x - a.x) + std::abs(b.y - a.y));
    bool has_interior = false;
    for (const auto& current : polygon) {
        const double after = side(a,b,current);
        if (after > edge_error) has_interior = true;
        if ((before >= 0) != (after >= 0)) {
            const double ratio = before/(before-after);
            result.push_back({previous.x+ratio*(current.x-previous.x),
                              previous.y+ratio*(current.y-previous.y)});
        }
        if (after >= 0) result.push_back(current);
        previous = current;
        before = after;
    }
    // При касании все вершины снаружи или на ребре с точностью поворота осей.
    if (!has_interior) return {};
    return result;
}

/** Суммирует площади треугольников с общей вершиной выпуклого контура.
 * @param polygon Вершины в порядке обхода.
 * @return Неотрицательная площадь; для отрезка и точки возвращается 0.
 */
double area(const Polygon& polygon) {
    double sum = 0;
    for (std::size_t i = 1; i + 1 < polygon.size(); ++i) {
        sum += side(polygon.front(), polygon[i], polygon[i + 1]);
    }
    return std::abs(sum)/2;
}

/** Проверяет данные, необходимые для построения прямоугольника.
 * @param object Объект в ENU с положительными габаритами.
 * @throws std::invalid_argument При неконечных значениях или нулевом размере.
 */
void validateShape(const Track3D& object) {
    const auto& p = object.position;
    const auto& s = object.size;
    if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z) ||
        !std::isfinite(object.heading) || !std::isfinite(s.x) || !std::isfinite(s.y) ||
        !std::isfinite(s.z) || s.x <= 0 || s.y <= 0 || s.z <= 0 ||
        !std::isfinite(s.x*s.z) || s.x*s.z == 0) {
        throw std::invalid_argument("Overlap requires finite coordinates, heading and positive dimensions");
    }
}
} // namespace

/** Сравнивает горизонтальные проекции объектов с учётом их поворота.
 * @param a Первый объект; position.z обозначает высоту основания.
 * @param b Второй объект в ENU.
 * @return IoU = Sпересечения/(Sa+Sb-Sпересечения), а также пересечение по высоте.
 */
Overlap footprintOverlap(const Track3D& a, const Track3D& b) {
    validateShape(a);
    validateShape(b);
    Overlap result;
    result.height_overlap = std::min(a.position.z+a.size.y,b.position.z+b.size.y) >
                            std::max(a.position.z,b.position.z);

    // Окружности вокруг прямоугольников позволяют пропустить явно далёкие пары.
    const double radius_a = std::hypot(a.size.x,a.size.z)/2;
    const double radius_b = std::hypot(b.size.x,b.size.z)/2;
    const double dx = a.position.x-b.position.x, dy = a.position.y-b.position.y;
    if (dx*dx+dy*dy >= (radius_a+radius_b)*(radius_a+radius_b)) return result;

    // В осях меньшего основания его вершины точны даже при большом общем курсе.
    // При одинаковом курсе второй контур тоже остаётся параллельным этим осям.
    const double area_a = a.size.x * a.size.z, area_b = b.size.x * b.size.z;
    const auto& first = area_a <= area_b ? a : b;
    const auto& second = area_a <= area_b ? b : a;
    const Vec3 delta = second.position - first.position;
    const double c = std::cos(first.heading), s = std::sin(first.heading);
    const Point center = {c*delta.x - s*delta.y, s*delta.x + c*delta.y};
    const double half_turn = std::acos(-1.0);
    const double angle = std::remainder(second.heading - first.heading, half_turn);
    const double scale = std::max({std::abs(delta.x), std::abs(delta.y),
                                   first.size.x, first.size.z, second.size.x, second.size.z});
    const double roundoff = (first.heading == 0 && angle == 0) ? 0
        : 4 * std::numeric_limits<double>::epsilon() * scale;
    const Polygon boundary = rectangle(second.size, center, angle);
    Polygon intersection = rectangle(first.size, {0, 0}, 0);
    for (std::size_t i = 0; i < boundary.size() && !intersection.empty(); ++i) {
        intersection = clip(intersection, boundary[i], boundary[(i + 1) % boundary.size()], roundoff);
    }
    result.intersection_area = std::min(area(intersection),std::min(area_a,area_b));
    const double union_area = std::max(area_a, area_b) + (std::min(area_a, area_b) - result.intersection_area);
    result.iou = result.intersection_area / union_area;
    return result;
}
