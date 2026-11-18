#pragma once

#include "fusion.hpp"

struct Overlap {
    double intersection_area = 0.0;
    double iou = 0.0;
    bool height_overlap = false;
};

/** Вычисляет площадь пересечения ориентированных прямоугольников в ENU.
 * @param a Первый объект с координатой основания, габаритами и курсом.
 * @param b Второй объект в той же системе координат.
 * @return Площадь в м², IoU в [0,1] и признак пересечения интервалов высоты.
 */
Overlap footprintOverlap(const Track3D& a, const Track3D& b);
