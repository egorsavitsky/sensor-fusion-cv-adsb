#include "overlap.hpp"
#include "geometry.hpp"
#include "report.hpp"
#include "simulation.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>

namespace {
/** Проверяет условие теста геометрии.
 * @param condition Ожидаемое истинное условие.
 * @param message Описание нарушения для вывода при ошибке.
 */
void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

/** Сравнивает два числа с заданной абсолютной погрешностью.
 * @param a Вычисленное значение.
 * @param b Аналитически ожидаемое значение.
 * @param tolerance Допустимое отклонение.
 * @param message Описание проверки.
 */
void close(double a, double b, double tolerance, const char* message) {
    check(std::abs(a-b) <= tolerance,message);
}

/** Создаёт прямоугольный объект в ENU для теста.
 * @param id Идентификатор объекта.
 * @param position Координаты основания в ENU.
 * @param size Ширина, высота и длина.
 * @param heading Курс в радианах.
 * @return Самолёт класса 0 на момент t=0.
 */
Track3D box(int id, Vec3 position = {}, Vec3 size = {2,2,2}, double heading = 0) {
    return {0,id,position,size,heading,0};
}

/** Проверяет IoU на случаях с известным аналитическим результатом. */
void analyticalGeometry() {
    const Track3D a = box(1);
    close(footprintOverlap(a,a).iou,1,1e-12,"Identical rectangles have IoU 1");
    const auto shifted = footprintOverlap(a,box(2,{1,0,0}));
    close(shifted.intersection_area,2,1e-12,"Half-overlap area");
    close(shifted.iou,1.0/3,1e-12,"Half-overlap IoU");
    close(footprintOverlap(a,box(2,{},{1,1,1})).iou,.25,1e-12,"Contained rectangle uses union, not smaller area");
    close(footprintOverlap(a,box(2,{}, {2,2,2},radians(45))).iou,1/std::sqrt(2.0),1e-12,"Square rotated by 45 degrees");
    close(footprintOverlap(box(1,{},{2,2,4}),box(2,{},{2,2,4},radians(90))).iou,1.0/3,1e-12,"Orientation changes long-rectangle IoU");
    check(footprintOverlap(a,box(2,{2,0,0})).iou == 0,"Touching edges do not have area");
    check(footprintOverlap(a,box(2,{20,0,0})).iou == 0,"Disjoint rectangles");
    const auto overhead = footprintOverlap(a,box(2,{0,0,10}));
    check(overhead.iou == 1 && !overhead.height_overlap,"Planar IoU and height filter are separate");
}

/** Проверяет симметрию и неизменность IoU при общем переносе и повороте. */
void invariants() {
    std::mt19937 random(77);
    std::uniform_real_distribution<double> position(-5,5), angle(-3,3), dimension(1,8);
    for (int trial=0; trial<80; ++trial) {
        Track3D a = box(1,{position(random),position(random),0},{dimension(random),2,dimension(random)},angle(random));
        Track3D b = box(2,{position(random),position(random),0},{dimension(random),2,dimension(random)},angle(random));
        const auto original = footprintOverlap(a,b);
        close(original.iou,footprintOverlap(b,a).iou,1e-10,"IoU symmetry");
        check(original.iou >= 0 && original.iou <= 1,"IoU stays in [0,1]");
        for (auto* item : {&a,&b}) {
            const double x = item->position.x, y = item->position.y, rotation = .7;
            item->position.x = std::cos(rotation)*x+std::sin(rotation)*y+1e8;
            item->position.y = -std::sin(rotation)*x+std::cos(rotation)*y-1e8;
            item->heading += rotation;
        }
        close(original.iou,footprintOverlap(a,b).iou,1e-7,"Common motion cannot change IoU");
    }
}

/** Проверяет порог, фильтры классов и высоты, а также неоднозначное назначение. */
void matching() {
    SensorFusion fusion;
    const Track3D a = box(1), b = box(2,{1,0,0});
    const double threshold = footprintOverlap(a,b).iou;
    auto result = fusion.associateByOverlap({a},{b},threshold);
    check(result.matches.size() == 1 && result.overlaps[0].assigned,"Threshold is inclusive");
    check(fusion.associateByOverlap({a},{b},std::nextafter(threshold,1.0)).matches.empty(),"Above-threshold rejection");
    check(fusion.associateByOverlap({a},{box(2,{2,0,0})},0).overlaps.empty(),"Zero threshold still rejects zero area");
    check(fusion.associateByOverlap({a},{a},1).matches.size() == 1,"Threshold 1 accepts identical objects");
    for (int degrees = 1; degrees < 180; ++degrees) {
        const auto rotated = box(1, {}, {2.2, 1.8, 4.8}, radians(degrees));
        check(fusion.associateByOverlap({rotated}, {rotated}, 1).matches.size() == 1,
              "Identical rotated footprints must pass threshold 1");
        auto touching = rotated;
        touching.position.x += rotated.size.x * std::cos(rotated.heading);
        touching.position.y -= rotated.size.x * std::sin(rotated.heading);
        check(footprintOverlap(rotated, touching).iou == 0, "Rotated edge contact has zero area");
    }
    const auto corner = box(2, {2 - 1e-7, 2 - 1e-7, 0});
    for (auto method : {Assignment::Greedy, Assignment::Optimal}) {
        const auto weak = fusion.associateByOverlap({a}, {corner}, 0, method);
        check(weak.matches.size() == 1 && weak.matches.front().overlap > 0,
              "A small shared corner must pass threshold 0");
        const auto inside = box(2, {}, {1e-8, 2, 1e-8});
        check(fusion.associateByOverlap({a}, {inside}, 0, method).matches.size() == 1,
              "A contained footprint remains eligible when IoU is below machine epsilon");
    }
    check(fusion.associateByOverlap({a},{box(2,{1.8,0,0})},0).matches.size() == 1 &&
          fusion.associateByOverlap({a},{box(2,{1.8,0,0})},.1).matches.empty(),"Thresholds 0 and 0.1 distinguish weak intersections");
    check(fusion.associateByOverlap({a},{box(2,{0,0,10})},0).matches.empty(),"Separate heights cannot be one object");
    Track3D vehicle = a;
    vehicle.id = 2;
    vehicle.object_class = 1;
    result = fusion.associateByOverlap({a},{vehicle},0);
    check(result.matches.empty() && result.overlaps.size() == 1 && !result.overlaps[0].same_class,"Class mismatch remains visible in overlap report");
    result = fusion.associateByOverlap({box(1,{},{40,3,50})},{box(2,{20,0,0},{40,3,50})});
    check(result.matches.size() == 1 && result.matches[0].distance == 20,"Overlap mode does not use old 15m gate");
    const std::vector<Track3D> cv = {box(1,{0,0,0},{15,2,15}),box(2,{10,0,0},{15,2,15})};
    const std::vector<Track3D> adsb = {box(101,{1,0,0},{15,2,15}),box(102,{-10,0,0},{15,2,15})};
    check(fusion.associateByOverlap(cv,adsb).matches.size() == 1,"Greedy overlap ambiguity");
    check(fusion.associateByOverlap(cv,adsb,.1,Assignment::Optimal).matches.size() == 2,"Global overlap assignment recovers two pairs");
    bool rejected = false;
    try { fusion.associateByOverlap({a},{b},1.1); } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected,"IoU threshold outside range is rejected");
}

/** Сравнивает глобальное назначение IoU с независимым полным перебором. */
void exhaustiveOverlap() {
    std::mt19937 random(91);
    std::uniform_real_distribution<double> position(-6,6), angle(-3,3);
    SensorFusion fusion;
    for (int trial=0; trial<60; ++trial) {
        const int n = trial%4, m = (trial/4)%4;
        std::vector<Track3D> cv, adsb;
        for (int i=0; i<n; ++i) cv.push_back(box(i,{position(random),position(random),0},{6,2,8},angle(random)));
        for (int j=0; j<m; ++j) adsb.push_back(box(100+j,{position(random),position(random),0},{6,2,8},angle(random)));
        int best_count = -1;
        double best_iou = -1;
        std::function<void(int,unsigned,int,double)> enumerate = [&](int row, unsigned used, int count, double sum) {
            if (row == n) {
                if (count > best_count || (count == best_count && sum > best_iou)) {
                    best_count = count;
                    best_iou = sum;
                }
                return;
            }
            enumerate(row+1,used,count,sum);
            for (int column=0; column<m; ++column) {
                const double iou = footprintOverlap(cv[row],adsb[column]).iou;
                if (!(used&(1u<<column)) && iou >= .1) enumerate(row+1,used|(1u<<column),count+1,sum+iou);
            }
        };
        enumerate(0,0,0,0);
        const auto result = fusion.associateByOverlap(cv,adsb,.1,Assignment::Optimal);
        double sum = 0;
        for (const auto& pair : result.matches) sum += pair.overlap;
        check(static_cast<int>(result.matches.size()) == best_count,"Global IoU pair count agrees with enumeration");
        close(sum,best_iou,1e-9,"Global IoU sum agrees with enumeration");
    }
}

/** Проверяет обратные преобразования, направление камеры и метаданные интерполяции. */
void coordinatesAndMetadata() {
    for (const auto& point : {Geodetic{0,0,0},Geodetic{radians(43),radians(40),27},
                             Geodetic{radians(-89.5),radians(179),500},Geodetic{radians(90),0,100}}) {
        const Vec3 original = toEcef(point);
        check(squaredDistance(original,toEcef(fromEcef(original))) < 1e-10,"WGS84 inverse round trip");
    }
    const Attitude mounting = {radians(10),radians(5),radians(-3)};
    check(squaredDistance(bodyToCamera(cameraToBody({1,2,3},mounting),mounting),{1,2,3}) < 1e-20,"Camera rotation inverse");
    Track3D left = box(1,{},{2,2,4},radians(359));
    Track3D right = box(1,{2,4,6},{4,4,6},radians(1));
    right.timestamp = 2;
    Track3D middle;
    check(interpolateAdsb(1,left,right,3,middle),"Shape interpolation succeeds");
    close(std::sin(middle.heading),0,1e-12,"Heading follows short path across north");
    check(squaredDistance(middle.size,{3,3,5}) < 1e-20 && middle.object_class == 0,"Dimensions and class survive synchronization");
    right.object_class = 1;
    check(!interpolateAdsb(1,left,right,3,middle),"Interpolation cannot mix classes");

    const Camera camera = {{0,1.1,15.5},{},true};
    const auto scene = twoObjectScenario(camera);
    check(scene.frames.size() == 12 && scene.records == 24,"Two-object scenario frame count");
    const LocalFrame local(scene.frames.front().own);
    for (const auto& frame : scene.frames) {
        const CameraTransform transform(frame.own,frame.attitude,camera,local);
        for (const auto& obstacle : frame.obstacles) {
            const auto measured = cameraTrack(obstacle,frame.timestamp,transform);
            const auto truth = referenceTrack(obstacle,frame.timestamp,local);
            check(squaredDistance(measured.position,truth.position) < 1,"Known camera mounting recovers moving objects");
            close(measured.heading,truth.heading,1e-9,"Camera direction agrees with independent reference heading");
        }
    }
}

/** Проверяет наличие размеров, направления и признаков пересечения в отчёте. */
void reporting() {
    std::ostringstream output;
    const auto a = box(1), b = box(2,{1,0,0});
    const auto pair = SensorFusion().associateByOverlap({a},{b}).overlaps[0];
    printOverlapPair(output,2,a,b,pair);
    const std::string text = output.str();
    check(text.find("t=2") != std::string::npos && text.find("size=(") != std::string::npos &&
          text.find("heading_deg=") != std::string::npos && text.find("overlap_positive=1") != std::string::npos,
          "Terminal pair includes requested fields");
    output.str("");
    writeOverlapPair(output,2,a,b,pair,true);
    const auto row = output.str();
    check(std::count(row.begin(),row.end(),',') == 26,"Overlap CSV has 27 columns");
}
} // namespace

/** Запускает проверки площади пересечения и связанных преобразований. */
void testOverlap() {
    analyticalGeometry();
    invariants();
    matching();
    exhaustiveOverlap();
    coordinatesAndMetadata();
    reporting();
    std::cout << "Overlap geometry, assignment and report tests passed\n";
}
