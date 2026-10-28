#include "geometry.hpp"

#include <cmath>
#include <stdexcept>

namespace {
double dot(const Vec3& a, const Vec3& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

struct Axes { Vec3 right, up, forward; };

Axes bodyAxes(const Attitude& attitude) {
    const double sy = std::sin(attitude.yaw), cy = std::cos(attitude.yaw);
    const double sp = std::sin(attitude.pitch), cp = std::cos(attitude.pitch);
    const double sr = std::sin(attitude.roll), cr = std::cos(attitude.roll);
    const Vec3 right = {cy, -sy, 0.0};
    const Vec3 up = {-sy*sp, -cy*sp, cp};
    return {right*cr - up*sr, right*sr + up*cr, {sy*cp, cy*cp, sp}};
}
} // namespace

double radians(double degrees) { return degrees * 3.14159265358979323846 / 180.0; }

Vec3 toEcef(const Geodetic& point) {
    if (!std::isfinite(point.latitude) || !std::isfinite(point.longitude) ||
        !std::isfinite(point.altitude) || std::abs(point.latitude) > radians(90) ||
        std::abs(point.longitude) > radians(180)) {
        throw std::invalid_argument("Invalid WGS84 coordinates (expected radians)");
    }
    const double a = 6378137.0, e2 = 6.6943799901413165e-3;
    const double s = std::sin(point.latitude), c = std::cos(point.latitude);
    const double n = a / std::sqrt(1.0 - e2*s*s);
    return {(n + point.altitude)*c*std::cos(point.longitude),
            (n + point.altitude)*c*std::sin(point.longitude),
            (n*(1.0-e2) + point.altitude)*s};
}

Vec3 enuToEcef(const Vec3& vector, const Geodetic& origin) {
    const double sl = std::sin(origin.longitude), cl = std::cos(origin.longitude);
    const double sp = std::sin(origin.latitude), cp = std::cos(origin.latitude);
    return {-sl*vector.x-sp*cl*vector.y+cp*cl*vector.z,
             cl*vector.x-sp*sl*vector.y+cp*sl*vector.z,
             cp*vector.y+sp*vector.z};
}

Vec3 ecefToEnu(const Vec3& vector, const Geodetic& origin) {
    const double sl = std::sin(origin.longitude), cl = std::cos(origin.longitude);
    const double sp = std::sin(origin.latitude), cp = std::cos(origin.latitude);
    return {-sl*vector.x+cl*vector.y,
            -sp*cl*vector.x-sp*sl*vector.y+cp*vector.z,
             cp*cl*vector.x+cp*sl*vector.y+sp*vector.z};
}

Vec3 bodyToEnu(const Vec3& vector, const Attitude& attitude) {
    const Axes axes = bodyAxes(attitude);
    return axes.right*vector.x + axes.up*vector.y + axes.forward*vector.z;
}

Vec3 enuToBody(const Vec3& vector, const Attitude& attitude) {
    const Axes axes = bodyAxes(attitude);
    return {dot(vector, axes.right), dot(vector, axes.up), dot(vector, axes.forward)};
}

Vec3 cameraToBody(const Vec3& vector, const Attitude& mounting) {
    const Vec3 enu = bodyToEnu(vector, mounting);
    return {enu.x, enu.z, enu.y};
}

LocalFrame::LocalFrame(const Geodetic& origin) : ecef_origin_(toEcef(origin)),
    east_(enuToEcef({1,0,0},origin)), north_(enuToEcef({0,1,0},origin)),
    up_(enuToEcef({0,0,1},origin)) {}

Vec3 LocalFrame::fromGeodetic(const Geodetic& point) const {
    return fromEcef(toEcef(point));
}

Vec3 LocalFrame::fromEcef(const Vec3& point) const {
    return vectorFromEcef(point - ecef_origin_);
}

Vec3 LocalFrame::vectorFromEcef(const Vec3& vector) const {
    return {dot(vector,east_), dot(vector,north_), dot(vector,up_)};
}

CameraTransform::CameraTransform(const Geodetic& own, const Attitude& attitude,
                                  const Camera& camera, const LocalFrame& frame) {
    auto rotate = [&](const Vec3& body) {
        return frame.vectorFromEcef(enuToEcef(bodyToEnu(body,attitude),own));
    };
    origin_ = frame.fromGeodetic(own) + rotate(camera.offset);
    x_ = rotate(cameraToBody({1,0,0},camera.mounting));
    y_ = rotate(cameraToBody({0,1,0},camera.mounting));
    z_ = rotate(cameraToBody({0,0,1},camera.mounting));
    if (camera.position_is_center) {
        height_correction_ = frame.vectorFromEcef(enuToEcef({0,0,0.5},own));
    }
}

Vec3 CameraTransform::apply(const Vec3& position, double object_height) const {
    return origin_ + x_*position.x + y_*position.y + z_*position.z - height_correction_*object_height;
}

Vec3 cameraPosition(const Vec3& position, double object_height,
                    const Geodetic& own, const Attitude& attitude,
                    const Camera& camera, const LocalFrame& frame) {
    return CameraTransform(own,attitude,camera,frame).apply(position,object_height);
}

Vec3 cameraOffset(const Vec3& position, double object_height,
                  const Geodetic& obstacle, const Geodetic& own,
                  const Attitude& attitude, const Camera& camera) {
    Vec3 relative = ecefToEnu(toEcef(obstacle) - toEcef(own), own);
    if (camera.position_is_center) relative.z += object_height / 2.0;
    return enuToBody(relative, attitude) - cameraToBody(position, camera.mounting);
}
