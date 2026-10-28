#pragma once

#include "fusion.hpp"

struct Geodetic {
    double latitude = 0.0; // радианы
    double longitude = 0.0;
    double altitude = 0.0; // метры над эллипсоидом WGS84
};

struct Attitude {
    double yaw = 0.0; // радианы; курс от севера по часовой стрелке
    double pitch = 0.0;
    double roll = 0.0;
};

struct Camera {
    Vec3 offset; // в осях ВС: вправо, вверх, вперёд
    Attitude mounting;
    bool position_is_center = true;
};

double radians(double degrees);
Vec3 toEcef(const Geodetic& point);
Vec3 enuToEcef(const Vec3& vector, const Geodetic& origin);
Vec3 ecefToEnu(const Vec3& vector, const Geodetic& origin);
Vec3 bodyToEnu(const Vec3& vector, const Attitude& attitude);
Vec3 enuToBody(const Vec3& vector, const Attitude& attitude);
Vec3 cameraToBody(const Vec3& vector, const Attitude& mounting);

class LocalFrame {
public:
    explicit LocalFrame(const Geodetic& origin);
    Vec3 fromGeodetic(const Geodetic& point) const;
    Vec3 fromEcef(const Vec3& point) const;
    Vec3 vectorFromEcef(const Vec3& vector) const;

private:
    Vec3 ecef_origin_;
    Vec3 east_, north_, up_;
};

class CameraTransform {
public:
    CameraTransform(const Geodetic& own, const Attitude& attitude,
                    const Camera& camera, const LocalFrame& frame);
    Vec3 apply(const Vec3& position, double object_height) const;

private:
    Vec3 origin_, x_, y_, z_, height_correction_;
};

Vec3 cameraPosition(const Vec3& position, double object_height,
                    const Geodetic& own, const Attitude& attitude,
                    const Camera& camera, const LocalFrame& frame);
Vec3 cameraOffset(const Vec3& position, double object_height,
                  const Geodetic& obstacle, const Geodetic& own,
                  const Attitude& attitude, const Camera& camera);
