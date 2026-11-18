#pragma once

#include "geometry.hpp"

#include <cstddef>
#include <istream>
#include <string>
#include <vector>

struct Obstacle {
    int id = 0;
    int object_class = 0;
    Vec3 position;
    Vec3 size;
    Geodetic reference;
    Vec3 direction = {0,0,1}; // направление в осях камеры
    double reference_heading = 0.0; // курс объекта из Obs.Yaw, радианы
};

struct Frame {
    double timestamp = 0.0;
    Geodetic own;
    Attitude attitude;
    std::vector<Obstacle> obstacles;
};

struct LogData {
    std::vector<Frame> frames;
    std::size_t records = 0;
    std::size_t empty_records = 0;
    std::size_t short_empty_records = 0;
};

LogData readLog(std::istream& input);
LogData readLog(const std::string& path);
Camera calibrateCamera(const std::vector<Frame>& frames, std::size_t count, Camera camera);
