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
    double height = 0.0;
    Geodetic reference;
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
