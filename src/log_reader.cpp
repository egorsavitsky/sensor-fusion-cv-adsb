#include "log_reader.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

namespace {
std::vector<std::string> split(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    for (;;) {
        const std::size_t end = line.find('\t', start);
        fields.push_back(line.substr(start, end - start));
        if (end == std::string::npos) break;
        start = end + 1;
    }
    return fields;
}

double number(const std::string& field) {
    std::size_t consumed = 0;
    const double value = std::stod(field, &consumed);
    if (consumed != field.size() || !std::isfinite(value)) {
        throw std::runtime_error("Invalid number: " + field);
    }
    return value;
}

int integer(const std::string& field) {
    std::size_t consumed = 0;
    const long long value = std::stoll(field, &consumed);
    if (consumed != field.size() || value < std::numeric_limits<int>::min() ||
        value > std::numeric_limits<int>::max()) {
        throw std::runtime_error("Invalid integer: " + field);
    }
    return static_cast<int>(value);
}

double median(std::vector<double>& values) {
    if (values.empty()) throw std::runtime_error("No obstacles in calibration frames");
    const std::size_t middle = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    const double high = values[middle];
    if (values.size() % 2 != 0) return high;
    return (*std::max_element(values.begin(), values.begin() + middle) + high) / 2.0;
}
} // namespace

LogData readLog(std::istream& input) {
    std::string line;
    if (!std::getline(input, line)) throw std::runtime_error("Empty log file");
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.compare(0, 3, "\xef\xbb\xbf") == 0) line.erase(0, 3);
    const auto header = split(line);
    std::map<std::string, std::size_t> columns;
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (header[i].empty() || !columns.emplace(header[i], i).second) {
            throw std::runtime_error("Empty or duplicate column name");
        }
    }
    const std::vector<std::string> required = {
        "ID", "Class", "Size.Y", "Position.X", "Position.Y", "Position.Z", "Timestamp(Ms)",
        "Obs.Lat(rad)", "Obs.Lon(rad)", "Obs.Alt(M)",
        "Own.Lat(rad)", "Own.Lon(rad)", "Own.Alt(M)",
        "Own.Pitch(Deg)", "Own.Roll(Deg)", "Own.Yaw(Deg)"};
    for (const auto& name : required) {
        if (columns.count(name) == 0) throw std::runtime_error("Missing column: " + name);
    }
    LogData result;
    double origin_ms = 0.0, last_ms = -1.0;
    std::set<int> frame_ids;
    bool frame_empty = false;
    std::size_t line_number = 1;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        try {
            if (line.empty()) throw std::runtime_error("Blank record");
            const auto fields = split(line);
            auto field = [&](const std::string& name) -> const std::string& {
                const std::size_t index = columns.at(name);
                if (index >= fields.size()) throw std::runtime_error("Missing field: " + name);
                return fields[index];
            };
            const bool empty = field("ID").empty();
            if (fields.size() > header.size() || (!empty && fields.size() != header.size())) {
                throw std::runtime_error("Wrong number of fields");
            }
            // В исходном логе пустые кадры не содержат двух последних полей скорости.
            if (fields.size() < header.size()) {
                for (std::size_t i = fields.size(); i < header.size(); ++i) {
                    if (header[i] != "Own.GroundSpeed(km/h)" && header[i] != "Own.RoC(m/s)") {
                        throw std::runtime_error("Truncated record");
                    }
                }
                ++result.short_empty_records;
            }
            const double time_ms = number(field("Timestamp(Ms)"));
            if (time_ms < 0.0 || time_ms < last_ms) throw std::runtime_error("Time is not ordered");
            const Geodetic own = {number(field("Own.Lat(rad)")), number(field("Own.Lon(rad)")),
                                  number(field("Own.Alt(M)"))};
            toEcef(own);
            const Attitude attitude = {radians(number(field("Own.Yaw(Deg)"))),
                                       radians(number(field("Own.Pitch(Deg)"))),
                                       radians(number(field("Own.Roll(Deg)")))};
            if (result.frames.empty() || time_ms != last_ms) {
                if (result.frames.empty()) origin_ms = time_ms;
                Frame frame;
                frame.timestamp = (time_ms - origin_ms) / 1000.0;
                frame.own = own;
                frame.attitude = attitude;
                result.frames.push_back(frame);
                frame_ids.clear();
                frame_empty = false;
            } else {
                const Frame& frame = result.frames.back();
                if (squaredDistance(toEcef(own), toEcef(frame.own)) > 1e-8 ||
                    std::abs(attitude.yaw-frame.attitude.yaw) > 1e-10 ||
                    std::abs(attitude.pitch-frame.attitude.pitch) > 1e-10 ||
                    std::abs(attitude.roll-frame.attitude.roll) > 1e-10) {
                    throw std::runtime_error("Inconsistent ownship data within frame");
                }
            }
            last_ms = time_ms;
            ++result.records;
            if (empty) {
                if (!result.frames.back().obstacles.empty() || frame_empty) {
                    throw std::runtime_error("Empty record mixed with obstacles");
                }
                frame_empty = true;
                ++result.empty_records;
                continue;
            }
            if (frame_empty) throw std::runtime_error("Obstacle mixed with empty frame record");
            Obstacle obstacle;
            obstacle.id = integer(field("ID"));
            obstacle.object_class = integer(field("Class"));
            obstacle.position = {number(field("Position.X")), number(field("Position.Y")),
                                 number(field("Position.Z"))};
            obstacle.height = number(field("Size.Y"));
            if (obstacle.height < 0.0) throw std::runtime_error("Negative object height");
            obstacle.reference = {number(field("Obs.Lat(rad)")), number(field("Obs.Lon(rad)")),
                                  number(field("Obs.Alt(M)"))};
            toEcef(obstacle.reference);
            if (!frame_ids.insert(obstacle.id).second) throw std::runtime_error("Duplicate obstacle ID");
            result.frames.back().obstacles.push_back(obstacle);
        } catch (const std::exception& error) {
            throw std::runtime_error("Log line " + std::to_string(line_number) + ": " + error.what());
        }
    }
    if (input.bad()) throw std::runtime_error("Error reading log file");
    if (result.frames.empty()) throw std::runtime_error("Log contains no frames");
    return result;
}

LogData readLog(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open log: " + path);
    return readLog(input);
}

Camera calibrateCamera(const std::vector<Frame>& frames, std::size_t count, Camera camera) {
    if (count == 0) return camera;
    if (count >= frames.size()) {
        throw std::invalid_argument("Calibration must leave at least one evaluation frame");
    }
    std::vector<double> x, y, z;
    for (std::size_t i = 0; i < count; ++i) {
        for (const auto& obstacle : frames[i].obstacles) {
            const Vec3 offset = cameraOffset(obstacle.position, obstacle.height, obstacle.reference,
                                              frames[i].own, frames[i].attitude, camera);
            x.push_back(offset.x);
            y.push_back(offset.y);
            z.push_back(offset.z);
        }
    }
    camera.offset = {median(x), median(y), median(z)};
    return camera;
}
