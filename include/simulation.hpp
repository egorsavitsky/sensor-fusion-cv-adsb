#pragma once

#include "log_reader.hpp"

#include <cstdint>
#include <map>
#include <vector>

struct SimulationConfig {
    double period = 1.0;
    double horizontal_sigma = 3.0;
    double vertical_sigma = 2.0;
    double delay = 0.2;
    double dropout = 0.1;
    double coverage = 0.8;
    std::uint32_t seed = 42;
};

struct AdsbMessage {
    Track3D track;
    double received_at;
};

struct ModelData {
    std::map<int, std::vector<AdsbMessage>> streams;
    std::map<int, int> reference_ids; // только для оценки, в SensorFusion не передаётся
    std::size_t messages = 0;
};

ModelData simulateAdsb(const std::vector<Frame>& frames, const LocalFrame& frame,
                       const SimulationConfig& config);

class AdsbSynchronizer {
public:
    AdsbSynchronizer(const ModelData& model, double max_age, double max_gap, bool offline);
    std::vector<Track3D> at(double time);

private:
    const ModelData& model_;
    std::map<int, std::size_t> cursors_;
    double max_age_;
    double max_gap_;
    bool offline_;
    double last_time_ = -1.0;
};
