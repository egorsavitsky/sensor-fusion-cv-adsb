#include "simulation.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

namespace {
void validate(const SimulationConfig& config) {
    if (!std::isfinite(config.period) || config.period <= 0.0 ||
        !std::isfinite(config.horizontal_sigma) || config.horizontal_sigma < 0.0 ||
        !std::isfinite(config.vertical_sigma) || config.vertical_sigma < 0.0 ||
        !std::isfinite(config.delay) || config.delay < 0.0 ||
        !std::isfinite(config.dropout) || config.dropout < 0.0 || config.dropout > 1.0 ||
        !std::isfinite(config.coverage) || config.coverage < 0.0 || config.coverage > 1.0) {
        throw std::invalid_argument("Invalid ADS-B simulation parameters");
    }
}
} // namespace

ModelData simulateAdsb(const std::vector<Frame>& frames, const LocalFrame& frame,
                       const SimulationConfig& config) {
    validate(config);
    std::map<int, std::vector<Track3D>> references;
    for (const auto& item : frames) {
        for (const auto& obstacle : item.obstacles) {
            references[obstacle.id].push_back({item.timestamp, obstacle.id,
                                              frame.fromGeodetic(obstacle.reference)});
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
            Vec3 position = truth[next].position;
            if (next + 1 < truth.size() && truth[next].timestamp != time) {
                const double ratio = (time-truth[next].timestamp) /
                                     (truth[next+1].timestamp-truth[next].timestamp);
                position = position + (truth[next+1].position-position)*ratio;
            }
            position.x += normal(noise_random) * config.horizontal_sigma;
            position.y += normal(noise_random) * config.horizontal_sigma;
            position.z += normal(noise_random) * config.vertical_sigma;
            if (uniform(loss_random) < config.dropout) continue;
            messages.push_back({{time, adsb_id, position}, time + config.delay});
            ++result.messages;
        }
    }
    return result;
}

AdsbSynchronizer::AdsbSynchronizer(const ModelData& model, double max_age,
                                   double max_gap, bool offline)
    : model_(model), max_age_(max_age), max_gap_(max_gap), offline_(offline) {
    if (!std::isfinite(max_age) || max_age < 0.0 || !std::isfinite(max_gap) || max_gap <= 0.0) {
        throw std::invalid_argument("Invalid synchronization limits");
    }
}

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
