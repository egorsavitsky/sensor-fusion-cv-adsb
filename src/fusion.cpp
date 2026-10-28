#include "fusion.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

Vec3 operator+(const Vec3& a, const Vec3& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 operator-(const Vec3& a, const Vec3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 operator*(const Vec3& a, double scale) {
    return {a.x * scale, a.y * scale, a.z * scale};
}

double squaredDistance(const Vec3& a, const Vec3& b) {
    const Vec3 d = a - b;
    return d.x * d.x + d.y * d.y + d.z * d.z;
}

namespace {
struct Candidate {
    std::size_t cv;
    std::size_t adsb;
    double cost;
};

void checkTracks(const std::vector<Track3D>& tracks) {
    std::vector<int> ids;
    ids.reserve(tracks.size());
    for (const auto& track : tracks) {
        const Vec3& p = track.position;
        if (!std::isfinite(track.timestamp) || !std::isfinite(p.x) ||
            !std::isfinite(p.y) || !std::isfinite(p.z)) {
            throw std::invalid_argument("Non-finite track coordinate or time");
        }
        ids.push_back(track.id);
        if (std::abs(track.timestamp - tracks.front().timestamp) > 1e-6) {
            throw std::invalid_argument("Tracks must refer to the same frame time");
        }
    }
    std::sort(ids.begin(), ids.end());
    if (std::adjacent_find(ids.begin(), ids.end()) != ids.end()) {
        throw std::invalid_argument("Duplicate track ID in one frame");
    }
}

void checkTimes(const std::vector<Track3D>& cv, const std::vector<Track3D>& adsb) {
    if (!cv.empty() && !adsb.empty() && std::abs(cv.front().timestamp-adsb.front().timestamp) > 1e-6) {
        throw std::invalid_argument("Synchronize sources before association");
    }
}

Association collect(const std::vector<Track3D>& cv, const std::vector<Track3D>& adsb,
                    const std::vector<Candidate>& selected, double limit) {
    Association result;
    std::vector<bool> used_cv(cv.size(), false), used_adsb(adsb.size(), false);
    for (const auto& pair : selected) {
        const double distance = std::sqrt(squaredDistance(cv[pair.cv].position,
                                                          adsb[pair.adsb].position));
        result.matches.push_back({cv[pair.cv].id, adsb[pair.adsb].id,
                                  distance, 1.0 - distance / limit});
        used_cv[pair.cv] = true;
        used_adsb[pair.adsb] = true;
    }
    for (std::size_t i = 0; i < cv.size(); ++i) {
        if (!used_cv[i]) result.cv_only.push_back(cv[i].id);
    }
    for (std::size_t i = 0; i < adsb.size(); ++i) {
        if (!used_adsb[i]) result.adsb_only.push_back(adsb[i].id);
    }
    return result;
}

std::vector<Candidate> greedy(std::vector<Candidate>& pairs, std::size_t n, std::size_t m) {
    std::sort(pairs.begin(), pairs.end(), [](const Candidate& a, const Candidate& b) {
        if (a.cost != b.cost) return a.cost < b.cost;
        if (a.cv != b.cv) return a.cv < b.cv;
        return a.adsb < b.adsb;
    });
    std::vector<bool> used_cv(n, false), used_adsb(m, false);
    std::vector<Candidate> selected;
    selected.reserve(std::min(n, m));
    for (const auto& pair : pairs) {
        if (!used_cv[pair.cv] && !used_adsb[pair.adsb]) {
            selected.push_back(pair);
            used_cv[pair.cv] = true;
            used_adsb[pair.adsb] = true;
        }
    }
    return selected;
}

// Венгерский алгоритм. Дополнительные столбцы позволяют оставить СТЗ без пары.
std::vector<Candidate> optimal(const std::vector<Track3D>& cv,
                               const std::vector<Track3D>& adsb, double limit_sq) {
    const std::size_t n = cv.size(), m = adsb.size() + n;
    if (n == 0) return {};
    const double unmatched = static_cast<double>(n + 1);
    const double forbidden = unmatched * 2.0;
    std::vector<std::vector<double>> cost(n, std::vector<double>(m, unmatched));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < adsb.size(); ++j) {
            const double d2 = squaredDistance(cv[i].position, adsb[j].position);
            cost[i][j] = d2 < limit_sq ? d2 / limit_sq : forbidden;
        }
    }
    std::vector<double> u(n + 1), v(m + 1);
    std::vector<std::size_t> p(m + 1), way(m + 1);
    for (std::size_t i = 1; i <= n; ++i) {
        p[0] = i;
        std::size_t column = 0;
        std::vector<double> best(m + 1, std::numeric_limits<double>::infinity());
        std::vector<bool> used(m + 1, false);
        do {
            used[column] = true;
            const std::size_t row = p[column];
            double delta = std::numeric_limits<double>::infinity();
            std::size_t next = 0;
            for (std::size_t j = 1; j <= m; ++j) {
                if (used[j]) continue;
                const double value = cost[row - 1][j - 1] - u[row] - v[j];
                if (value < best[j]) {
                    best[j] = value;
                    way[j] = column;
                }
                if (best[j] < delta) {
                    delta = best[j];
                    next = j;
                }
            }
            for (std::size_t j = 0; j <= m; ++j) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    best[j] -= delta;
                }
            }
            column = next;
        } while (p[column] != 0);
        do {
            const std::size_t previous = way[column];
            p[column] = p[previous];
            column = previous;
        } while (column != 0);
    }
    std::vector<Candidate> selected;
    for (std::size_t j = 1; j <= adsb.size(); ++j) {
        if (p[j] != 0 && cost[p[j] - 1][j - 1] < 1.0) {
            selected.push_back({p[j] - 1, j - 1, cost[p[j] - 1][j - 1]});
        }
    }
    return selected;
}
} // namespace

bool interpolateAdsb(double time, const Track3D& left, const Track3D& right,
                     double max_gap, Track3D& result) {
    if (!std::isfinite(time) || !std::isfinite(max_gap) || max_gap <= 0.0 ||
        left.id != right.id || !std::isfinite(left.timestamp) ||
        !std::isfinite(right.timestamp) || right.timestamp < left.timestamp ||
        time < left.timestamp || time > right.timestamp ||
        right.timestamp - left.timestamp > max_gap) return false;
    if (!std::isfinite(squaredDistance(left.position, right.position))) return false;
    if (left.timestamp == right.timestamp) {
        if (squaredDistance(left.position, right.position) != 0.0) return false;
        result = left;
        return true;
    }
    const double ratio = (time - left.timestamp) / (right.timestamp - left.timestamp);
    result = {time, left.id, left.position + (right.position - left.position) * ratio};
    return true;
}

SensorFusion::SensorFusion(double max_distance) : max_distance_(max_distance) {
    if (!std::isfinite(max_distance) || max_distance <= 0.0 ||
        !std::isfinite(max_distance * max_distance) || max_distance * max_distance == 0.0) {
        throw std::invalid_argument("Distance threshold must be finite and positive");
    }
}

Association SensorFusion::associateTracks(const std::vector<Track3D>& cv,
                                          const std::vector<Track3D>& adsb,
                                          Assignment method) const {
    checkTracks(cv);
    checkTracks(adsb);
    checkTimes(cv, adsb);
    const double limit_sq = max_distance_ * max_distance_;
    if (method == Assignment::Optimal) {
        return collect(cv, adsb, optimal(cv, adsb, limit_sq), max_distance_);
    }
    std::vector<Candidate> pairs;
    for (std::size_t i = 0; i < cv.size(); ++i) {
        for (std::size_t j = 0; j < adsb.size(); ++j) {
            const double dx = cv[i].position.x - adsb[j].position.x;
            if (std::abs(dx) >= max_distance_) continue;
            const double dy = cv[i].position.y - adsb[j].position.y;
            if (std::abs(dy) >= max_distance_) continue;
            const double dz = cv[i].position.z - adsb[j].position.z;
            if (std::abs(dz) >= max_distance_) continue;
            const double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 < limit_sq) pairs.push_back({i, j, d2});
        }
    }
    return collect(cv, adsb, greedy(pairs, cv.size(), adsb.size()), max_distance_);
}

Association SensorFusion::associateReference(const std::vector<Track3D>& cv,
                                             const std::vector<Track3D>& adsb) const {
    checkTracks(cv);
    checkTracks(adsb);
    checkTimes(cv, adsb);
    std::vector<Candidate> pairs;
    for (std::size_t i = 0; i < cv.size(); ++i) {
        for (std::size_t j = 0; j < adsb.size(); ++j) {
            const Vec3 d = cv[i].position - adsb[j].position;
            const double distance = std::sqrt(std::pow(d.x, 2) + std::pow(d.y, 2) +
                                              std::pow(d.z, 2));
            if (distance < max_distance_) {
                pairs.push_back({i, j, distance});
            }
        }
    }
    return collect(cv, adsb, greedy(pairs, cv.size(), adsb.size()), max_distance_);
}
