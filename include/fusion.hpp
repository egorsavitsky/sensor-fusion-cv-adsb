#pragma once

#include <vector>

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

Vec3 operator+(const Vec3& a, const Vec3& b);
Vec3 operator-(const Vec3& a, const Vec3& b);
Vec3 operator*(const Vec3& a, double scale);
double squaredDistance(const Vec3& a, const Vec3& b);

struct Track3D {
    double timestamp = 0.0; // секунды от начала записи
    int id = 0;
    Vec3 position;
};

struct Match {
    int cv_id;
    int adsb_id;
    double distance;
    double score; // оценка близости, не статистическая вероятность
};

struct Association {
    std::vector<Match> matches;
    std::vector<int> cv_only;
    std::vector<int> adsb_only;
};

enum class Assignment { Greedy, Optimal };

bool interpolateAdsb(double time, const Track3D& left, const Track3D& right,
                     double max_gap, Track3D& result);

class SensorFusion {
public:
    explicit SensorFusion(double max_distance = 15.0);
    Association associateTracks(const std::vector<Track3D>& cv,
                                const std::vector<Track3D>& adsb,
                                Assignment method = Assignment::Greedy) const;
    Association associateReference(const std::vector<Track3D>& cv,
                                   const std::vector<Track3D>& adsb) const;

private:
    double max_distance_;
};
