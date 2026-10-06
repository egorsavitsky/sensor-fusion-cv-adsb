#include <iostream>
#include <windows.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

struct Track3D {
    double timestamp;
    int id;
    double x, y, z;
};

struct FusionResult {
    int cv_id;
    int adsb_id;
    double probability; // От 0.0 до 1.0
};

class SensorFusion {
private:
    const double MAX_DISTANCE = 15.0; // порог погрешности

public:
    // интерполяция
    Track3D interpolateAdsb(double target_t, const Track3D& t0, const Track3D& t1) {
        if (t1.timestamp == t0.timestamp) return t0;

        double ratio = (target_t - t0.timestamp) / (t1.timestamp - t0.timestamp);
        
        Track3D result;
        result.timestamp = target_t;
        result.id = t0.id;
        result.x = t0.x + ratio * (t1.x - t0.x);
        result.y = t0.y + ratio * (t1.y - t0.y);
        result.z = t0.z + ratio * (t1.z - t0.z);
        
        return result;
    }

    // расчет вероятности
    double calculateProbability(const Track3D& cv_track, const Track3D& adsb_track) {
        double distance = std::sqrt(
            std::pow(cv_track.x - adsb_track.x, 2) +
            std::pow(cv_track.y - adsb_track.y, 2) +
            std::pow(cv_track.z - adsb_track.z, 2)
        );

        if (distance >= MAX_DISTANCE) {
            return 0.0;
        }
        
        return 1.0 - (distance / MAX_DISTANCE);
    }

    // ассоциация треков жадным алгоритмом
    std::vector<FusionResult> associateTracks(const std::vector<Track3D>& cv_tracks, const std::vector<Track3D>& adsb_tracks) {
        std::vector<FusionResult> results;
        
        struct Pair {
            int cv_idx;
            int adsb_idx;
            double prob;
        };
        std::vector<Pair> all_pairs;

        // пострение матрицы вероятностей
        for (size_t i = 0; i < cv_tracks.size(); ++i) {
            for (size_t j = 0; j < adsb_tracks.size(); ++j) {
                double prob = calculateProbability(cv_tracks[i], adsb_tracks[j]);
                if (prob > 0.0) {
                    all_pairs.push_back({(int)i, (int)j, prob});
                }
            }
        }

        std::sort(all_pairs.begin(), all_pairs.end(), [](const Pair& a, const Pair& b) {
            return a.prob > b.prob;
        });

        std::vector<bool> cv_assigned(cv_tracks.size(), false);
        std::vector<bool> adsb_assigned(adsb_tracks.size(), false);

        // идем по самым вероятным парам и связываем свободные объекты
        for (const auto& pair : all_pairs) {
            if (!cv_assigned[pair.cv_idx] && !adsb_assigned[pair.adsb_idx]) {
                results.push_back({
                    cv_tracks[pair.cv_idx].id, 
                    adsb_tracks[pair.adsb_idx].id, 
                    pair.prob
                });
                cv_assigned[pair.cv_idx] = true;
                adsb_assigned[pair.adsb_idx] = true;
            }
        }

        return results;
    }
};


int main() {
    SetConsoleOutputCP(CP_UTF8);
    SensorFusion fusion;

    std::cout << "Ассоциация треков данных от АЗН-В и СТЗ:\n\n";

    // тестовые данные от AЗН-В
    Track3D adsb_101_t1 = {1.0, 101, 100.0, 50.0, 0.0};
    Track3D adsb_102_t1 = {1.0, 102, -50.0, 20.0, 0.0};
    
    // Те же два самолета в момент времени t = 2.0
    Track3D adsb_101_t2 = {2.0, 101, 110.0, 50.0, 0.0};
    Track3D adsb_102_t2 = {2.0, 102, -50.0, 30.0, 0.0};

    // тестовые данные от СТЗ
    double current_time = 1.5;
    std::vector<Track3D> cv_tracks = {
        {current_time, 1, 104.5, 49.5, 0.0},
        {current_time, 2, -49.0, 24.0, 0.0},
        {current_time, 3, 500.0, 500.0, 0.0}
    };

    std::vector<Track3D> interpolated_adsb;
    interpolated_adsb.push_back(fusion.interpolateAdsb(current_time, adsb_101_t1, adsb_101_t2));
    interpolated_adsb.push_back(fusion.interpolateAdsb(current_time, adsb_102_t1, adsb_102_t2));

    std::cout << "[Ожидаемые позиции АЗН-В на t=" << current_time << "]\n";
    for (const auto& track : interpolated_adsb) {
        std::cout << "ADSB_ID: " << track.id 
                  << " | X: " << track.x << ", Y: " << track.y << ", Z: " << track.z << "\n";
    }

    std::vector<FusionResult> results = fusion.associateTracks(cv_tracks, interpolated_adsb);

    std::cout << "\n[Результаты слияния]\n";
    for (const auto& res : results) {
        std::cout << "CV_Track_" << res.cv_id << " <== ассоциирован с ==> ADSB_Track_" << res.adsb_id 
                  << " | Вероятность: " << std::fixed << std::setprecision(1) << (res.probability * 100.0) << "%\n";
    }

    return 0;
}