#include "fusion.hpp"
#include "geometry.hpp"
#include "log_reader.hpp"
#include "simulation.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _WIN32
#include <windows.h>
#endif

namespace {
using Clock = std::chrono::steady_clock;

struct Options {
    std::string input, summary, matches;
    std::size_t calibration_frames = 1000;
    std::size_t limit_frames = 0;
    double max_distance = 15.0;
    double max_age = 2.5;
    double max_gap = 2.5;
    bool offline = false;
    bool benchmark = false;
    Assignment assignment = Assignment::Greedy;
    Camera camera;
    SimulationConfig simulation;
};

double numeric(const std::string& text) {
    std::size_t used = 0;
    const double value = std::stod(text, &used);
    if (used != text.size() || !std::isfinite(value)) throw std::invalid_argument("Invalid number: " + text);
    return value;
}

std::size_t count(const std::string& text) {
    const double value = numeric(text);
    if (value < 0.0 || value != std::floor(value) || value > 1000000000.0) {
        throw std::invalid_argument("Expected a nonnegative integer: " + text);
    }
    return static_cast<std::size_t>(value);
}

void help() {
    std::cout <<
        "Usage: fusion_core --input log.txt [options]\n"
        "Without arguments: small 3D example.\n\n"
        "  --calibration-frames N  Estimate camera offset on the first N frames (1000)\n"
        "  --limit-frames N        Evaluate only the first N frames; 0 means all\n"
        "  --camera-right M        Manual offset to the right (0)\n"
        "  --camera-up M           Manual upward offset (0)\n"
        "  --camera-forward M      Manual forward offset (0)\n"
        "  --camera-yaw DEG        Fixed camera yaw relative to ownship (0)\n"
        "  --camera-pitch DEG      Fixed camera pitch (0)\n"
        "  --camera-roll DEG       Fixed camera roll (0)\n"
        "  --position-anchor MODE  center or base; center subtracts Size.Y/2 (center)\n"
        "  --assignment MODE       greedy or optimal (greedy)\n"
        "  --max-distance M        Strict association distance limit (15)\n"
        "  --max-age SEC           Maximum age of a held measurement (2.5)\n"
        "  --max-gap SEC           Maximum interpolation interval (2.5)\n"
        "  --offline               Use future samples for offline interpolation\n"
        "  --period SEC            Synthetic ADS-B message period (1)\n"
        "  --horizontal-sigma M    Synthetic horizontal noise per axis (3)\n"
        "  --vertical-sigma M      Synthetic vertical noise (2)\n"
        "  --delay SEC             Synthetic delivery delay (0.2)\n"
        "  --dropout FRACTION      Synthetic message loss, 0..1 (0.1)\n"
        "  --coverage FRACTION     Fraction of objects with synthetic ADS-B (0.8)\n"
        "  --seed N                Random seed (42)\n"
        "  --benchmark             Compare reference and optimized greedy code\n"
        "  --summary FILE          Write settings and metrics as CSV\n"
        "  --matches FILE          Write matches and unmatched objects as CSV\n";
}

Options parse(int argc, char** argv) {
    Options options;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        if (key == "--offline") { options.offline = true; continue; }
        if (key == "--benchmark") { options.benchmark = true; continue; }
        if (i + 1 >= argc) throw std::invalid_argument("Missing value for " + key);
        const std::string value = argv[++i];
        if (key == "--input") options.input = value;
        else if (key == "--summary") options.summary = value;
        else if (key == "--matches") options.matches = value;
        else if (key == "--calibration-frames") options.calibration_frames = count(value);
        else if (key == "--limit-frames") options.limit_frames = count(value);
        else if (key == "--seed") options.simulation.seed = static_cast<std::uint32_t>(count(value));
        else if (key == "--max-distance") options.max_distance = numeric(value);
        else if (key == "--max-age") options.max_age = numeric(value);
        else if (key == "--max-gap") options.max_gap = numeric(value);
        else if (key == "--period") options.simulation.period = numeric(value);
        else if (key == "--horizontal-sigma") options.simulation.horizontal_sigma = numeric(value);
        else if (key == "--vertical-sigma") options.simulation.vertical_sigma = numeric(value);
        else if (key == "--delay") options.simulation.delay = numeric(value);
        else if (key == "--dropout") options.simulation.dropout = numeric(value);
        else if (key == "--coverage") options.simulation.coverage = numeric(value);
        else if (key == "--camera-right") options.camera.offset.x = numeric(value);
        else if (key == "--camera-up") options.camera.offset.y = numeric(value);
        else if (key == "--camera-forward") options.camera.offset.z = numeric(value);
        else if (key == "--camera-yaw") options.camera.mounting.yaw = radians(numeric(value));
        else if (key == "--camera-pitch") options.camera.mounting.pitch = radians(numeric(value));
        else if (key == "--camera-roll") options.camera.mounting.roll = radians(numeric(value));
        else if (key == "--position-anchor" && value == "center") options.camera.position_is_center = true;
        else if (key == "--position-anchor" && value == "base") options.camera.position_is_center = false;
        else if (key == "--assignment" && value == "greedy") options.assignment = Assignment::Greedy;
        else if (key == "--assignment" && value == "optimal") options.assignment = Assignment::Optimal;
        else throw std::invalid_argument("Unknown option or value: " + key + " " + value);
    }
    if (options.input.empty()) throw std::invalid_argument("Specify --input log.txt");
    if ((!options.summary.empty() && options.summary == options.input) ||
        (!options.matches.empty() && options.matches == options.input) ||
        (!options.summary.empty() && options.summary == options.matches)) {
        throw std::invalid_argument("Input and output paths must be different");
    }
    for (const auto& output : {options.summary, options.matches}) {
        if (!output.empty() && std::ifstream(output).good()) {
            throw std::invalid_argument("Output already exists; choose a new filename: " + output);
        }
    }
    return options;
}

double milliseconds(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now()-start).count();
}

struct Snapshot { std::vector<Track3D> cv, adsb; };
struct Benchmark { double reference_ms = 0.0, optimized_ms = 0.0; std::size_t frames = 0; };
volatile std::size_t benchmark_sink = 0;

std::vector<std::pair<int,int>> pairs(const Association& association) {
    std::vector<std::pair<int,int>> result;
    for (const auto& match : association.matches) result.emplace_back(match.cv_id, match.adsb_id);
    std::sort(result.begin(), result.end());
    return result;
}

Benchmark benchmark(const SensorFusion& fusion, const std::vector<Snapshot>& samples) {
    for (const auto& sample : samples) {
        if (pairs(fusion.associateTracks(sample.cv, sample.adsb)) !=
            pairs(fusion.associateReference(sample.cv, sample.adsb))) {
            throw std::runtime_error("Benchmark: reference and optimized assignments differ");
        }
    }
    std::vector<double> reference, optimized;
    for (int repeat = 0; repeat < 5; ++repeat) {
        for (int turn = 0; turn < 2; ++turn) {
            const bool baseline = (repeat + turn) % 2 == 0;
            const auto start = Clock::now();
            for (int pass = 0; pass < 20; ++pass) {
                for (const auto& sample : samples) {
                    const Association result = baseline ? fusion.associateReference(sample.cv, sample.adsb)
                                                         : fusion.associateTracks(sample.cv, sample.adsb);
                    benchmark_sink += result.matches.size();
                }
            }
            (baseline ? reference : optimized).push_back(milliseconds(start)/20.0);
        }
    }
    std::sort(reference.begin(), reference.end());
    std::sort(optimized.begin(), optimized.end());
    return {reference[2], optimized[2], samples.size()};
}

struct Metrics {
    std::size_t frames = 0, observations = 0, expected = 0, matches = 0, correct = 0;
    std::size_t cooperative = 0;
    std::size_t cv_only = 0, adsb_only = 0, switches = 0;
    double raw_error_sq = 0.0, error_sq = 0.0;
    double read_ms = 0.0, model_ms = 0.0, geometry_ms = 0.0, sync_ms = 0.0, association_ms = 0.0;
    std::vector<double> errors;
};

double ratio(std::size_t numerator, std::size_t denominator) {
    return denominator == 0 ? std::numeric_limits<double>::quiet_NaN()
                            : static_cast<double>(numerator) / denominator;
}

double rmse(double sum, std::size_t n) {
    return n == 0 ? std::numeric_limits<double>::quiet_NaN() : std::sqrt(sum / n);
}

double percentile95(std::vector<double>& values) {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    const std::size_t index = static_cast<std::size_t>(std::ceil(0.95*values.size())) - 1;
    std::nth_element(values.begin(), values.begin()+index, values.end());
    return values[index];
}

const Track3D& findTrack(const std::vector<Track3D>& tracks, int id) {
    for (const auto& track : tracks) if (track.id == id) return track;
    throw std::logic_error("Association refers to an unknown track");
}

void writeObservation(std::ostream& out, double time, const char* status,
                      const Track3D* cv, const Track3D* adsb, const Match* match, bool correct) {
    out << time << ',' << status << ',';
    if (cv) out << cv->id;
    out << ',';
    if (adsb) out << adsb->id;
    out << ',';
    if (match) out << match->distance;
    out << ',';
    if (match) out << match->score;
    out << ',';
    if (match) out << correct;
    for (const auto* track : {cv, adsb}) {
        out << ',';
        if (track) out << track->position.x << ',' << track->position.y << ',' << track->position.z;
        else out << ",,";
    }
    out << '\n';
}

void demo() {
    Track3D midpoint;
    const Track3D left = {1.0, 101, {100.0, 50.0, 10.0}};
    const Track3D right = {2.0, 101, {110.0, 50.0, 30.0}};
    if (!interpolateAdsb(1.5, left, right, 2.0, midpoint)) throw std::runtime_error("Demo interpolation failed");
    const std::vector<Track3D> cv = {{1.5, 1, {104.5, 49.5, 20.5}},
                                   {1.5, 2, {-49.0, 24.0, 8.0}},
                                   {1.5, 3, {500.0, 500.0, 40.0}}};
    const std::vector<Track3D> adsb = {midpoint, {1.5, 102, {-50.0, 25.0, 10.0}}};
    const Association result = SensorFusion().associateTracks(cv, adsb);
    std::cout << "3D example; coordinates in metres\n";
    for (const auto& match : result.matches) {
        std::cout << "CV " << match.cv_id << " -> ADS-B " << match.adsb_id
                  << ": distance=" << match.distance << ", score=" << match.score << '\n';
    }
    std::cout << "CV only: " << result.cv_only.size() << "\n";
}

void run(const Options& options) {
    const auto read_start = Clock::now();
    const LogData log = readLog(options.input);
    Metrics metrics;
    metrics.read_ms = milliseconds(read_start);
    const std::size_t end = options.limit_frames == 0 ? log.frames.size()
                            : std::min(options.limit_frames, log.frames.size());
    if (options.calibration_frames >= end) {
        throw std::invalid_argument("Calibration must leave evaluation frames; reduce --calibration-frames");
    }
    const LocalFrame local(log.frames.front().own);
    const Camera camera = calibrateCamera(log.frames, options.calibration_frames, options.camera);
    const SensorFusion fusion(options.max_distance);
    const auto model_start = Clock::now();
    const ModelData model = simulateAdsb(log.frames, local, options.simulation);
    std::set<int> cooperative_ids;
    for (const auto& entry : model.reference_ids) cooperative_ids.insert(entry.second);
    metrics.model_ms = milliseconds(model_start);
    AdsbSynchronizer synchronizer(model, options.max_age, options.max_gap, options.offline);
    std::ofstream matches;
    if (!options.matches.empty()) {
        matches.open(options.matches);
        if (!matches) throw std::runtime_error("Cannot write matches: " + options.matches);
        matches << "timestamp_s,status,cv_id,adsb_id,distance_m,score,correct,cv_east_m,cv_north_m,cv_up_m,adsb_east_m,adsb_north_m,adsb_up_m\n"
                << std::fixed << std::setprecision(6);
    }
    std::map<int,int> previous;
    std::vector<Snapshot> samples;
    for (std::size_t i = options.calibration_frames; i < end; ++i) {
        const Frame& frame = log.frames[i];
        const auto geometry_start = Clock::now();
        const CameraTransform corrected(frame.own,frame.attitude,camera,local);
        const CameraTransform uncalibrated(frame.own,frame.attitude,options.camera,local);
        std::vector<Track3D> cv;
        cv.reserve(frame.obstacles.size());
        for (const auto& obstacle : frame.obstacles) {
            const Vec3 position = corrected.apply(obstacle.position, obstacle.height);
            const Vec3 raw = uncalibrated.apply(obstacle.position, obstacle.height);
            const Vec3 reference = local.fromGeodetic(obstacle.reference);
            const double error = squaredDistance(position, reference);
            metrics.error_sq += error;
            metrics.raw_error_sq += squaredDistance(raw, reference);
            metrics.errors.push_back(std::sqrt(error));
            cv.push_back({frame.timestamp, obstacle.id, position});
        }
        metrics.geometry_ms += milliseconds(geometry_start);
        const auto sync_start = Clock::now();
        const auto adsb = synchronizer.at(frame.timestamp);
        metrics.sync_ms += milliseconds(sync_start);
        std::set<int> available;
        for (const auto& track : adsb) available.insert(model.reference_ids.at(track.id));
        for (const auto& track : cv) {
            metrics.expected += available.count(track.id);
            metrics.cooperative += cooperative_ids.count(track.id);
        }
        const auto association_start = Clock::now();
        const Association result = fusion.associateTracks(cv, adsb, options.assignment);
        metrics.association_ms += milliseconds(association_start);
        std::map<int,int> current;
        for (const auto& match : result.matches) {
            const bool correct = model.reference_ids.at(match.adsb_id) == match.cv_id;
            metrics.correct += correct;
            const auto last = previous.find(match.cv_id);
            if (last != previous.end() && last->second != match.adsb_id) ++metrics.switches;
            current[match.cv_id] = match.adsb_id;
            if (matches.is_open()) {
                writeObservation(matches,frame.timestamp,"matched",&findTrack(cv,match.cv_id),
                                  &findTrack(adsb,match.adsb_id),&match,correct);
            }
        }
        if (matches.is_open()) {
            for (int id : result.cv_only) {
                writeObservation(matches,frame.timestamp,"cv_only",&findTrack(cv,id),nullptr,nullptr,false);
            }
            for (int id : result.adsb_only) {
                writeObservation(matches,frame.timestamp,"adsb_only",nullptr,&findTrack(adsb,id),nullptr,false);
            }
        }
        previous.swap(current);
        ++metrics.frames;
        metrics.observations += cv.size();
        metrics.matches += result.matches.size();
        metrics.cv_only += result.cv_only.size();
        metrics.adsb_only += result.adsb_only.size();
        const std::size_t stride = (end-options.calibration_frames+499)/500;
        if (options.benchmark && (i-options.calibration_frames)%stride == 0) samples.push_back({cv, adsb});
    }
    if (matches.is_open()) {
        matches.flush();
        if (!matches) throw std::runtime_error("Error writing matches");
    }
    const double p95 = percentile95(metrics.errors);
    const Benchmark timing = options.benchmark ? benchmark(fusion, samples) : Benchmark{};
    const char* method = options.assignment == Assignment::Greedy ? "greedy" : "optimal";
    const char* mode = options.offline ? "offline" : "causal";
    std::cout << std::fixed << std::setprecision(4)
              << "Source: CV log + SYNTHETIC ADS-B (seed " << options.simulation.seed << ")\n"
              << "Frames: " << log.frames.size() << ", records: " << log.records
              << ", empty: " << log.empty_records << ", short empty: " << log.short_empty_records << '\n'
              << "Calibration frames: " << options.calibration_frames << " (excluded from evaluation)\n"
              << "Camera offset [right, up, forward], m: " << camera.offset.x << ", " << camera.offset.y << ", " << camera.offset.z << '\n'
              << "Evaluated frames: " << metrics.frames << ", CV observations: " << metrics.observations << '\n'
              << "Synthetic objects: " << model.reference_ids.size() << ", messages: " << model.messages << '\n'
              << "Assignment: " << method << ", synchronization: " << mode << '\n'
              << "Geometry RMSE before/after calibration, m: " << rmse(metrics.raw_error_sq, metrics.observations)
              << " / " << rmse(metrics.error_sq, metrics.observations) << ", p95: " << p95 << '\n'
              << "Matches: " << metrics.matches << ", correct: " << metrics.correct
              << ", wrong: " << metrics.matches-metrics.correct << '\n'
              << "Available reference pairs: " << metrics.expected
              << ", precision: " << ratio(metrics.correct, metrics.matches)
              << ", recall: " << ratio(metrics.correct, metrics.expected) << '\n'
              << "Cooperative CV observations: " << metrics.cooperative
              << ", availability: " << ratio(metrics.expected,metrics.cooperative)
              << ", end-to-end recall: " << ratio(metrics.correct,metrics.cooperative) << '\n'
              << "Unmatched CV/ADS-B observations: " << metrics.cv_only << " / " << metrics.adsb_only
              << ", ID switches between adjacent frames: " << metrics.switches << '\n'
              << "Read/model/geometry/sync/association, ms: " << metrics.read_ms << " / " << metrics.model_ms
              << " / " << metrics.geometry_ms << " / " << metrics.sync_ms << " / " << metrics.association_ms << '\n'
              << "Association per frame, ms: " << metrics.association_ms/metrics.frames << '\n';
    if (options.offline) std::cout << "Offline mode ignores delivery delay and can use future measurements.\n";
    if (options.benchmark) {
        std::cout << "Greedy benchmark: " << timing.frames << " frames, 20 passes/run, median of 5 runs\n"
                  << "Reference / optimized, ms: " << timing.reference_ms << " / " << timing.optimized_ms
                  << ", ratio: " << timing.reference_ms/timing.optimized_ms << '\n';
    }
    if (!options.summary.empty()) {
        std::ofstream out(options.summary);
        if (!out) throw std::runtime_error("Cannot write summary: " + options.summary);
        out << "source,assignment,sync,seed,calibration_frames,evaluation_frames,observations,period_s,horizontal_sigma_m,vertical_sigma_m,delay_s,dropout,coverage,max_distance_m,max_age_s,max_gap_s,camera_right_m,camera_up_m,camera_forward_m,mount_yaw_rad,mount_pitch_rad,mount_roll_rad,position_is_center,synthetic_objects,synthetic_messages,matches,correct,wrong,expected,precision,recall,cv_only,adsb_only,id_switches,raw_rmse_m,rmse_m,p95_m,read_ms,model_ms,geometry_ms,sync_ms,association_ms,benchmark_frames,reference_ms,optimized_ms,cooperative_observations,availability,end_to_end_recall\n";
        out << std::setprecision(12) << "synthetic_adsb," << method << ',' << mode << ',' << options.simulation.seed
            << ',' << options.calibration_frames << ',' << metrics.frames << ',' << metrics.observations
            << ',' << options.simulation.period << ',' << options.simulation.horizontal_sigma << ',' << options.simulation.vertical_sigma
            << ',' << options.simulation.delay << ',' << options.simulation.dropout << ',' << options.simulation.coverage
            << ',' << options.max_distance << ',' << options.max_age << ',' << options.max_gap
            << ',' << camera.offset.x << ',' << camera.offset.y << ',' << camera.offset.z
            << ',' << camera.mounting.yaw << ',' << camera.mounting.pitch << ',' << camera.mounting.roll
            << ',' << camera.position_is_center << ',' << model.reference_ids.size() << ',' << model.messages
            << ',' << metrics.matches << ',' << metrics.correct << ',' << metrics.matches-metrics.correct << ',' << metrics.expected
            << ',' << ratio(metrics.correct,metrics.matches) << ',' << ratio(metrics.correct,metrics.expected)
            << ',' << metrics.cv_only << ',' << metrics.adsb_only << ',' << metrics.switches
            << ',' << rmse(metrics.raw_error_sq,metrics.observations) << ',' << rmse(metrics.error_sq,metrics.observations)
            << ',' << p95 << ',' << metrics.read_ms << ',' << metrics.model_ms << ',' << metrics.geometry_ms
            << ',' << metrics.sync_ms << ',' << metrics.association_ms << ',' << timing.frames
            << ',' << timing.reference_ms << ',' << timing.optimized_ms
            << ',' << metrics.cooperative << ',' << ratio(metrics.expected,metrics.cooperative)
            << ',' << ratio(metrics.correct,metrics.cooperative) << '\n';
        out.flush();
        if (!out) throw std::runtime_error("Error writing summary");
    }
}
} // namespace

int main(int argc, char** argv) {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    try {
        if (argc == 1) demo();
        else if (argc == 2 && std::string(argv[1]) == "--help") help();
        else run(parse(argc, argv));
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
