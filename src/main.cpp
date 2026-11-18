#include "fusion.hpp"
#include "geometry.hpp"
#include "log_reader.hpp"
#include "simulation.hpp"
#include "report.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
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
#else
#include <fcntl.h>
#include <unistd.h>
#endif

namespace {
using Clock = std::chrono::steady_clock;

struct Options {
    std::string input, summary, matches, pairs_file, model_log, cv_log;
    std::size_t calibration_frames = 0;
    std::size_t limit_frames = 0;
    double max_distance = 15.0;
    double max_age = 2.5;
    double max_gap = 2.5;
    bool offline = false;
    bool benchmark = false;
    bool quiet = false;
    bool all_overlaps = false;
    bool overlap_metric = true;
    double overlap_threshold = 0.1;
    Assignment assignment = Assignment::Greedy;
    Camera camera = {{0,1.1,15.5},{},true};
    SimulationConfig simulation;
};

/** Читает конечное вещественное значение параметра.
 * @param text Строковое значение параметра командной строки.
 * @return Числовое значение аргумента.
 */
double numeric(const std::string& text) {
    std::size_t used = 0;
    const double value = std::stod(text, &used);
    if (used != text.size() || !std::isfinite(value)) throw std::invalid_argument("Invalid number: " + text);
    return value;
}

/** Читает неотрицательный целочисленный параметр.
 * @param text Строковое значение параметра командной строки.
 * @return Количество кадров или значение seed.
 */
std::size_t count(const std::string& text) {
    const double value = numeric(text);
    if (value < 0.0 || value != std::floor(value) || value > 1000000000.0) {
        throw std::invalid_argument("Expected a nonnegative integer: " + text);
    }
    return static_cast<std::size_t>(value);
}

/** Печатает справку командной строки. */
void help() {
    std::cout <<
        "Usage: fusion_core [--demo | --input log.txt] [options]\n"
        "Without input: 12 frames with two moving objects (aircraft and vehicle).\n\n"
        "  --metric MODE          overlap or distance (overlap)\n"
        "  --overlap-threshold N  Minimum positive IoU for output and matching, 0..1 (0.1)\n"
        "  --quiet                Print totals without individual pairs\n"
        "  --all-overlaps         Print positive pairs even below threshold\n"
        "  --pairs FILE           Write ALL positive overlaps, including rejected pairs\n"
        "  --model-log FILE       Write synthetic ADS-B times, sizes and directions\n"
        "  --cv-log FILE          Write normalized CV observations\n"
        "  --calibration-frames N  Estimate offset instead of known camera mounting (0)\n"
        "  --limit-frames N        Evaluate only the first N frames; 0 means all\n"
        "  --camera-right M        Manual offset to the right (0)\n"
        "  --camera-up M           Upward offset: block 1.1 + camera 0 (1.1)\n"
        "  --camera-forward M      Forward offset: block 15.5 + camera 0 (15.5)\n"
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
        "  --horizontal-sigma M    Horizontal noise per axis (log 3, demo 0.5)\n"
        "  --vertical-sigma M      Vertical noise (log 2, demo 0.1)\n"
        "  --heading-sigma DEG     Synthetic heading noise (0)\n"
        "  --delay SEC             Synthetic delivery delay (0.2)\n"
        "  --dropout FRACTION      Message loss, 0..1 (log 0.1, demo 0)\n"
        "  --coverage FRACTION     Objects with ADS-B (log 0.8, demo 1)\n"
        "  --seed N                Random seed (42)\n"
        "  --benchmark             Compare reference and optimized greedy code\n"
        "  --summary FILE          Write settings and metrics as CSV\n"
        "  --matches FILE          Write matches and unmatched objects as CSV\n";
}

/** Разбирает параметры запуска и проверяет пути выходных файлов.
 * @param argc Число аргументов командной строки.
 * @param argv Массив аргументов командной строки.
 * @return Проверенные настройки запуска.
 */
Options parse(int argc, char** argv) {
    Options options;
    std::set<std::string> supplied;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        supplied.insert(key);
        if (key == "--demo") { continue; }
        if (key == "--quiet") { options.quiet = true; continue; }
        if (key == "--all-overlaps") { options.all_overlaps = true; continue; }
        if (key == "--offline") { options.offline = true; continue; }
        if (key == "--benchmark") { options.benchmark = true; continue; }
        if (i + 1 >= argc) throw std::invalid_argument("Missing value for " + key);
        const std::string value = argv[++i];
        if (key == "--input") options.input = value;
        else if (key == "--summary") options.summary = value;
        else if (key == "--matches") options.matches = value;
        else if (key == "--pairs") options.pairs_file = value;
        else if (key == "--model-log") options.model_log = value;
        else if (key == "--cv-log") options.cv_log = value;
        else if (key == "--overlap-threshold") options.overlap_threshold = numeric(value);
        else if (key == "--metric" && value == "overlap") options.overlap_metric = true;
        else if (key == "--metric" && value == "distance") options.overlap_metric = false;
        else if (key == "--calibration-frames") options.calibration_frames = count(value);
        else if (key == "--limit-frames") options.limit_frames = count(value);
        else if (key == "--seed") options.simulation.seed = static_cast<std::uint32_t>(count(value));
        else if (key == "--max-distance") options.max_distance = numeric(value);
        else if (key == "--max-age") options.max_age = numeric(value);
        else if (key == "--max-gap") options.max_gap = numeric(value);
        else if (key == "--period") options.simulation.period = numeric(value);
        else if (key == "--horizontal-sigma") options.simulation.horizontal_sigma = numeric(value);
        else if (key == "--vertical-sigma") options.simulation.vertical_sigma = numeric(value);
        else if (key == "--heading-sigma") options.simulation.heading_sigma = radians(numeric(value));
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
    if (supplied.count("--demo") && !options.input.empty()) {
        throw std::invalid_argument("Choose --demo or --input, not both");
    }
    if (options.input.empty()) {
        if (!supplied.count("--horizontal-sigma")) options.simulation.horizontal_sigma = 0.5;
        if (!supplied.count("--vertical-sigma")) options.simulation.vertical_sigma = 0.1;
        if (!supplied.count("--coverage")) options.simulation.coverage = 1;
        if (!supplied.count("--dropout")) options.simulation.dropout = 0;
    }
    if (options.overlap_threshold < 0 || options.overlap_threshold > 1) {
        throw std::invalid_argument("Overlap threshold must be in [0,1]");
    }
    std::set<std::string> outputs;
    for (const auto& output : {options.summary,options.matches,options.pairs_file,options.model_log,options.cv_log}) {
        if (!output.empty() && (output == options.input || !outputs.insert(output).second)) {
            throw std::invalid_argument("Input and output paths must be different");
        }
        if (!output.empty() && std::ifstream(output).good()) {
            throw std::invalid_argument("Output already exists; choose a new filename: " + output);
        }
    }
    return options;
}

/** Измеряет интервал монотонных часов.
 * @param start Начальная отметка монотонных часов.
 * @return Прошедшее время в миллисекундах.
 */
double milliseconds(Clock::time_point start) {
    return std::chrono::duration<double, std::milli>(Clock::now()-start).count();
}

struct Snapshot { std::vector<Track3D> cv, adsb; };
struct Benchmark { double reference_ms = 0.0, optimized_ms = 0.0; std::size_t frames = 0, passes = 0; };
volatile std::size_t benchmark_sink = 0;

/** Выделяет упорядоченные ID соответствий для сравнения вариантов расчёта.
 * @param association Результат назначения пар.
 * @return Пары ID, отсортированные независимо от порядка выдачи.
 */
std::vector<std::pair<int,int>> pairs(const Association& association) {
    std::vector<std::pair<int,int>> result;
    for (const auto& match : association.matches) result.emplace_back(match.cv_id, match.adsb_id);
    std::sort(result.begin(), result.end());
    return result;
}

/** Сравнивает исходный и оптимизированный расчёты на одинаковых кадрах.
 * @param fusion Алгоритм сопоставления с выбранным порогом расстояния.
 * @param samples Одинаковые наборы входных объектов для обоих вариантов расчёта.
 * @return Медианные времена одного прохода и число кадров.
 */
Benchmark benchmark(const SensorFusion& fusion, const std::vector<Snapshot>& samples) {
    for (const auto& sample : samples) {
        if (pairs(fusion.associateTracks(sample.cv, sample.adsb)) !=
            pairs(fusion.associateReference(sample.cv, sample.adsb))) {
            throw std::runtime_error("Benchmark: reference and optimized assignments differ");
        }
    }
    // На коротком примере число повторений повышается до измеримого интервала часов.
    std::size_t passes = 20;
    for (;;) {
        const auto start = Clock::now();
        for (std::size_t pass=0; pass<passes; ++pass) {
            for (const auto& sample : samples) benchmark_sink += fusion.associateReference(sample.cv,sample.adsb).matches.size();
        }
        if (milliseconds(start) >= 20 || passes >= 327680) break;
        passes *= 2;
    }
    std::vector<double> reference, optimized;
    for (int repeat = 0; repeat < 5; ++repeat) {
        for (int turn = 0; turn < 2; ++turn) {
            const bool baseline = (repeat + turn) % 2 == 0;
            const auto start = Clock::now();
            for (std::size_t pass = 0; pass < passes; ++pass) {
                for (const auto& sample : samples) {
                    const Association result = baseline ? fusion.associateReference(sample.cv, sample.adsb)
                                                         : fusion.associateTracks(sample.cv, sample.adsb);
                    benchmark_sink += result.matches.size();
                }
            }
            (baseline ? reference : optimized).push_back(milliseconds(start)/passes);
        }
    }
    std::sort(reference.begin(), reference.end());
    std::sort(optimized.begin(), optimized.end());
    return {reference[2], optimized[2], samples.size(),passes};
}

struct Metrics {
    std::size_t frames = 0, observations = 0, expected = 0, matches = 0, correct = 0;
    std::size_t cooperative = 0;
    std::size_t adsb_observations = 0, positive_pairs = 0, threshold_pairs = 0, eligible_pairs = 0;
    std::set<int> cv_ids, adsb_ids;
    std::size_t cv_only = 0, adsb_only = 0, switches = 0;
    double raw_error_sq = 0.0, error_sq = 0.0;
    double read_ms = 0.0, model_ms = 0.0, geometry_ms = 0.0, sync_ms = 0.0, association_ms = 0.0;
    std::vector<double> errors;
};

/** Вычисляет долю для показателей точности, полноты и доступности.
 * @param numerator Количество успешных событий.
 * @param denominator Общее количество событий.
 * @return Отношение чисел или NaN при нулевом знаменателе.
 */
double ratio(std::size_t numerator, std::size_t denominator) {
    return denominator == 0 ? std::numeric_limits<double>::quiet_NaN()
                            : static_cast<double>(numerator) / denominator;
}

/** Вычисляет среднеквадратическую трёхмерную ошибку.
 * @param sum Сумма квадратов трёхмерных ошибок.
 * @param n Число наблюдений, вошедших в сумму квадратов.
 * @return RMSE либо NaN для пустой выборки.
 */
double rmse(double sum, std::size_t n) {
    return n == 0 ? std::numeric_limits<double>::quiet_NaN() : std::sqrt(sum / n);
}

/** Находит 95-й процентиль ошибок без полной сортировки.
 * @param values Числовой набор; порядок элементов может измениться.
 * @return Граница, не превышаемая 95 процентами ошибок, либо NaN.
 */
double percentile95(std::vector<double>& values) {
    if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
    const std::size_t index = static_cast<std::size_t>(std::ceil(0.95*values.size())) - 1;
    std::nth_element(values.begin(), values.begin()+index, values.end());
    return values[index];
}

/** Резервирует имена всех CSV до начала записи.
 * @param paths Запрошенные файлы отчёта. Пустой путь отключает выгрузку.
 */
void reserveOutputs(const std::vector<std::string>& paths) {
    std::vector<const std::string*> created;
    created.reserve(paths.size());
    for (const auto& path : paths) {
        if (path.empty()) continue;
        // ОС проверяет и создаёт файл за один вызов, в том числе при двух одновременных запусках.
#ifdef _WIN32
        HANDLE file = CreateFileA(path.c_str(), GENERIC_WRITE, 0, nullptr, CREATE_NEW, FILE_ATTRIBUTE_NORMAL, nullptr);
        const bool reserved = file != INVALID_HANDLE_VALUE;
        if (reserved) CloseHandle(file);
#else
        const int file = open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL, 0666);
        const bool reserved = file >= 0;
        if (reserved) close(file);
#endif
        if (reserved) {
            created.push_back(&path);
            continue;
        }
        // Если два пути указывают на один файл, до этой точки записаны только пустые файлы.
        for (const auto* empty_file : created) std::remove(empty_file->c_str());
        throw std::runtime_error("Cannot create new CSV: " + path);
    }
}

/** Записывает заголовок в заранее созданный CSV.
 * @param stream Выходной поток, который будет открыт.
 * @param path Зарезервированный путь либо пустая строка.
 * @param header Первая строка с названиями столбцов.
 */
void openCsv(std::ofstream& stream, const std::string& path, const std::string& header) {
    if (path.empty()) return;
    stream.open(path);
    if (!stream) throw std::runtime_error("Cannot write CSV: "+path);
    stream << header << '\n' << std::setprecision(std::numeric_limits<double>::max_digits10);
}

/** Проверяет завершение записи CSV.
 * @param stream Открытый поток либо отключённая выгрузка.
 * @throws std::runtime_error При ошибке записи на диск.
 */
void finishCsv(std::ofstream& stream) {
    if (!stream.is_open()) return;
    stream.flush();
    if (!stream) throw std::runtime_error("Error writing CSV");
}

/** Выполняет моделирование, сопоставление, вывод пар и подсчёт итогов.
 * @param options Параметры источников, камеры, сопоставления и выгрузки.
 */
void run(const Options& options) {
    // Исходный лог либо воспроизводимый пример с двумя физическими объектами.
    const auto read_start = Clock::now();
    const LogData log = options.input.empty() ? twoObjectScenario(options.camera) : readLog(options.input);
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
    // Таблица опорных ID остаётся в оценке качества и не передаётся сопоставлению.
    const auto model_start = Clock::now();
    const ModelData model = simulateAdsb(log.frames, local, options.simulation);
    std::set<int> cooperative_ids;
    for (const auto& entry : model.reference_ids) cooperative_ids.insert(entry.second);
    metrics.model_ms = milliseconds(model_start);
    AdsbSynchronizer synchronizer(model, options.max_age, options.max_gap, options.offline);
    reserveOutputs({options.summary, options.matches, options.pairs_file, options.model_log, options.cv_log});
    std::ofstream matches, overlaps, cv_log, model_log;
    openCsv(matches,options.matches,"timestamp_s,status,"+trackCsvHeader("cv_")+","+trackCsvHeader("adsb_")+
                                     ",distance_m,score,iou,intersection_m2,correct");
    openCsv(overlaps,options.pairs_file,"timestamp_s,"+trackCsvHeader("cv_")+","+trackCsvHeader("adsb_")+
                                      ",intersection_m2,iou,overlap_positive,above_threshold,height_overlap,same_class,assigned,correct");
    openCsv(cv_log,options.cv_log,"timestamp_s,"+trackCsvHeader(""));
    openCsv(model_log,options.model_log,"measurement_time_s,received_at_s,"+trackCsvHeader(""));
    if (model_log.is_open()) {
        for (const auto& stream : model.streams) {
            for (const auto& message : stream.second) {
                model_log << message.track.timestamp << ',' << message.received_at << ',';
                writeTrackCsv(model_log,&message.track);
                model_log << '\n';
            }
        }
    }
    std::cout << std::defaultfloat << std::setprecision(6)
              << "Source: " << (options.input.empty() ? "TWO-OBJECT SCENARIO" : "CV log")
              << " + SYNTHETIC ADS-B (seed " << options.simulation.seed << ")\n";
    std::map<int,int> previous;
    std::vector<Snapshot> samples;
    for (std::size_t i = options.calibration_frames; i < end; ++i) {
        const Frame& frame = log.frames[i];
        const auto geometry_start = Clock::now();
        const CameraTransform configured(frame.own, frame.attitude, options.camera, local);
        const CameraTransform corrected = options.calibration_frames == 0 ? configured
            : CameraTransform(frame.own, frame.attitude, camera, local);
        std::vector<Track3D> cv;
        cv.reserve(frame.obstacles.size());
        for (const auto& obstacle : frame.obstacles) {
            const Track3D track = cameraTrack(obstacle,frame.timestamp,corrected);
            const Vec3 position = track.position;
            const Vec3 reference = local.fromGeodetic(obstacle.reference);
            const double error = squaredDistance(position, reference);
            metrics.error_sq += error;
            metrics.raw_error_sq += options.calibration_frames == 0 ? error
                : squaredDistance(configured.apply(obstacle.position, obstacle.size.y), reference);
            metrics.errors.push_back(std::sqrt(error));
            cv.push_back(track);
            if (cv_log.is_open()) {
                cv_log << frame.timestamp << ',';
                writeTrackCsv(cv_log,&track);
                cv_log << '\n';
            }
            metrics.cv_ids.insert(track.id);
        }
        metrics.geometry_ms += milliseconds(geometry_start);
        const auto sync_start = Clock::now();
        const auto adsb = synchronizer.at(frame.timestamp);
        metrics.sync_ms += milliseconds(sync_start);
        std::set<int> available;
        for (const auto& track : adsb) {
            available.insert(model.reference_ids.at(track.id));
            metrics.adsb_ids.insert(track.id);
        }
        for (const auto& track : cv) {
            metrics.expected += available.count(track.id);
            metrics.cooperative += cooperative_ids.count(track.id);
        }
        const auto association_start = Clock::now();
        const Association result = options.overlap_metric
            ? fusion.associateByOverlap(cv,adsb,options.overlap_threshold,options.assignment)
            : fusion.associateTracks(cv,adsb,options.assignment);
        metrics.association_ms += milliseconds(association_start);
        for (const auto& pair : result.overlaps) {
            ++metrics.positive_pairs;
            metrics.threshold_pairs += pair.above_threshold;
            metrics.eligible_pairs += pair.above_threshold && pair.height_overlap && pair.same_class;
            const bool show_pair = !options.quiet && (pair.above_threshold || options.all_overlaps);
            if (!show_pair && !overlaps.is_open()) continue;
            const auto& left = findTrack(cv,pair.cv_id);
            const auto& right = findTrack(adsb,pair.adsb_id);
            if (show_pair) {
                printOverlapPair(std::cout,frame.timestamp,left,right,pair);
            }
            if (overlaps.is_open()) writeOverlapPair(overlaps,frame.timestamp,left,right,pair,
                                                      model.reference_ids.at(pair.adsb_id) == pair.cv_id);
        }
        std::map<int,int> current;
        // Правильность назначения проверяется после работы алгоритма по скрытым опорным ID.
        for (const auto& match : result.matches) {
            const bool correct = model.reference_ids.at(match.adsb_id) == match.cv_id;
            metrics.correct += correct;
            const auto last = previous.find(match.cv_id);
            if (last != previous.end() && last->second != match.adsb_id) ++metrics.switches;
            current[match.cv_id] = match.adsb_id;
            if (matches.is_open()) {
                writeAssociationRow(matches,frame.timestamp,"matched",&findTrack(cv,match.cv_id),
                                  &findTrack(adsb,match.adsb_id),&match,correct);
            }
        }
        if (matches.is_open()) {
            for (int id : result.cv_only) {
                writeAssociationRow(matches,frame.timestamp,"cv_only",&findTrack(cv,id),nullptr,nullptr,false);
            }
            for (int id : result.adsb_only) {
                writeAssociationRow(matches,frame.timestamp,"adsb_only",nullptr,&findTrack(adsb,id),nullptr,false);
            }
        }
        previous.swap(current);
        ++metrics.frames;
        metrics.observations += cv.size();
        metrics.adsb_observations += adsb.size();
        metrics.matches += result.matches.size();
        metrics.cv_only += result.cv_only.size();
        metrics.adsb_only += result.adsb_only.size();
        const std::size_t stride = (end-options.calibration_frames+499)/500;
        if (options.benchmark && (i-options.calibration_frames)%stride == 0) samples.push_back({cv, adsb});
    }
    finishCsv(matches);
    finishCsv(overlaps);
    finishCsv(cv_log);
    finishCsv(model_log);
    const double p95 = percentile95(metrics.errors);
    const Benchmark timing = options.benchmark ? benchmark(fusion, samples) : Benchmark{};
    const char* method = options.assignment == Assignment::Greedy ? "greedy" : "optimal";
    const char* mode = options.offline ? "offline" : "causal";
    std::cout << std::fixed << std::setprecision(4)
              << "Frames: " << log.frames.size() << ", records: " << log.records
              << ", empty: " << log.empty_records << ", short empty: " << log.short_empty_records << '\n'
              << "Calibration frames: " << options.calibration_frames << " (excluded from evaluation)\n"
              << "Camera offset [right, up, forward], m: " << camera.offset.x << ", " << camera.offset.y << ", " << camera.offset.z << '\n'
              << "Evaluated frames: " << metrics.frames << ", CV observations: " << metrics.observations << '\n'
              << "Synthetic objects: " << model.reference_ids.size() << ", messages: " << model.messages << '\n'
              << "Assignment: " << method << ", synchronization: " << mode << '\n'
              << "Metric: " << (options.overlap_metric ? "footprint IoU" : "distance")
              << ", overlap threshold: " << options.overlap_threshold << '\n'
              << "Positive overlap pairs: " << metrics.positive_pairs
              << ", at threshold: " << metrics.threshold_pairs << ", eligible: " << metrics.eligible_pairs << '\n'
              << "CV/ADS-B observations: " << metrics.observations << " / " << metrics.adsb_observations
              << ", unique IDs: " << metrics.cv_ids.size() << " / " << metrics.adsb_ids.size() << '\n'
              << "Geometry RMSE configured/final camera, m: " << rmse(metrics.raw_error_sq, metrics.observations)
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
        std::cout << "Greedy distance benchmark: " << timing.frames << " frames, " << timing.passes << " passes/run, median of 5 runs\n"
                  << "Reference / optimized, ms: " << timing.reference_ms << " / " << timing.optimized_ms
                  << ", ratio: " << timing.reference_ms/timing.optimized_ms << '\n';
    }
    if (!options.summary.empty()) {
        std::ofstream out(options.summary);
        if (!out) throw std::runtime_error("Cannot write summary: " + options.summary);
        out << "source,assignment,sync,seed,calibration_frames,evaluation_frames,observations,period_s,horizontal_sigma_m,vertical_sigma_m,delay_s,dropout,coverage,max_distance_m,max_age_s,max_gap_s,camera_right_m,camera_up_m,camera_forward_m,mount_yaw_rad,mount_pitch_rad,mount_roll_rad,position_is_center,synthetic_objects,synthetic_messages,matches,correct,wrong,expected,precision,recall,cv_only,adsb_only,id_switches,raw_rmse_m,rmse_m,p95_m,read_ms,model_ms,geometry_ms,sync_ms,association_ms,benchmark_frames,reference_ms,optimized_ms,cooperative_observations,availability,end_to_end_recall,metric,overlap_threshold,positive_pairs,threshold_pairs,eligible_pairs,adsb_observations,unique_cv,unique_adsb,heading_sigma_rad,benchmark_passes\n";
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
            << ',' << ratio(metrics.correct,metrics.cooperative)
            << ',' << (options.overlap_metric ? "overlap" : "distance") << ',' << options.overlap_threshold
            << ',' << metrics.positive_pairs << ',' << metrics.threshold_pairs << ',' << metrics.eligible_pairs
            << ',' << metrics.adsb_observations << ',' << metrics.cv_ids.size() << ',' << metrics.adsb_ids.size()
            << ',' << options.simulation.heading_sigma << ',' << timing.passes << '\n';
        out.flush();
        if (!out) throw std::runtime_error("Error writing summary");
    }
}
} // namespace

/** Запускает программу и преобразует исключения в сообщение об ошибке.
 * @param argc Число аргументов командной строки.
 * @param argv Массив аргументов командной строки.
 * @return 0 при успешном завершении, 1 при ошибке.
 */
int main(int argc, char** argv) {
#ifdef _WIN32
    SetConsoleOutputCP(CP_UTF8);
#endif
    try {
        if (argc == 2 && std::string(argv[1]) == "--help") help();
        else run(parse(argc, argv));
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
