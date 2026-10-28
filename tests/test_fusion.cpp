#include "fusion.hpp"
#include "geometry.hpp"
#include "log_reader.hpp"
#include "simulation.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
int checks = 0;

void require(bool condition, const std::string& message) {
    ++checks;
    if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const std::string& message) {
    require(std::abs(actual-expected) <= tolerance, message);
}

void near(const Vec3& actual, const Vec3& expected, double tolerance, const std::string& message) {
    require(squaredDistance(actual, expected) <= tolerance*tolerance, message);
}

template<class Function>
void rejects(Function function, const std::string& message) {
    bool thrown = false;
    try { function(); } catch (const std::exception&) { thrown = true; }
    require(thrown, message);
}

void interpolation() {
    const Track3D a = {1.0, 7, {0, 10, 20}}, b = {3.0, 7, {10, 20, 60}};
    Track3D result;
    require(interpolateAdsb(2, a, b, 3, result), "3D interpolation");
    near(result.position, {5,15,40}, 1e-12, "Interpolation includes height");
    require(result.id == 7 && result.timestamp == 2, "Interpolation identity and time");
    require(!interpolateAdsb(0, a, b, 3, result), "No implicit backward extrapolation");
    require(!interpolateAdsb(4, a, b, 3, result), "No implicit forward extrapolation");
    require(!interpolateAdsb(2, a, b, 1, result), "Long gap rejected");
    require(!interpolateAdsb(2, b, a, 3, result), "Unordered samples rejected");
    Track3D other = b;
    other.id = 8;
    require(!interpolateAdsb(2, a, other, 3, result), "Different source objects rejected");
    require(interpolateAdsb(1, a, a, 3, result), "Repeated identical sample");
    other = a;
    other.position.z += 1;
    require(!interpolateAdsb(1, a, other, 3, result), "Conflicting equal-time samples rejected");
    require(!interpolateAdsb(std::numeric_limits<double>::quiet_NaN(), a, b, 3, result), "NaN time rejected");
}

void association() {
    SensorFusion fusion;
    const std::vector<Track3D> cv = {{0, 10, {0,0,0}}, {0, 20, {100,100,0}}};
    auto result = fusion.associateTracks(cv, {{0, 900, {0,0,20}}, {0, 901, {100,100,2}}});
    require(result.matches.size() == 1 && result.matches[0].cv_id == 20 &&
            result.matches[0].adsb_id == 901, "Vertical separation gates association");
    require(result.cv_only == std::vector<int>{10} && result.adsb_only == std::vector<int>{900}, "Unmatched objects retained");
    near(result.matches[0].distance, 2, 1e-12, "Match distance");
    require(fusion.associateTracks({cv[0]}, {{0,1,{15,0,0}}}).matches.empty(), "Strict 15m boundary");
    require(fusion.associateTracks({cv[0]}, {{0,1,{std::nextafter(15.0,0.0),0,0}}}).matches.size() == 1, "Inside distance boundary");
    require(fusion.associateTracks({}, {}).matches.empty(), "Empty inputs");
    result = fusion.associateTracks({{0,99,{0,0,0}}}, {{0,99,{100,0,0}}, {0,-7,{0,0,0}}});
    require(result.matches.size() == 1 && result.matches[0].adsb_id == -7, "ID equality cannot reveal reference answer");
    rejects([&] { fusion.associateTracks({cv[0],cv[0]}, {}); }, "Duplicate IDs rejected");
    rejects([&] { fusion.associateTracks({cv[0]}, {{1,1,{0,0,0}}}); }, "Unsynchronized sources rejected");
    rejects([&] { fusion.associateTracks({{0,1,{std::numeric_limits<double>::infinity(),0,0}}}, {}); }, "Invalid coordinate rejected");
    rejects([] { SensorFusion invalid(0); }, "Zero gate rejected");

    const std::vector<Track3D> close_cv = {{0,1,{0,0,0}}, {0,2,{10,0,0}}};
    const std::vector<Track3D> close_adsb = {{0,101,{1,0,0}}, {0,102,{-10,0,0}}};
    require(fusion.associateTracks(close_cv, close_adsb).matches.size() == 1, "Greedy ambiguous example");
    require(fusion.associateTracks(close_cv, close_adsb, Assignment::Optimal).matches.size() == 2, "Global assignment recovers two pairs");
}

void optimalAgainstEnumeration() {
    std::mt19937 random(5);
    std::uniform_real_distribution<double> coordinate(-20,20);
    SensorFusion fusion;
    for (int trial = 0; trial < 100; ++trial) {
        const int n = trial % 5, m = (trial / 5) % 5;
        std::vector<Track3D> cv, adsb;
        for (int i=0; i<n; ++i) cv.push_back({0,i,{coordinate(random),coordinate(random),coordinate(random)}});
        for (int i=0; i<m; ++i) adsb.push_back({0,100+i,{coordinate(random),coordinate(random),coordinate(random)}});
        int best_count = -1;
        double best_cost = std::numeric_limits<double>::infinity();
        std::function<void(int,unsigned,int,double)> search = [&](int i, unsigned used, int matched, double cost) {
            if (i == n) {
                if (matched > best_count || (matched == best_count && cost < best_cost)) {
                    best_count = matched;
                    best_cost = cost;
                }
                return;
            }
            search(i+1,used,matched,cost);
            for (int j=0; j<m; ++j) {
                const double d2 = squaredDistance(cv[i].position,adsb[j].position);
                if (!(used & (1u<<j)) && d2 < 225) search(i+1,used|(1u<<j),matched+1,cost+d2);
            }
        };
        search(0,0,0,0);
        const auto result = fusion.associateTracks(cv,adsb,Assignment::Optimal);
        double cost = 0;
        for (const auto& match : result.matches) cost += match.distance*match.distance;
        require(static_cast<int>(result.matches.size()) == best_count, "Optimal cardinality agrees with enumeration");
        near(cost,best_cost,1e-8,"Optimal squared distance agrees with enumeration");
        const auto reference = fusion.associateReference(cv,adsb);
        const auto optimized = fusion.associateTracks(cv,adsb);
        require(reference.matches.size() == optimized.matches.size(), "Optimization retains greedy pair count");
        for (std::size_t i=0; i<reference.matches.size(); ++i) {
            require(reference.matches[i].cv_id == optimized.matches[i].cv_id &&
                    reference.matches[i].adsb_id == optimized.matches[i].adsb_id, "Optimization retains greedy IDs");
        }
    }
}

void geometry() {
    const Geodetic origin = {0,0,10};
    near(toEcef(origin), {6378147,0,0}, 1e-7, "WGS84 equator");
    near(toEcef({radians(90),0,0}), {0,0,6356752.314245179}, 1e-6, "WGS84 pole");
    const LocalFrame local(origin);
    near(local.fromGeodetic(origin), {0,0,0}, 1e-9, "Local origin");
    near(local.fromGeodetic({0,0,17}), {0,0,7}, 1e-8, "Ellipsoidal height to ENU up");
    const Geodetic arbitrary = {radians(43),radians(40),27};
    near(ecefToEnu(enuToEcef({13,25,70},arbitrary),arbitrary), {13,25,70}, 1e-10, "ECEF/ENU inverse");
    near(bodyToEnu({0,0,1}, {radians(90),0,0}), {1,0,0}, 1e-12, "Heading east");
    near(bodyToEnu({0,0,1}, {0,radians(90),0}), {0,0,1}, 1e-12, "Positive pitch raises nose");
    near(bodyToEnu({1,0,0}, {0,0,radians(90)}), {0,0,-1}, 1e-12, "Positive roll lowers right wing");
    const Attitude attitude = {radians(123),radians(20),radians(-30)};
    near(enuToBody(bodyToEnu({4,5,6},attitude),attitude), {4,5,6}, 1e-10, "Body transform inverse");
    near(cameraToBody({4,5,6},{}), {4,5,6}, 1e-12, "Zero mounting rotation");
    Camera camera;
    camera.offset = {1,2,3};
    near(cameraPosition({10,5,100},2,origin,{},camera,local), {11,103,6}, 1e-7, "Camera translation and object anchor");
    near(cameraOffset({0,5,0},2,{0,0,17},origin,{},Camera{}), {0,3,0}, 1e-7, "Recover camera translation");
    rejects([] { toEcef({43,40,0}); }, "Degrees cannot be interpreted as radians");

    Frame first;
    first.own = origin;
    first.obstacles.push_back({7,0,{0,5,0},2,{0,0,17}});
    Frame outlier = first;
    outlier.obstacles[0].position.y = 1000;
    const Camera calibrated = calibrateCamera({first,first,outlier},2,Camera{});
    near(calibrated.offset,{0,3,0},1e-7,"Calibration excludes evaluation frame");
    rejects([&] { calibrateCamera({first,first},2,Camera{}); }, "Calibration must leave held-out data");
}

const std::string header = "ID\tClass\tClassProbability\tSize.X\tSize.Y\tSize.Z\tPosition.X\tPosition.Y\tPosition.Z\tDirection.X\tDirection.Y\tDirection.Z\tSpatialProbability\tTimestamp(Ms)\tObs.Lat(rad)\tObs.Lon(rad)\tObs.Alt(M)\tObs.Pitch(Deg)\tObs.Roll(Deg)\tObs.Yaw(Deg)\tOwn.Lat(rad)\tOwn.Lon(rad)\tOwn.Alt(M)\tOwn.Pitch(Deg)\tOwn.Roll(Deg)\tOwn.Yaw(Deg)\tOwn.GroundSpeed(km/h)\tOwn.RoC(m/s)";

std::string record(int id, int timestamp) {
    return std::to_string(id)+"\t0\t0.8\t2\t2\t2\t10\t5\t100\t0\t0\t1\t0.8\t"+
           std::to_string(timestamp)+"\t0\t0\t17\t0\t0\t0\t0\t0\t10\t0\t0\t0\t10\t0";
}

LogData parse(const std::string& text) {
    std::istringstream input(text);
    return readLog(input);
}

void reader() {
    const std::string empty = std::string(13,'\t')+"2000"+std::string(7,'\t')+"0\t0\t10\t0\t0\t0";
    const auto log = parse("\xef\xbb\xbf"+header+"\r\n"+record(-7,1000)+"\r\n"+empty+"\r\n");
    require(log.frames.size() == 2 && log.frames[0].obstacles.size() == 1 && log.frames[1].obstacles.empty(), "Read obstacle and short empty frame");
    require(log.short_empty_records == 1 && log.empty_records == 1, "Short records counted");
    near(log.frames[0].timestamp,0,1e-12,"Recording-relative time");
    near(log.frames[1].timestamp,1,1e-12,"Milliseconds converted to seconds");
    rejects([&] { parse(header+"\n"+record(1,2000)+"\n"+record(2,1000)); }, "Unordered log time rejected");
    rejects([&] { parse(header+"\n"+record(1,1000)+"\n"+record(1,1000)); }, "Duplicate per-frame obstacle rejected");
    rejects([&] { parse("ID\tID\n1\t2"); }, "Duplicate header rejected");
    rejects([&] { parse("ID\tTimestamp(Ms)\n1\t1000"); }, "Missing geometry columns rejected");
    std::string invalid = record(1,1000);
    invalid.replace(invalid.find("\t10\t5\t100"),10,"\tnan\t5\t100");
    rejects([&] { parse(header+"\n"+invalid); }, "Non-finite log value rejected");
    rejects([&] { parse(header+"\n"+record(1,1000)+"\textra"); }, "Extra field rejected");
    rejects([&] { parse(header+"\n"+empty+"\n"+record(1,2000)); }, "Empty frame cannot contain obstacle record");
}

void synchronization() {
    ModelData model;
    model.streams[987] = {{{0,987,{0,0,0}},1}, {{2,987,{2,4,6}},3}};
    AdsbSynchronizer causal(model,2.5,3,false);
    require(causal.at(.5).empty(), "Delayed message unavailable before arrival");
    auto tracks = causal.at(1);
    require(tracks.size() == 1, "Message available on arrival");
    near(tracks[0].position,{0,0,0},1e-12,"Causal mode does not use future sample");
    tracks = causal.at(2);
    near(tracks[0].position,{0,0,0},1e-12,"Still waiting for delayed second sample");
    tracks = causal.at(3);
    near(tracks[0].position,{2,4,6},1e-12,"Second sample received");
    require(causal.at(5).empty(), "Stale measurement discarded");
    rejects([&] { causal.at(4); }, "Backward synchronization rejected");
    AdsbSynchronizer offline(model,2.5,3,true);
    tracks = offline.at(1);
    near(tracks[0].position,{1,2,3},1e-12,"Offline interpolation uses all three axes");
    AdsbSynchronizer gap(model,.5,1,true);
    require(gap.at(1).empty(), "Long gap cannot be bridged");
}

void simulation() {
    Frame first, last;
    first.own = {0,0,10};
    first.obstacles.push_back({17,0,{0,0,0},2,{0,0,17}});
    last = first;
    last.timestamp = 2;
    const std::vector<Frame> frames = {first,last};
    const LocalFrame local(first.own);
    SimulationConfig config;
    config.coverage = 1;
    config.dropout = 0;
    config.horizontal_sigma = 0;
    config.vertical_sigma = 0;
    const auto model = simulateAdsb(frames,local,config);
    require(model.messages == 3 && model.streams.size() == 1, "Synthetic message frequency");
    const int id = model.streams.begin()->first;
    require(id != 17 && model.reference_ids.at(id) == 17, "Independent synthetic IDs");
    near(model.streams.at(id)[1].track.position,{0,0,7},1e-7,"Synthetic positions from reference geometry");
    near(model.streams.at(id)[1].received_at,1.2,1e-12,"Measurement time differs from delivery time");
    config.horizontal_sigma = 3;
    const auto a = simulateAdsb(frames,local,config), b = simulateAdsb(frames,local,config);
    near(a.streams.at(id)[1].track.position,b.streams.at(id)[1].track.position,0,"Fixed seed is repeatable");
    config.coverage = 0;
    require(simulateAdsb(frames,local,config).streams.empty(), "CV-only scenario");
    config.coverage = 1;
    config.dropout = 1;
    require(simulateAdsb(frames,local,config).messages == 0, "Complete packet loss");
    config.period = 0;
    rejects([&] { simulateAdsb(frames,local,config); }, "Invalid simulation period rejected");
    config.period = 1e-300;
    rejects([&] { simulateAdsb(frames,local,config); }, "Unbounded message generation rejected");

    first.obstacles.clear();
    for (int i=0; i<10; ++i) first.obstacles.push_back({i,0,{},2,{0,0,17.0+i*20}});
    last = first;
    last.timestamp = 20;
    config.period = 1;
    config.coverage = .8;
    config.dropout = 0;
    const auto complete = simulateAdsb({first,last},local,config);
    config.dropout = .35;
    const auto partial = simulateAdsb({first,last},local,config);
    require(complete.reference_ids == partial.reference_ids,"Changing loss cannot change object coverage");
    for (const auto& stream : partial.streams) {
        for (const auto& message : stream.second) {
            const auto& original = complete.streams.at(stream.first);
            const auto found = std::find_if(original.begin(),original.end(),[&](const AdsbMessage& candidate) {
                return candidate.track.timestamp == message.track.timestamp;
            });
            require(found != original.end(),"Dropped stream is a subset of scheduled messages");
            near(message.track.position,found->track.position,0,"Shared samples retain identical noise");
        }
    }
}
} // namespace

int main() {
    try {
        interpolation();
        association();
        optimalAgainstEnumeration();
        geometry();
        reader();
        synchronization();
        simulation();
        std::cout << "Passed " << checks << " checks\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Test failed: " << error.what() << '\n';
        return 1;
    }
}
