#include <cmath>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/json.h"
#include "core/metrics.h"

using namespace dmw;
using Catch::Matchers::WithinAbs;

TEST_CASE("json: parses the 3DTeethLand landmark shape", "[json]") {
    const auto r = parse_json(R"({"version":"1.1","key":"A_lower","objects":[
        {"key":"uuid_0","class":"Cusp","coord":[-2.47,-19.08,-91.05]},
        {"key":"uuid_1","class":"Mesial","coord":[1e-3, 2, 3]}]})");
    REQUIRE(r.ok());
    const JsonValue* objs = r.value.find("objects");
    REQUIRE(objs);
    REQUIRE(objs->array.size() == 2);
    CHECK(objs->array[0].find("class")->string == "Cusp");
    CHECK(objs->array[0].find("coord")->array[2].number == -91.05);
    CHECK(objs->array[1].find("coord")->array[0].number == 1e-3);
    CHECK(r.value.find("missing") == nullptr);
}

TEST_CASE("json: literals, escapes, nesting, empty containers", "[json]") {
    const auto r = parse_json(R"( {"a":[true,false,null,[],{}],"s":"x\"y\\z\n","u":"A"} )");
    REQUIRE(r.ok());
    const auto& a = r.value.find("a")->array;
    CHECK(a[0].boolean);
    CHECK(a[2].type == JsonValue::Type::Null);
    CHECK(a[3].array.empty());
    CHECK(r.value.find("s")->string == "x\"y\\z\n");
    CHECK(r.value.find("u")->string == "A");
}

TEST_CASE("json: malformed input is rejected with an offset", "[json]") {
    for (const char* bad : {"{\"a\":}", "[1,2", "{\"a\" 1}", "tru", "\"open", "[1] x", "{\"a\":1,}"}) {
        const auto r = parse_json(bad);
        CHECK_FALSE(r.ok());
        CHECK(r.error.rfind("offset", 0) == 0);
    }
}

TEST_CASE("metrics: one-to-one greedy matching within tolerance", "[metrics]") {
    // Two detections near one ground-truth point: only the closer one matches.
    const std::vector<Vec3> gt{{0, 0, 0}, {10, 0, 0}};
    const std::vector<Vec3> det{{0.2, 0, 0}, {0.1, 0, 0}, {10.5, 0, 0}, {50, 0, 0}};
    const auto r = match_points(det, gt, 1.0);
    CHECK(r.true_positives == 2);
    CHECK_THAT(r.precision(), WithinAbs(0.5, 1e-12));
    CHECK_THAT(r.recall(), WithinAbs(1.0, 1e-12));
    CHECK_THAT(r.f1(), WithinAbs(2 * 0.5 / 1.5, 1e-12));
    REQUIRE(r.matched_distances.size() == 2);
    CHECK_THAT(r.matched_distances[0], WithinAbs(0.1, 1e-12));  // the closer detection won
    CHECK_THAT(r.matched_distances[1], WithinAbs(0.5, 1e-12));
}

TEST_CASE("metrics: pooling across scans is a micro-average", "[metrics]") {
    MatchResult total;
    total.add(match_points(std::vector<Vec3>{{0, 0, 0}}, std::vector<Vec3>{{0, 0, 0}, {5, 0, 0}}, 1.0));
    total.add(match_points(std::vector<Vec3>{{9, 9, 9}}, std::vector<Vec3>{}, 1.0));
    CHECK(total.true_positives == 1);
    CHECK(total.detections == 2);
    CHECK(total.ground_truth == 2);
    CHECK_THAT(total.precision(), WithinAbs(0.5, 1e-12));
}
