// Pipeline smoke test: proves the build + test toolchain works end to end.
// Contains no geometry; it is replaced by real tests in Milestone 2.
#include <catch2/catch_test_macros.hpp>

#include "core/arith.h"

TEST_CASE("add returns the sum of two integers", "[smoke]") {
    REQUIRE(dmw::add(2, 3) == 5);
    REQUIRE(dmw::add(-4, 4) == 0);
    REQUIRE(dmw::add(0, 0) == 0);
}
