// MATLAB facade — assumptions tests.
//
// Validates the side-table assumption registry: `assume` /
// `assumeAlso` register MATLAB-style property strings, `refresh`
// returns a fresh Symbol carrying the registered mask, and downstream
// simplification picks up the assumption.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/matlab/matlab.hpp>

using namespace sympp;

namespace {

bool contains(const std::vector<std::string>& v, const std::string& s) {
    return std::find(v.begin(), v.end(), s) != v.end();
}

}  // namespace

TEST_CASE("matlab::assume positive then refresh yields positive symbol",
          "[15m][matlab][assumptions]") {
    auto x = matlab::sym("x_pos");
    matlab::clearAssumptions(x);
    matlab::assume(x, "positive");
    x = matlab::refresh(x);
    auto props = matlab::assumptions(x);
    REQUIRE(contains(props, "positive"));
    matlab::clearAssumptions(x);
}

TEST_CASE("matlab::assume positive: refreshed symbol is positive / real / nonzero",
          "[15m][matlab][assumptions]") {
    auto x = matlab::sym("xpa");
    matlab::clearAssumptions(x);
    matlab::assume(x, "positive");
    x = matlab::refresh(x);
    REQUIRE(is_positive(x) == true);
    REQUIRE(is_real(x) == true);
    REQUIRE(is_nonzero(x) == true);
    matlab::clearAssumptions(x);
}

TEST_CASE("matlab::assumeAlso merges into existing mask",
          "[15m][matlab][assumptions]") {
    auto x = matlab::sym("xn");
    matlab::clearAssumptions(x);
    matlab::assume(x, "integer");
    matlab::assumeAlso(x, "positive");
    auto props = matlab::assumptions(x);
    REQUIRE(contains(props, "integer"));
    REQUIRE(contains(props, "positive"));
    matlab::clearAssumptions(x);
}

TEST_CASE("matlab::assume replaces (does not merge)",
          "[15m][matlab][assumptions]") {
    auto x = matlab::sym("xr");
    matlab::clearAssumptions(x);
    matlab::assume(x, "integer");
    matlab::assume(x, "real");  // replace, not merge
    auto props = matlab::assumptions(x);
    REQUIRE(contains(props, "real"));
    REQUIRE_FALSE(contains(props, "integer"));
    matlab::clearAssumptions(x);
}

TEST_CASE("matlab::clearAssumptions wipes the registered mask",
          "[15m][matlab][assumptions]") {
    auto x = matlab::sym("xclr");
    matlab::assume(x, "positive");
    matlab::clearAssumptions(x);
    REQUIRE(matlab::assumptions(x).empty());
}

TEST_CASE("matlab::refresh on unregistered symbol returns input unchanged",
          "[15m][matlab][assumptions]") {
    auto x = matlab::sym("never_registered_xyz");
    matlab::clearAssumptions(x);
    auto y = matlab::refresh(x);
    REQUIRE(y == x);
}

TEST_CASE("matlab::assume rejects unsupported property",
          "[15m][matlab][assumptions]") {
    auto x = matlab::sym("xprime");
    REQUIRE_THROWS_AS(matlab::assume(x, "prime"), std::runtime_error);
}
