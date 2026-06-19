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
    auto x = matlab::sym("xherm");
    REQUIRE_THROWS_AS(matlab::assume(x, "hermitian"), std::runtime_error);
}

TEST_CASE("matlab::assume supports the prime property",
          "[15m][matlab][assumptions]") {
    auto x = matlab::sym("xprime");
    matlab::clearAssumptions(x);
    matlab::assume(x, "prime");
    x = matlab::refresh(x);
    auto props = matlab::assumptions(x);
    REQUIRE(contains(props, "prime"));
    REQUIRE(is_prime(x) == true);
    REQUIRE(is_integer(x) == true);   // prime ⇒ integer, positive
    REQUIRE(is_positive(x) == true);
    matlab::clearAssumptions(x);
}

TEST_CASE("matlab::assume supports the composite property",
          "[15m][matlab][assumptions]") {
    auto x = matlab::sym("xcomp");
    matlab::clearAssumptions(x);
    matlab::assume(x, "composite");
    x = matlab::refresh(x);
    auto props = matlab::assumptions(x);
    REQUIRE(contains(props, "composite"));
    REQUIRE(is_composite(x) == true);
    REQUIRE(is_integer(x) == true);   // composite ⇒ integer, positive, ¬prime
    REQUIRE(is_positive(x) == true);
    REQUIRE(is_prime(x) == false);
    matlab::clearAssumptions(x);
}

TEST_CASE("matlab::assume supports the irrational property",
          "[15m][matlab][assumptions]") {
    auto x = matlab::sym("xirr");
    matlab::clearAssumptions(x);
    matlab::assume(x, "irrational");
    x = matlab::refresh(x);
    auto props = matlab::assumptions(x);
    REQUIRE(contains(props, "irrational"));
    REQUIRE(is_irrational(x) == true);
    REQUIRE(is_real(x) == true);       // irrational ⇒ real ∧ ¬rational
    REQUIRE(is_rational(x) == false);
    REQUIRE(is_integer(x) == false);
    matlab::clearAssumptions(x);
}

TEST_CASE("matlab::assume supports extended_real / infinite",
          "[15m][matlab][assumptions]") {
    auto e = matlab::sym("xext");
    matlab::clearAssumptions(e);
    matlab::assume(e, "extended_real");
    e = matlab::refresh(e);
    REQUIRE(contains(matlab::assumptions(e), "extended_real"));
    REQUIRE(is_extended_real(e) == true);
    REQUIRE(is_imaginary(e) == false);   // extended_real ⇒ ¬imaginary
    matlab::clearAssumptions(e);

    auto f = matlab::sym("xinf");
    matlab::clearAssumptions(f);
    matlab::assume(f, "infinite");
    f = matlab::refresh(f);
    REQUIRE(contains(matlab::assumptions(f), "infinite"));
    REQUIRE(is_infinite(f) == true);
    REQUIRE(is_finite(f) == false);      // infinite ⟺ ¬finite
    REQUIRE(is_real(f) == false);
    matlab::clearAssumptions(f);
}

TEST_CASE("matlab::assume supports algebraic / transcendental",
          "[15m][matlab][assumptions]") {
    auto a = matlab::sym("xalg");
    matlab::clearAssumptions(a);
    matlab::assume(a, "algebraic");
    a = matlab::refresh(a);
    REQUIRE(contains(matlab::assumptions(a), "algebraic"));
    REQUIRE(is_algebraic(a) == true);
    REQUIRE(is_complex(a) == true);    // algebraic ⇒ complex, finite
    REQUIRE(is_transcendental(a) == false);
    matlab::clearAssumptions(a);

    auto t = matlab::sym("xtrans");
    matlab::clearAssumptions(t);
    matlab::assume(t, "transcendental");
    t = matlab::refresh(t);
    REQUIRE(contains(matlab::assumptions(t), "transcendental"));
    REQUIRE(is_transcendental(t) == true);
    REQUIRE(is_algebraic(t) == false);  // transcendental ⇒ ¬algebraic, ¬rational
    REQUIRE(is_rational(t) == false);
    matlab::clearAssumptions(t);
}
