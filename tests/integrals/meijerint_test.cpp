// Meijer-G Mellin-transform / definite integration (meijerint Phase 3),
// cross-checked against SymPy.

#include <optional>

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/hypergeometric.hpp>
#include <sympp/integrals/meijerint.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;

namespace {
auto& oracle() { return sympp::testing::Oracle::instance(); }
}  // namespace

TEST_CASE("meijerint: ∫₀^∞ G dx via the Mellin master formula", "[meijerint][oracle]") {
    auto x = symbol("x");
    auto half = rational(1, 2);

    // G^{1,0}_{0,1}([],[],[b],[], x) = xᵇ·e^{−x}; ∫₀^∞ = Γ(b+1).
    auto I0 = meijerg_integrate_0_inf(meijerg({}, {}, {integer(0)}, {}, x), x);
    REQUIRE(I0.has_value());
    REQUIRE(oracle().equivalent((*I0)->str(), "1"));  // ∫ e^{−x} = 1

    auto Ih = meijerg_integrate_0_inf(meijerg({}, {}, {half}, {}, x), x);
    REQUIRE(Ih.has_value());
    REQUIRE(oracle().equivalent((*Ih)->str(), "sqrt(pi)/2"));  // Γ(3/2)

    auto I1 = meijerg_integrate_0_inf(meijerg({}, {}, {integer(1)}, {}, x), x);
    REQUIRE(I1.has_value());
    REQUIRE(oracle().equivalent((*I1)->str(), "1"));  // ∫ x·e^{−x} = 1

    // Scaled argument η = 2: ∫₀^∞ e^{−2x} dx = 1/2.
    auto Is = meijerg_integrate_0_inf(meijerg({}, {}, {integer(0)}, {}, mul(integer(2), x)), x);
    REQUIRE(Is.has_value());
    REQUIRE(oracle().equivalent((*Is)->str(), "1/2"));
}

TEST_CASE("meijerint: Mellin transform matches SymPy", "[meijerint][oracle]") {
    auto x = symbol("x"), s = symbol("s");
    auto half = rational(1, 2);

    // M[x^{1/2}·e^{−x}](s) = Γ(s + 1/2).
    auto mt = meijerg_mellin_transform(meijerg({}, {}, {half}, {}, x), x, s);
    REQUIRE(mt.has_value());
    REQUIRE(oracle().equivalent((*mt)->str(), "gamma(s + 1/2)"));

    // M[1/(1+x)](s) = Γ(s)·Γ(1−s) = π/sin(π·s).
    auto mt2 = meijerg_mellin_transform(meijerg({integer(0)}, {}, {integer(0)}, {}, x), x, s);
    REQUIRE(mt2.has_value());
    REQUIRE(oracle().equivalent((*mt2)->str(), "pi/sin(pi*s)"));
}

TEST_CASE("meijerint: divergent and non-Meijer inputs are rejected", "[meijerint]") {
    auto x = symbol("x");
    // ∫₀^∞ 1/(1+x) dx diverges → no finite value (Gamma pole in the formula).
    auto div = meijerg_integrate_0_inf(meijerg({integer(0)}, {}, {integer(0)}, {}, x), x);
    REQUIRE_FALSE(div.has_value());

    // A non-Meijer-G argument is rejected.
    REQUIRE_FALSE(meijerg_integrate_0_inf(x * x, x).has_value());
    // Argument not of the form η·x (here x²) → rejected.
    REQUIRE_FALSE(
        meijerg_integrate_0_inf(meijerg({}, {}, {integer(0)}, {}, pow(x, integer(2))), x)
            .has_value());
}
