// Meijer-G Mellin-transform / definite integration (meijerint Phase 3),
// cross-checked against SymPy.

#include <optional>

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/hypergeometric.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/integrals/integrate.hpp>
#include <sympp/integrals/meijerint.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/simplify/hyperexpand.hpp>

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

// ----- Phase 4: function → Meijer-G recognition ------------------------------

TEST_CASE("meijerint: to_meijerg round-trips through hyperexpand", "[meijerint][oracle]") {
    auto x = symbol("x");
    auto a = symbol("a");
    auto& o = oracle();
    auto roundtrip = [&](const Expr& f, const char* expected) {
        auto g = to_meijerg(f, x);
        REQUIRE(g.has_value());
        REQUIRE(o.equivalent(hyperexpand(*g)->str(), expected));
    };
    roundtrip(exp(mul(integer(-1), x)), "exp(-x)");          // e^{−x}
    roundtrip(exp(mul(integer(-2), x)), "exp(-2*x)");        // e^{−2x}
    roundtrip(mul(pow(x, integer(2)), exp(mul(integer(-1), x))),
              "x**2*exp(-x)");                                // x²·e^{−x}
    roundtrip(mul(pow(x, a), exp(mul(integer(-1), x))),
              "x**a*exp(-x)");                                // xᵃ·e^{−x}
    roundtrip(pow(integer(1) + x, integer(-1)), "1/(1 + x)");  // 1/(1+x)

    // Unrecognized → nullopt.
    REQUIRE_FALSE(to_meijerg(sin(x), x).has_value());
}

TEST_CASE("meijerint: ∫₀^∞ f dx via recognition equals SymPy", "[meijerint][oracle]") {
    auto x = symbol("x");
    auto a = symbol("a");
    auto& o = oracle();

    // ∫₀^∞ e^{−x} dx = 1.
    auto i1 = meijerg_integrate_0_inf_of(exp(mul(integer(-1), x)), x);
    REQUIRE(i1.has_value());
    REQUIRE(o.equivalent((*i1)->str(), "1"));

    // ∫₀^∞ x²·e^{−x} dx = Γ(3) = 2.
    auto i2 = meijerg_integrate_0_inf_of(mul(pow(x, integer(2)), exp(mul(integer(-1), x))), x);
    REQUIRE(i2.has_value());
    REQUIRE(o.equivalent((*i2)->str(), "2"));

    // ∫₀^∞ xᵃ·e^{−x} dx = Γ(a+1)  (the Gamma function from Meijer-G).
    auto ia = meijerg_integrate_0_inf_of(mul(pow(x, a), exp(mul(integer(-1), x))), x);
    REQUIRE(ia.has_value());
    REQUIRE(o.equivalent((*ia)->str(), "gamma(a + 1)"));
}

// ----- Wired into the general integrate(expr, var, 0, oo) ---------------------

TEST_CASE("integrate: ∫₀^∞ routes a Meijer-G through the Mellin formula",
          "[meijerint][oracle][integrate]") {
    auto x = symbol("x");
    auto a = symbol("a");
    auto oo = S::Infinity();
    auto& o = oracle();

    // A bare Meijer-G integrand — only reachable via the Meijer-G wiring.
    REQUIRE(o.equivalent(
        integrate(meijerg({}, {}, {integer(0)}, {}, x), x, S::Zero(), oo)->str(), "1"));
    REQUIRE(o.equivalent(
        integrate(meijerg({}, {}, {rational(1, 2)}, {}, x), x, S::Zero(), oo)->str(),
        "sqrt(pi)/2"));
    // x²·e^{−x} and xᵃ·e^{−x} via the general integrate definite path.
    REQUIRE(o.equivalent(
        integrate(mul(pow(x, integer(2)), exp(mul(integer(-1), x))), x, S::Zero(), oo)->str(),
        "2"));
    REQUIRE(o.equivalent(
        integrate(mul(pow(x, a), exp(mul(integer(-1), x))), x, S::Zero(), oo)->str(),
        "gamma(a + 1)"));
}

// ----- Table growth: ∫₀^∞ G(η·xᶜ) via the substitution (Gaussian/Dirichlet/Fresnel)

TEST_CASE("integrate: classic ∫₀^∞ via Meijer-G η·xᶜ substitution",
          "[meijerint][oracle][integrate]") {
    auto x = symbol("x");
    auto oo = S::Infinity();
    auto& o = oracle();
    auto I = [&](const Expr& f) { return integrate(f, x, S::Zero(), oo)->str(); };

    // Gaussian: ∫₀^∞ e^{−x²} = √π/2,  ∫₀^∞ x·e^{−x²} = 1/2,  ∫₀^∞ x²·e^{−x²} = √π/4.
    REQUIRE(o.equivalent(I(exp(mul(integer(-1), pow(x, integer(2))))), "sqrt(pi)/2"));
    REQUIRE(o.equivalent(I(mul(x, exp(mul(integer(-1), pow(x, integer(2)))))), "1/2"));
    REQUIRE(o.equivalent(I(mul(pow(x, integer(2)), exp(mul(integer(-1), pow(x, integer(2)))))),
                         "sqrt(pi)/4"));
    // Dirichlet: ∫₀^∞ sin(x)/x = π/2.
    REQUIRE(o.equivalent(I(mul(sin(x), pow(x, integer(-1)))), "pi/2"));
    // Fresnel: ∫₀^∞ cos(x²) = ∫₀^∞ sin(x²) = √(2π)/4.
    REQUIRE(o.equivalent(I(cos(pow(x, integer(2)))), "sqrt(2)*sqrt(pi)/4"));
    REQUIRE(o.equivalent(I(sin(pow(x, integer(2)))), "sqrt(2)*sqrt(pi)/4"));
    // Scaled Gaussian: ∫₀^∞ e^{−2x²} = √(2π)/4 = sqrt(2*pi)/4.
    REQUIRE(o.equivalent(I(exp(mul(integer(-2), pow(x, integer(2))))), "sqrt(2)*sqrt(pi)/4"));
}

TEST_CASE("meijerint: divergent oscillatory integrals are not (wrongly) valued",
          "[meijerint]") {
    auto x = symbol("x");
    // ∫₀^∞ sin(x) and ∫₀^∞ cos(x) diverge (oscillate) — the Meijer-G route must
    // abstain (strip condition fails), so the recognizer returns nullopt.
    REQUIRE_FALSE(meijerg_integrate_0_inf_of(sin(x), x).has_value());
    REQUIRE_FALSE(meijerg_integrate_0_inf_of(cos(x), x).has_value());
    // ∫₀^∞ x·sin(x) also diverges.
    REQUIRE_FALSE(meijerg_integrate_0_inf_of(mul(x, sin(x)), x).has_value());
}
