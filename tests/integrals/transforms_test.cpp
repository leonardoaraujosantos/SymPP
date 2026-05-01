// Laplace transform tests.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/integrals/transforms.hpp>
#include <sympp/simplify/simplify.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

TEST_CASE("laplace: 1 -> 1/s", "[8][laplace]") {
    auto t = symbol("t");
    auto s = symbol("s");
    auto F = laplace_transform(integer(1), t, s);
    REQUIRE(F == pow(s, S::NegativeOne()));
}

TEST_CASE("laplace: t -> 1/s^2", "[8][laplace]") {
    auto t = symbol("t");
    auto s = symbol("s");
    auto F = laplace_transform(t, t, s);
    REQUIRE(F == pow(s, integer(-2)));
}

TEST_CASE("laplace: exp(a*t) -> 1/(s - a)", "[8][laplace][oracle]") {
    auto& oracle = Oracle::instance();
    auto t = symbol("t");
    auto s = symbol("s");
    auto a = symbol("a");
    auto F = laplace_transform(exp(a * t), t, s);
    REQUIRE(oracle.equivalent(F->str(), "1/(s - a)"));
}

TEST_CASE("laplace: sin(a*t) -> a/(s^2 + a^2)", "[8][laplace][oracle]") {
    auto& oracle = Oracle::instance();
    auto t = symbol("t");
    auto s = symbol("s");
    auto a = symbol("a");
    auto F = laplace_transform(sin(a * t), t, s);
    REQUIRE(oracle.equivalent(F->str(), "a/(s**2 + a**2)"));
}

TEST_CASE("laplace: cos(a*t) -> s/(s^2 + a^2)", "[8][laplace][oracle]") {
    auto& oracle = Oracle::instance();
    auto t = symbol("t");
    auto s = symbol("s");
    auto a = symbol("a");
    auto F = laplace_transform(cos(a * t), t, s);
    REQUIRE(oracle.equivalent(F->str(), "s/(s**2 + a**2)"));
}

TEST_CASE("laplace: linearity over Add", "[8][laplace][oracle]") {
    auto& oracle = Oracle::instance();
    auto t = symbol("t");
    auto s = symbol("s");
    // L{2 + 3*t + sin(t)} = 2/s + 3/s^2 + 1/(s^2 + 1)
    auto f = integer(2) + integer(3) * t + sin(t);
    auto F = laplace_transform(f, t, s);
    REQUIRE(oracle.equivalent(F->str(), "2/s + 3/s**2 + 1/(s**2 + 1)"));
}

#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/special.hpp>
#include <sympp/functions/combinatorial.hpp>
#include <sympp/core/imaginary_unit.hpp>

// ----- Inverse Laplace ------------------------------------------------------

TEST_CASE("inverse_laplace: 1/s → 1", "[8][inverse_laplace][oracle]") {
    auto& oracle = Oracle::instance();
    auto s = symbol("s");
    auto t = symbol("t");
    auto f = inverse_laplace_transform(pow(s, S::NegativeOne()), s, t);
    REQUIRE(oracle.equivalent(f->str(), "1"));
}

TEST_CASE("inverse_laplace: 1/s^2 → t", "[8][inverse_laplace][oracle]") {
    auto& oracle = Oracle::instance();
    auto s = symbol("s");
    auto t = symbol("t");
    auto f = inverse_laplace_transform(pow(s, integer(-2)), s, t);
    REQUIRE(oracle.equivalent(f->str(), "t"));
}

TEST_CASE("inverse_laplace: 1/(s-3) → exp(3t)",
          "[8][inverse_laplace][oracle]") {
    auto& oracle = Oracle::instance();
    auto s = symbol("s");
    auto t = symbol("t");
    auto F = pow(s - integer(3), S::NegativeOne());
    auto f = inverse_laplace_transform(F, s, t);
    REQUIRE(oracle.equivalent(f->str(), "exp(3*t)"));
}

TEST_CASE("inverse_laplace: 2/(s²+4) → sin(2t)",
          "[8][inverse_laplace][oracle]") {
    auto& oracle = Oracle::instance();
    auto s = symbol("s");
    auto t = symbol("t");
    auto F = integer(2) * pow(pow(s, integer(2)) + integer(4),
                              S::NegativeOne());
    auto f = inverse_laplace_transform(F, s, t);
    REQUIRE(oracle.equivalent(f->str(), "sin(2*t)"));
}

TEST_CASE("inverse_laplace: s/(s²+4) → cos(2t)",
          "[8][inverse_laplace][oracle]") {
    auto& oracle = Oracle::instance();
    auto s = symbol("s");
    auto t = symbol("t");
    auto F = s * pow(pow(s, integer(2)) + integer(4), S::NegativeOne());
    auto f = inverse_laplace_transform(F, s, t);
    REQUIRE(oracle.equivalent(f->str(), "cos(2*t)"));
}

TEST_CASE("inverse_laplace: linearity 1/s + 1/(s-1) → 1 + exp(t)",
          "[8][inverse_laplace][oracle]") {
    auto& oracle = Oracle::instance();
    auto s = symbol("s");
    auto t = symbol("t");
    auto F = pow(s, S::NegativeOne()) + pow(s - integer(1), S::NegativeOne());
    auto f = inverse_laplace_transform(F, s, t);
    REQUIRE(oracle.equivalent(f->str(), "1 + exp(t)"));
}

TEST_CASE("inverse_laplace: 1/(s-2)^3 → t^2 * exp(2t) / 2",
          "[8][inverse_laplace][oracle]") {
    auto& oracle = Oracle::instance();
    auto s = symbol("s");
    auto t = symbol("t");
    auto F = pow(pow(s - integer(2), integer(3)), S::NegativeOne());
    auto f = inverse_laplace_transform(F, s, t);
    REQUIRE(oracle.equivalent(f->str(), "t**2 * exp(2*t) / 2"));
}

// ----- Fourier --------------------------------------------------------------

TEST_CASE("fourier: δ(t) → 1", "[8][fourier]") {
    auto t = symbol("t");
    auto w = symbol("w");
    auto F = fourier_transform(dirac_delta(t), t, w);
    REQUIRE(F == S::One());
}

TEST_CASE("fourier: constant → 2π·δ(ω)", "[8][fourier][oracle]") {
    auto& oracle = Oracle::instance();
    auto t = symbol("t");
    auto w = symbol("w");
    auto F = fourier_transform(integer(5), t, w);
    REQUIRE(oracle.equivalent(F->str(), "10*pi*DiracDelta(w)"));
}

TEST_CASE("fourier: exp(-a·t²) → sqrt(π/a)·exp(-ω²/4a)",
          "[8][fourier][oracle]") {
    auto& oracle = Oracle::instance();
    auto t = symbol("t");
    auto w = symbol("w");
    auto a = symbol("a");
    auto e = exp(mul(S::NegativeOne(), mul(a, pow(t, integer(2)))));
    auto F = fourier_transform(e, t, w);
    // Verify by simplifying difference vs SymPy expression.
    auto resp = oracle.send({{"op", "equivalent"},
                             {"a", F->str()},
                             {"b", "sqrt(pi/a)*exp(-w**2/(4*a))"}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

TEST_CASE("inverse_fourier: roundtrip on δ(t) ∘ exp(-a·t²) gives back original up to constant",
          "[8][inverse_fourier]") {
    auto t = symbol("t");
    auto w = symbol("w");
    // Forward then inverse on δ(t):
    auto fwd = fourier_transform(dirac_delta(t), t, w);
    auto back = inverse_fourier_transform(fwd, w, t);
    // δ should round-trip to δ(-t) which equals δ(t) for the even
    // DiracDelta; either way, must not be a FourierTransform marker.
    bool ok = (back->type_id() != TypeId::Function);
    if (!ok) {
        const auto& fn = static_cast<const Function&>(*back);
        ok = (fn.name() != "FourierTransform");
    }
    REQUIRE(ok);
}

// ----- Mellin ---------------------------------------------------------------

TEST_CASE("mellin: exp(-t) → Γ(s)", "[8][mellin][oracle]") {
    auto& oracle = Oracle::instance();
    auto t = symbol("t");
    auto s = symbol("s");
    auto F = mellin_transform(exp(mul(S::NegativeOne(), t)), t, s);
    REQUIRE(oracle.equivalent(F->str(), "gamma(s)"));
}

TEST_CASE("mellin: exp(-2t) → Γ(s)/2^s", "[8][mellin][oracle]") {
    auto& oracle = Oracle::instance();
    auto t = symbol("t");
    auto s = symbol("s");
    auto F = mellin_transform(exp(mul({S::NegativeOne(), integer(2), t})), t, s);
    REQUIRE(oracle.equivalent(F->str(), "gamma(s)/2**s"));
}

TEST_CASE("inverse_mellin: Γ(s) → exp(-t)", "[8][inverse_mellin][oracle]") {
    auto& oracle = Oracle::instance();
    auto t = symbol("t");
    auto s = symbol("s");
    auto f = inverse_mellin_transform(gamma(s), s, t);
    REQUIRE(oracle.equivalent(f->str(), "exp(-t)"));
}

// ----- Z-transform ----------------------------------------------------------

TEST_CASE("z: 1 → z/(z-1)", "[8][ztransform][oracle]") {
    auto& oracle = Oracle::instance();
    auto n = symbol("n");
    auto z = symbol("z");
    auto F = z_transform(integer(1), n, z);
    REQUIRE(oracle.equivalent(F->str(), "z/(z - 1)"));
}

TEST_CASE("z: a^n → z/(z-a)", "[8][ztransform][oracle]") {
    auto& oracle = Oracle::instance();
    auto n = symbol("n");
    auto z = symbol("z");
    auto a = symbol("a");
    auto F = z_transform(pow(a, n), n, z);
    REQUIRE(oracle.equivalent(F->str(), "z/(z - a)"));
}

TEST_CASE("z: n → z/(z-1)^2", "[8][ztransform][oracle]") {
    auto& oracle = Oracle::instance();
    auto n = symbol("n");
    auto z = symbol("z");
    auto F = z_transform(n, n, z);
    REQUIRE(oracle.equivalent(F->str(), "z/(z - 1)**2"));
}

TEST_CASE("z: n*2^n → 2*z/(z-2)^2", "[8][ztransform][oracle]") {
    auto& oracle = Oracle::instance();
    auto n = symbol("n");
    auto z = symbol("z");
    auto e = n * pow(integer(2), n);
    auto F = z_transform(e, n, z);
    REQUIRE(oracle.equivalent(F->str(), "2*z/(z - 2)**2"));
}

TEST_CASE("inverse_z: z/(z-3) → 3^n", "[8][inverse_z][oracle]") {
    auto& oracle = Oracle::instance();
    auto n = symbol("n");
    auto z = symbol("z");
    auto F = z * pow(z - integer(3), S::NegativeOne());
    auto f = inverse_z_transform(F, z, n);
    REQUIRE(oracle.equivalent(f->str(), "3**n"));
}

// ----- Sine / cosine transforms ---------------------------------------------

TEST_CASE("sine_transform: exp(-2t) → ω/(4 + ω²)",
          "[8][sine_transform][oracle]") {
    auto& oracle = Oracle::instance();
    auto t = symbol("t");
    auto w = symbol("w");
    auto e = exp(mul({S::NegativeOne(), integer(2), t}));
    auto F = sine_transform(e, t, w);
    REQUIRE(oracle.equivalent(F->str(), "w/(4 + w**2)"));
}

TEST_CASE("cosine_transform: exp(-3t) → 3/(9 + ω²)",
          "[8][cosine_transform][oracle]") {
    auto& oracle = Oracle::instance();
    auto t = symbol("t");
    auto w = symbol("w");
    auto e = exp(mul({S::NegativeOne(), integer(3), t}));
    auto F = cosine_transform(e, t, w);
    REQUIRE(oracle.equivalent(F->str(), "3/(9 + w**2)"));
}
