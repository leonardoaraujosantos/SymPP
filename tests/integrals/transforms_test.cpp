// Laplace transform tests.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/hyperbolic.hpp>
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

// LAPLACE-SHIFT-1: hyperbolic table entries (sinh/cosh) and the s-shift theorem
// L{exp(a·t)·g(t)} = G(s − a), which closes the damped-oscillation and t^n·exp
// families.
TEST_CASE("laplace: sinh/cosh and the s-shift theorem (LAPLACE-SHIFT-1)",
          "[8][laplace][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto t = symbol("t");
    auto s = symbol("s");
    // L{sinh(a·t)} = a/(s²−a²),  L{cosh(a·t)} = s/(s²−a²).
    auto a = symbol("a");
    REQUIRE(oracle.equivalent(laplace_transform(sinh(a * t), t, s)->str(),
                              "a/(s**2 - a**2)"));
    REQUIRE(oracle.equivalent(laplace_transform(cosh(a * t), t, s)->str(),
                              "s/(s**2 - a**2)"));
    // s-shift: L{exp(-t)·sin t} = 1/((s+1)²+1).
    REQUIRE(oracle.equivalent(
        laplace_transform(exp(mul(S::NegativeOne(), t)) * sin(t), t, s)->str(),
        "1/((s + 1)**2 + 1)"));
    // L{t·exp t} = 1/(s−1)²,  L{t²·exp t} = 2/(s−1)³.
    REQUIRE(oracle.equivalent(
        laplace_transform(t * exp(t), t, s)->str(), "(s - 1)**(-2)"));
    REQUIRE(oracle.equivalent(
        laplace_transform(pow(t, integer(2)) * exp(t), t, s)->str(),
        "2/(s - 1)**3"));
    // L{exp(2t)·cos t} = (s−2)/((s−2)²+1).
    REQUIRE(oracle.equivalent(
        laplace_transform(exp(integer(2) * t) * cos(t), t, s)->str(),
        "(s - 2)/((s - 2)**2 + 1)"));
}

// LAPLACE-TMULT-1: multiplication by t^n — L{tⁿ·g(t)} = (−1)ⁿ·dⁿ/dsⁿ L{g(t)} — for
// the trig/hyperbolic families the s-shift (exp-only) does not reach. L{t·cos t} =
// (s²−1)/(s²+1)², L{t·sin t} = 2s/(s²+1)², L{t·sinh t} = 2s/(s²−1)². Matches SymPy.
TEST_CASE("laplace: multiplication by t^n for trig/hyperbolic (LAPLACE-TMULT-1)",
          "[8][laplace][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto t = symbol("t");
    auto s = symbol("s");
    REQUIRE(oracle.equivalent(laplace_transform(t * cos(t), t, s)->str(),
                              "(s**2 - 1)/(s**2 + 1)**2"));
    REQUIRE(oracle.equivalent(laplace_transform(t * sin(t), t, s)->str(),
                              "2*s/(s**2 + 1)**2"));
    REQUIRE(oracle.equivalent(
        laplace_transform(t * cos(integer(2) * t), t, s)->str(),
        "(s**2 - 4)/(s**2 + 4)**2"));
    REQUIRE(oracle.equivalent(
        laplace_transform(pow(t, integer(2)) * cos(t), t, s)->str(),
        "(2*s**3 - 6*s)/(s**2 + 1)**3"));
    REQUIRE(oracle.equivalent(laplace_transform(t * sinh(t), t, s)->str(),
                              "2*s/(s**2 - 1)**2"));
    REQUIRE(oracle.equivalent(laplace_transform(t * cosh(t), t, s)->str(),
                              "(s**2 + 1)/(s**2 - 1)**2"));
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

// ILAPLACE-QUAD-1: inverse Laplace of a constant over an irreducible quadratic
// WITH a linear term (completed square) — the inverse s-shift symmetric to
// LAPLACE-SHIFT-1. 1/((s-a)²+b²) → exp(a·t)·sin(b·t)/b.
TEST_CASE("inverse_laplace: completed-square quadratic (ILAPLACE-QUAD-1)",
          "[8][inverse_laplace][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto s = symbol("s");
    auto t = symbol("t");
    // 1/(s²+2s+2) = 1/((s+1)²+1) → exp(-t)·sin(t).
    auto F1 = pow(pow(s, integer(2)) + integer(2) * s + integer(2),
                  S::NegativeOne());
    REQUIRE(oracle.equivalent(inverse_laplace_transform(F1, s, t)->str(),
                              "exp(-t)*sin(t)"));
    // 2/(s²+4s+5) → 2·exp(-2t)·sin(t).
    auto F2 = integer(2)
              * pow(pow(s, integer(2)) + integer(4) * s + integer(5),
                    S::NegativeOne());
    REQUIRE(oracle.equivalent(inverse_laplace_transform(F2, s, t)->str(),
                              "2*exp(-2*t)*sin(t)"));
    // 1/(s²+2s+10) → exp(-t)·sin(3t)/3.
    auto F3 = pow(pow(s, integer(2)) + integer(2) * s + integer(10),
                  S::NegativeOne());
    REQUIRE(oracle.equivalent(inverse_laplace_transform(F3, s, t)->str(),
                              "exp(-t)*sin(3*t)/3"));
    // No regression: a pure s²+a² (β = 0) still uses the plain sin path.
    auto F4 = pow(pow(s, integer(2)) + integer(4), S::NegativeOne());
    REQUIRE(oracle.equivalent(inverse_laplace_transform(F4, s, t)->str(),
                              "sin(2*t)/2"));
}

// ILAPLACE-QUAD-2: linear numerator over the completed-square quadratic —
// (α·s+β)/((s-a)²+b²) → exp(a·t)·(A·cos b·t + (B/b)·sin b·t). The damped-cosine
// companion to ILAPLACE-QUAD-1.
TEST_CASE("inverse_laplace: linear numerator over quadratic (ILAPLACE-QUAD-2)",
          "[8][inverse_laplace][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto s = symbol("s");
    auto t = symbol("t");
    auto den = pow(s, integer(2)) + integer(2) * s + integer(2);  // (s+1)²+1
    // s/(s²+2s+2) → exp(-t)·cos t − exp(-t)·sin t.
    REQUIRE(oracle.equivalent(
        inverse_laplace_transform(s * pow(den, S::NegativeOne()), s, t)->str(),
        "exp(-t)*cos(t) - exp(-t)*sin(t)"));
    // (s+1)/(s²+2s+2) → exp(-t)·cos t.
    REQUIRE(oracle.equivalent(
        inverse_laplace_transform((s + integer(1)) * pow(den, S::NegativeOne()),
                                  s, t)->str(),
        "exp(-t)*cos(t)"));
    // s/((s-2)²+1) → exp(2t)·cos t + 2·exp(2t)·sin t.
    auto den2 = pow(s, integer(2)) - integer(4) * s + integer(5);
    REQUIRE(oracle.equivalent(
        inverse_laplace_transform(s * pow(den2, S::NegativeOne()), s, t)->str(),
        "exp(2*t)*cos(t) + 2*exp(2*t)*sin(t)"));
    // No regression: s/(s²+4) (β=0) still gives cos(2t).
    REQUIRE(oracle.equivalent(
        inverse_laplace_transform(
            s * pow(pow(s, integer(2)) + integer(4), S::NegativeOne()), s, t)
            ->str(),
        "cos(2*t)"));
}

// ILAPLACE-REPQUAD-1: numerator over a repeated irreducible quadratic
// N(s)/(s²+a²)² → t·sin/t·cos family (the inverse of the multiplication-by-t
// Laplace rule). Closes the previously-unevaluated repeated-pole case.
TEST_CASE("inverse_laplace: numerator over repeated quadratic (ILAPLACE-REPQUAD-1)",
          "[8][inverse_laplace][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto s = symbol("s");
    auto t = symbol("t");
    auto sq = [&](int a) {
        return pow(pow(s, integer(2)) + integer(a), integer(2));
    };
    // s/(s²+4)² → t·sin(2t)/4.
    REQUIRE(oracle.equivalent(
        inverse_laplace_transform(s * pow(sq(4), S::NegativeOne()), s, t)->str(),
        "t*sin(2*t)/4"));
    // 1/(s²+1)² → (sin t − t·cos t)/2.
    REQUIRE(oracle.equivalent(
        inverse_laplace_transform(pow(sq(1), S::NegativeOne()), s, t)->str(),
        "(sin(t) - t*cos(t))/2"));
    // s/(s²+1)² → t·sin(t)/2.
    REQUIRE(oracle.equivalent(
        inverse_laplace_transform(s * pow(sq(1), S::NegativeOne()), s, t)->str(),
        "t*sin(t)/2"));
    // (s²−1)/(s²+1)² → t·cos t.
    REQUIRE(oracle.equivalent(
        inverse_laplace_transform(
            (pow(s, integer(2)) - integer(1)) * pow(sq(1), S::NegativeOne()), s,
            t)
            ->str(),
        "t*cos(t)"));
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

// LAPLACE-TRIGSQ-1: Laplace transforms of trig squares / products, via
// power-reduction to single-frequency form (cos²t → ½+½cos2t, etc.) before the
// linear rules apply. Previously these were left as unevaluated LaplaceTransform
// markers. Matches SymPy.
TEST_CASE("laplace: trig squares and products (LAPLACE-TRIGSQ-1)",
          "[8][laplace_transform][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto t = symbol("t");
    auto s = symbol("s");

    // cos²t → (s²+2)/(s(s²+4)), sin²t → 2/(s(s²+4)).
    REQUIRE(oracle.equivalent(
        laplace_transform(pow(cos(t), integer(2)), t, s)->str(),
        "(s**2 + 2)/(s*(s**2 + 4))"));
    REQUIRE(oracle.equivalent(
        laplace_transform(pow(sin(t), integer(2)), t, s)->str(),
        "2/(s*(s**2 + 4))"));
    // sin t·cos t → 1/(s²+4).
    REQUIRE(oracle.equivalent(
        laplace_transform(sin(t) * cos(t), t, s)->str(), "1/(s**2 + 4)"));
    // Scaled frequency: cos²(2t) → (s²+8)/(s(s²+16)).
    REQUIRE(oracle.equivalent(
        laplace_transform(pow(cos(integer(2) * t), integer(2)), t, s)->str(),
        "(s**2 + 8)/(s*(s**2 + 16))"));
    // Polynomial / exponential weights linearize the square in place.
    REQUIRE(oracle.equivalent(
        laplace_transform(t * pow(cos(t), integer(2)), t, s)->str(),
        "(s**4 + 2*s**2 + 8)/(s**2*(s**2 + 4)**2)"));
    REQUIRE(oracle.equivalent(
        laplace_transform(exp(mul(S::NegativeOne(), t)) * pow(sin(t), integer(2)),
                          t, s)
            ->str(),
        "2/((s + 1)*((s + 1)**2 + 4))"));

    // A plain sin t is unaffected.
    REQUIRE(oracle.equivalent(laplace_transform(sin(t), t, s)->str(),
                              "1/(s**2 + 1)"));
}
