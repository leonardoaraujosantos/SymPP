// integrate tests.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <string>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/functions/special.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/integrals/integrate.hpp>
#include <sympp/simplify/simplify.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- Polynomials ----------------------------------------------------------

TEST_CASE("integrate: constant", "[7][integrate]") {
    auto x = symbol("x");
    REQUIRE(integrate(integer(5), x) == integer(5) * x);
}

TEST_CASE("integrate: x", "[7][integrate]") {
    auto x = symbol("x");
    auto r = integrate(x, x);
    REQUIRE(r == mul(rational(1, 2), pow(x, integer(2))));
}

TEST_CASE("integrate: x^n", "[7][integrate][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(pow(x, integer(3)), x);
    // x^4 / 4
    REQUIRE(oracle.equivalent(r->str(), "x**4 / 4"));
}

TEST_CASE("integrate: 1/x = log(x)", "[7][integrate]") {
    auto x = symbol("x");
    auto r = integrate(pow(x, S::NegativeOne()), x);
    REQUIRE(r == log(x));
}

TEST_CASE("integrate: linearity", "[7][integrate][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // ∫ (3*x^2 + 2*x + 1) dx = x^3 + x^2 + x
    auto e = integer(3) * pow(x, integer(2)) + integer(2) * x + integer(1);
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(r->str(), "x**3 + x**2 + x"));
}

// ----- exp / sin / cos ------------------------------------------------------

TEST_CASE("integrate: exp(x)", "[7][integrate]") {
    auto x = symbol("x");
    REQUIRE(integrate(exp(x), x) == exp(x));
}

TEST_CASE("integrate: sin(x) = -cos(x)", "[7][integrate]") {
    auto x = symbol("x");
    REQUIRE(integrate(sin(x), x) == mul(S::NegativeOne(), cos(x)));
}

TEST_CASE("integrate: cos(x) = sin(x)", "[7][integrate]") {
    auto x = symbol("x");
    REQUIRE(integrate(cos(x), x) == sin(x));
}

TEST_CASE("integrate: exp(2*x) = exp(2*x) / 2",
          "[7][integrate][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(exp(integer(2) * x), x);
    REQUIRE(oracle.equivalent(r->str(), "exp(2*x) / 2"));
}

TEST_CASE("integrate: sin(3*x) = -cos(3*x)/3",
          "[7][integrate][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(sin(integer(3) * x), x);
    REQUIRE(oracle.equivalent(r->str(), "-cos(3*x) / 3"));
}

// ----- Definite integration -------------------------------------------------

TEST_CASE("integrate: definite x from 0 to 1 = 1/2",
          "[7][integrate]") {
    auto x = symbol("x");
    auto r = integrate(x, x, integer(0), integer(1));
    REQUIRE(r == rational(1, 2));
}

TEST_CASE("integrate: definite x^2 from 0 to 3 = 9",
          "[7][integrate]") {
    auto x = symbol("x");
    auto r = integrate(pow(x, integer(2)), x, integer(0), integer(3));
    REQUIRE(r == integer(9));
}

TEST_CASE("integrate: definite cos(x) from 0 to pi = 0",
          "[7][integrate]") {
    auto x = symbol("x");
    auto r = integrate(cos(x), x, S::Zero(), S::Pi());
    REQUIRE(r == integer(0));
}

// ----- Round-trip property: diff(integrate(f, x), x) == f --------------------

TEST_CASE("integrate / diff round-trip on polynomial",
          "[7][integrate][diff]") {
    auto x = symbol("x");
    auto f = pow(x, integer(3)) + integer(2) * x + integer(7);
    auto F = integrate(f, x);
    auto fp = diff(F, x);
    REQUIRE(simplify(fp - f) == S::Zero());
}

// ----- Hyperbolic integration (regression, issue #9 / INT-2) -----------------
// ∫sinh(ax+b) and ∫cosh(ax+b) used to fall through the table and return the
// unevaluated Integral(...) marker.
TEST_CASE("integrate: ∫sinh(x) dx = cosh(x)", "[7][integrate][hyperbolic][regression]") {
    auto x = symbol("x");
    REQUIRE(integrate(sinh(x), x) == cosh(x));
}

TEST_CASE("integrate: ∫cosh(x) dx = sinh(x)", "[7][integrate][hyperbolic][regression]") {
    auto x = symbol("x");
    REQUIRE(integrate(cosh(x), x) == sinh(x));
}

TEST_CASE("integrate: ∫sinh(2x+1) dx = cosh(2x+1)/2",
          "[7][integrate][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(sinh(integer(2) * x + integer(1)), x);
    REQUIRE(oracle.equivalent(r->str(), "cosh(2*x + 1)/2"));
}

TEST_CASE("integrate: ∫cosh(3x) dx = sinh(3x)/3",
          "[7][integrate][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(cosh(integer(3) * x), x);
    REQUIRE(oracle.equivalent(r->str(), "sinh(3*x)/3"));
}

// ----- Basic trig integration (regression, INT-3) ----------------------------
// ∫tan, ∫1/cos² (sec²) and ∫1/sin² (csc²) of an affine argument used to fall
// through the table and return the unevaluated Integral(...) marker. Forms are
// cross-checked against SymPy via the oracle (SymPP emits sin/cos rather than
// the tan/cot SymPy uses; the oracle confirms they are equivalent).
TEST_CASE("integrate: ∫tan(x) dx = -log(cos(x))",
          "[7][integrate][trig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(tan(x), x);
    REQUIRE(oracle.equivalent(r->str(), "-log(cos(x))"));
}

TEST_CASE("integrate: ∫1/cos(x)^2 dx = tan(x)",
          "[7][integrate][trig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(pow(cos(x), integer(-2)), x);
    REQUIRE(oracle.equivalent(r->str(), "tan(x)"));
}

TEST_CASE("integrate: ∫1/sin(x)^2 dx = -cot(x)",
          "[7][integrate][trig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(pow(sin(x), integer(-2)), x);
    REQUIRE(oracle.equivalent(r->str(), "-cot(x)"));
}

TEST_CASE("integrate: ∫tan(2x+1) dx = -log(cos(2x+1))/2",
          "[7][integrate][trig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(tan(integer(2) * x + integer(1)), x);
    REQUIRE(oracle.equivalent(r->str(), "-log(cos(2*x + 1))/2"));
}

TEST_CASE("integrate: ∫1/cos(3x)^2 dx = tan(3x)/3",
          "[7][integrate][trig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(pow(cos(integer(3) * x), integer(-2)), x);
    REQUIRE(oracle.equivalent(r->str(), "tan(3*x)/3"));
}

// ----- Polynomial × log integration by parts (regression, INT-4) -------------
// ∫xⁿ·log(ax+b) dx fell through: by-parts only used sin/cos/exp as the dv
// factor, never log. Now u = log(ax+b), dv = poly dx. Verified by the oracle
// and by differentiation (forms are valid up to an additive constant).
TEST_CASE("integrate: ∫x*log(x) dx = x^2*log(x)/2 - x^2/4",
          "[7][integrate][byparts][log][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(x * log(x), x);
    REQUIRE(r->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(r->str(), "x**2*log(x)/2 - x**2/4"));
}

TEST_CASE("integrate: ∫x^2*log(x) dx = x^3*log(x)/3 - x^3/9",
          "[7][integrate][byparts][log][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(pow(x, integer(2)) * log(x), x);
    REQUIRE(oracle.equivalent(r->str(), "x**3*log(x)/3 - x**3/9"));
}

TEST_CASE("integrate: ∫(x+1)*log(x) dx",
          "[7][integrate][byparts][log][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate((x + integer(1)) * log(x), x);
    REQUIRE(oracle.equivalent(r->str(),
                              "(x**2/2 + x)*log(x) - x**2/4 - x"));
}

TEST_CASE("integrate: ∫x*log(2x+1) dx (affine log argument)",
          "[7][integrate][byparts][log][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto f = x * log(integer(2) * x + integer(1));
    auto F = integrate(f, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    // SymPP's antiderivative differs from SymPy's by an additive constant
    // (log(x+1/2) vs log(2x+1)), so compare derivatives, not the primitives.
    REQUIRE(oracle.equivalent(diff(F, x)->str(), "x*log(2*x + 1)"));
}

// ----- Arctangent of an irreducible quadratic (regression, INT-5) ------------
// ∫1/(a·x²+b·x+c) dx with negative discriminant used to fall through (the
// rational path only decomposes denominators with real roots). It now returns
// the arctangent closed form. Reducible denominators (e.g. 1/(x²−1)) keep
// going through try_rational and are unaffected.
TEST_CASE("integrate: ∫1/(1+x^2) dx = atan(x)",
          "[7][integrate][arctan][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(pow(pow(x, integer(2)) + integer(1), integer(-1)), x);
    REQUIRE(r->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(r->str(), "atan(x)"));
}

TEST_CASE("integrate: ∫1/(x^2+4) dx = atan(x/2)/2",
          "[7][integrate][arctan][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(pow(pow(x, integer(2)) + integer(4), integer(-1)), x);
    REQUIRE(oracle.equivalent(r->str(), "atan(x/2)/2"));
}

TEST_CASE("integrate: ∫1/(4x^2+9) dx = atan(2x/3)/6",
          "[7][integrate][arctan][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto den = integer(4) * pow(x, integer(2)) + integer(9);
    auto r = integrate(pow(den, integer(-1)), x);
    REQUIRE(oracle.equivalent(r->str(), "atan(2*x/3)/6"));
}

TEST_CASE("integrate: ∫1/(x^2+2x+5) dx = atan((x+1)/2)/2 (completed square)",
          "[7][integrate][arctan][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto den = pow(x, integer(2)) + integer(2) * x + integer(5);
    auto r = integrate(pow(den, integer(-1)), x);
    REQUIRE(r->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(r->str(), "atan((x+1)/2)/2"));
}

TEST_CASE("integrate: ∫1/(x^2-1) dx stays a log (reducible, unaffected)",
          "[7][integrate][arctan][regression]") {
    auto x = symbol("x");
    auto r = integrate(pow(pow(x, integer(2)) - integer(1), integer(-1)), x);
    // Real roots -> partial fractions / logs, NOT atan.
    REQUIRE(r->str().find("atan") == std::string::npos);
    REQUIRE(r->str().find("log") != std::string::npos);
}

// ----- Partial fractions over irreducible factors (INT-17) -------------------
// apart() now decomposes over irreducible quadratics (squarefree), so rational
// integrals with such factors close: ∫1/(x³+1), ∫1/(x⁴-1).
TEST_CASE("integrate: ∫1/(x^3+1) via partial fractions (log + atan)",
          "[7][integrate][rational][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(pow(pow(x, integer(3)) + integer(1), integer(-1)), x);
    REQUIRE(r->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(
        r->str(),
        "log(x+1)/3 - log(x**2-x+1)/6 + sqrt(3)*atan(2*sqrt(3)*x/3 - sqrt(3)/3)/3"));
}

TEST_CASE("integrate: ∫1/(x^4-1) via partial fractions",
          "[7][integrate][rational][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(pow(pow(x, integer(4)) - integer(1), integer(-1)), x);
    REQUIRE(r->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(
        r->str(), "log(x-1)/4 - log(x+1)/4 - atan(x)/2"));
}

// ----- Linear over irreducible quadratic (INT-16) ----------------------------
// ∫(p·x+q)/(a·x²+b·x+c) splits into (p/2a)·log(quadratic) + remainder·atan.
TEST_CASE("integrate: ∫(linear)/(irreducible quadratic)",
          "[7][integrate][rational][arctan][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto sq = [&](Expr e) { return pow(e, integer(2)); };
    // Numerator ∝ d(denom): pure log, no atan term.
    REQUIRE(oracle.equivalent(
        integrate((x + integer(1)) / (sq(x) + integer(2) * x + integer(5)), x)
            ->str(),
        "log(x**2 + 2*x + 5)/2"));
    // log + atan.
    REQUIRE(oracle.equivalent(
        integrate((integer(2) * x + integer(3)) / (sq(x) + integer(1)), x)
            ->str(),
        "log(x**2 + 1) + 3*atan(x)"));
    REQUIRE(oracle.equivalent(
        integrate((integer(3) * x + integer(5)) / (sq(x) + integer(4)), x)
            ->str(),
        "3*log(x**2 + 4)/2 + 5*atan(x/2)/2"));
    // Completed-square denominator with a linear numerator.
    REQUIRE(oracle.equivalent(
        integrate(x / (sq(x) + integer(2) * x + integer(5)), x)->str(),
        "log(x**2 + 2*x + 5)/2 - atan(x/2 + 1/2)/2"));
    // Cleaned-up log (previously emitted the spurious 1/2*log(2*(x²+1))).
    REQUIRE(oracle.equivalent(integrate(x / (sq(x) + integer(1)), x)->str(),
                              "log(x**2 + 1)/2"));
}

// ----- Arc-sine / arc-sinh of a quadratic radicand (regression, INT-6) -------
// ∫1/√(a·x²+c) dx (pure quadratic radicand, c>0) used to fall through. It now
// returns asinh (a>0) or asin (a<0). The x²−1 case (acosh/log) stays deferred.
TEST_CASE("integrate: ∫1/sqrt(1-x^2) dx = asin(x)",
          "[7][integrate][invtrig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(pow(integer(1) - pow(x, integer(2)), rational(-1, 2)), x);
    REQUIRE(r->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(r->str(), "asin(x)"));
}

TEST_CASE("integrate: ∫1/sqrt(9-4x^2) dx = asin(2x/3)/2",
          "[7][integrate][invtrig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto rad = integer(9) - integer(4) * pow(x, integer(2));
    auto r = integrate(pow(rad, rational(-1, 2)), x);
    REQUIRE(oracle.equivalent(r->str(), "asin(2*x/3)/2"));
}

TEST_CASE("integrate: ∫1/sqrt(x^2+1) dx = asinh(x)",
          "[7][integrate][invtrig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(pow(pow(x, integer(2)) + integer(1), rational(-1, 2)), x);
    REQUIRE(r->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(r->str(), "asinh(x)"));
}

TEST_CASE("integrate: ∫1/sqrt(4x^2+9) dx = asinh(2x/3)/2",
          "[7][integrate][invtrig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto rad = integer(4) * pow(x, integer(2)) + integer(9);
    auto r = integrate(pow(rad, rational(-1, 2)), x);
    REQUIRE(oracle.equivalent(r->str(), "asinh(2*x/3)/2"));
}

// ----- Logarithmic form of 1/sqrt(a*x^2+c), c<0 (regression, INT-7) ----------
// The a>0, c<0 case: ∫1/√(a·x²+c) = log(√a·x + √(a·x²+c))/√a (acosh-equivalent,
// the form SymPy prints). Verified by differentiation (valid up to a constant).
TEST_CASE("integrate: ∫1/sqrt(x^2-1) dx = log(x + sqrt(x^2-1))",
          "[7][integrate][invtrig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto f = pow(pow(x, integer(2)) - integer(1), rational(-1, 2));
    auto F = integrate(f, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    // SymPP's simplify can't reduce the radical algebra to 0, so confirm the
    // derivative equals the integrand through the oracle instead.
    REQUIRE(oracle.equivalent(diff(F, x)->str(), "1/sqrt(x**2 - 1)"));
}

TEST_CASE("integrate: ∫1/sqrt(4x^2-9) dx = log(2x + sqrt(4x^2-9))/2",
          "[7][integrate][invtrig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto rad = integer(4) * pow(x, integer(2)) - integer(9);
    auto r = integrate(pow(rad, rational(-1, 2)), x);
    REQUIRE(oracle.equivalent(r->str(), "log(2*x + sqrt(4*x**2 - 9))/2"));
}

// ----- √(quadratic) in the numerator via by-parts (regression, INT-20) -------
// ∫√(a·x²+c) dx (pure quadratic radicand) used to fall through. By parts it
// reduces to (x/2)·√(a·x²+c) + (c/2)·∫1/√(a·x²+c), reusing the INT-6/7 result:
// asin (a<0), asinh (a>0,c>0) or log (a>0,c<0). Verified by differentiation,
// since simplify can't reduce the radical algebra of the antiderivative to 0.
TEST_CASE("integrate: ∫sqrt(1-x^2) dx = x*sqrt(1-x^2)/2 + asin(x)/2",
          "[7][integrate][invtrig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(integer(1) - pow(x, integer(2)), rational(1, 2));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫sqrt(4-x^2) dx = x*sqrt(4-x^2)/2 + 2*asin(x/2)",
          "[7][integrate][invtrig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(integer(4) - pow(x, integer(2)), rational(1, 2));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫sqrt(x^2+1) dx = x*sqrt(x^2+1)/2 + asinh(x)/2",
          "[7][integrate][invtrig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(pow(x, integer(2)) + integer(1), rational(1, 2));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫sqrt(2x^2+3) dx (asinh form, non-unit coefficients)",
          "[7][integrate][invtrig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(integer(2) * pow(x, integer(2)) + integer(3), rational(1, 2));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫sqrt(x^2-1) dx = x*sqrt(x^2-1)/2 - log(x+sqrt(x^2-1))/2",
          "[7][integrate][invtrig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(pow(x, integer(2)) - integer(1), rational(1, 2));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

// ----- Polynomial × (linear)^rational via u-sub (regression, INT-21) ---------
// ∫P(x)·(a·x+b)^r dx (r a non-integer rational) used to fall through. The
// substitution u = a·x+b turns it into Σ cₖ·u^(k+r), integrated term-by-term.
// SymPP can't reduce the radical antiderivative to SymPy's printed form, so
// each case is verified by differentiation against the oracle.
TEST_CASE("integrate: ∫x*sqrt(x+1) dx (algebraic u-sub)",
          "[7][integrate][algebraic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = x * pow(x + integer(1), rational(1, 2));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫x*sqrt(2x+3) dx (non-unit slope)",
          "[7][integrate][algebraic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = x * pow(integer(2) * x + integer(3), rational(1, 2));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫x^2*sqrt(x+1) dx (degree-2 polynomial factor)",
          "[7][integrate][algebraic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(x, integer(2)) * pow(x + integer(1), rational(1, 2));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫x/sqrt(x+1) dx (negative fractional exponent)",
          "[7][integrate][algebraic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = x * pow(x + integer(1), rational(-1, 2));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫x*(x+1)^(1/3) dx (cube root)",
          "[7][integrate][algebraic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = x * pow(x + integer(1), rational(1, 3));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

// ----- Rational function of exp(x) via u = exp(x) (regression, INT-22) -------
// ∫1/(1+exp(x)) and friends fell through: try_heurisch found g = exp(x) but its
// inner integration was table-only, so the substituted rational integrand
// 1/(u·(1+u)) wasn't decomposed. heurisch now hands a substituted-but-rational
// integrand to try_rational (partial fractions). Verified by differentiation.
TEST_CASE("integrate: ∫1/(1+exp(x)) dx = x - log(1+exp(x))",
          "[7][integrate][heurisch][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(integer(1) + exp(x), integer(-1));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫exp(x)/(1+exp(x)) dx = log(1+exp(x))",
          "[7][integrate][heurisch][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = exp(x) * pow(integer(1) + exp(x), integer(-1));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫1/(exp(x)-1) dx (rational in exp, pole at 0)",
          "[7][integrate][heurisch][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(exp(x) - integer(1), integer(-1));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫1/(1+exp(2x)) dx (rational in exp(2x))",
          "[7][integrate][heurisch][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(integer(1) + exp(integer(2) * x), integer(-1));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

// ----- Polynomial × exp × sin/cos via cyclic by-parts (regression, INT-23) ---
// ∫P(x)·e^(a·x)·sin/cos(g·x) fell through: the pure cyclic case bails on the
// polynomial factor, and the single-function by-parts bails because u = x·sin(x)
// isn't polynomial. By parts with u = P(x), dv = e·trig (the cyclic closed form)
// lowers deg(P) each step down to the bare cyclic base case.
TEST_CASE("integrate: ∫x*exp(x)*sin(x) dx (poly × exp × trig)",
          "[7][integrate][byparts][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = x * exp(x) * sin(x);
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫x*exp(x)*cos(x) dx (poly × exp × cos)",
          "[7][integrate][byparts][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = x * exp(x) * cos(x);
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫x^2*exp(x)*sin(x) dx (degree-2 poly, recurses twice)",
          "[7][integrate][byparts][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(x, integer(2)) * exp(x) * sin(x);
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫x*exp(2x)*sin(3x) dx (non-unit exp and trig rates)",
          "[7][integrate][byparts][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = x * exp(integer(2) * x) * sin(integer(3) * x);
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

// ----- cot / sec / csc antiderivatives (regression, INT-24) ------------------
// With cot/sec/csc now distinct function types (TRIG-RECIP), the integrate_term
// affine-function switch gains their table entries: ∫cot=log(sin),
// ∫sec=(log(sin+1)-log(sin-1))/2, ∫csc=(log(cos-1)-log(cos+1))/2, each scaled
// by 1/a for an affine argument. Verified by differentiation against the oracle.
TEST_CASE("integrate: ∫cot(x) dx = log(sin(x))",
          "[7][integrate][reciprocal][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = cot(x);
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫sec(x) dx (logarithmic form)",
          "[7][integrate][reciprocal][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = sec(x);
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫csc(x) dx (logarithmic form)",
          "[7][integrate][reciprocal][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = csc(x);
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫cot(2x+1) dx (affine argument scaling)",
          "[7][integrate][reciprocal][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = cot(integer(2) * x + integer(1));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫sec(3x) dx (affine argument scaling)",
          "[7][integrate][reciprocal][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = sec(integer(3) * x);
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

// ----- Reciprocal-trio powers (regression, INT-25) ---------------------------
// ∫sec²/csc²/cot² rewrite to the 1/cos², 1/sin² table cases (cot² via the
// Pythagorean cot²=csc²−1); ∫cotⁿ uses the reduction
// ∫cotⁿ=−cot^(n-1)/((n-1)a)−∫cot^(n-2), the sign-flipped analogue of ∫tanⁿ.
// Verified by differentiation against the oracle.
TEST_CASE("integrate: ∫sec(x)^2 dx = tan(x)",
          "[7][integrate][reciprocal][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(sec(x), integer(2));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫csc(x)^2 dx = -cot(x)",
          "[7][integrate][reciprocal][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(csc(x), integer(2));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫cot(x)^2 dx = -cot(x) - x (Pythagorean)",
          "[7][integrate][reciprocal][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(cot(x), integer(2));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫cot(x)^3 dx (reduction to ∫cot)",
          "[7][integrate][reciprocal][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(cot(x), integer(3));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫cot(x)^4 dx (reduction to ∫cot²)",
          "[7][integrate][reciprocal][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(cot(x), integer(4));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫sec(2x)^2 dx (affine argument scaling)",
          "[7][integrate][reciprocal][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(sec(integer(2) * x), integer(2));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

// ----- ∫secⁿ / ∫cscⁿ via by-parts reduction (regression, INT-27) -------------
// ∫secⁿ = sec^(n-2)·tan/((n-1)a) + (n-2)/(n-1)·∫sec^(n-2) (and the −cot analogue
// for csc), recursing to the ∫sec table case (INT-24) / ∫sec² square (INT-25).
// Verified by differentiation against the oracle.
TEST_CASE("integrate: ∫sec(x)^3 dx (by-parts reduction)",
          "[7][integrate][reciprocal][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(sec(x), integer(3));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫sec(x)^4 dx (reduction to ∫sec²)",
          "[7][integrate][reciprocal][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(sec(x), integer(4));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫csc(x)^3 dx (by-parts reduction)",
          "[7][integrate][reciprocal][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(csc(x), integer(3));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫csc(x)^4 dx (reduction to ∫csc²)",
          "[7][integrate][reciprocal][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(csc(x), integer(4));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫sec(2x)^3 dx (affine argument scaling)",
          "[7][integrate][reciprocal][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(sec(integer(2) * x), integer(3));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

// ----- Hyperbolic reciprocal antiderivatives & squares (regression, INT-26) --
// With coth/sech/csch now function types (HYP-RECIP): ∫coth=log(sinh),
// ∫sech=atan(sinh) (Gudermannian), ∫csch=log(tanh(u/2)); the squares rewrite to
// the 1/cosh², 1/sinh² table cases (coth²=csch²+1). Verified by differentiation.
TEST_CASE("integrate: ∫coth(x) dx = log(sinh(x))",
          "[7][integrate][reciprocal][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = coth(x);
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫sech(x) dx = atan(sinh(x))",
          "[7][integrate][reciprocal][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = sech(x);
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫csch(x) dx = log(tanh(x/2))",
          "[7][integrate][reciprocal][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = csch(x);
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫sech(x)^2 dx = tanh(x)",
          "[7][integrate][reciprocal][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(sech(x), integer(2));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫csch(x)^2 dx = -coth(x)",
          "[7][integrate][reciprocal][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(csch(x), integer(2));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫coth(x)^2 dx = x - coth(x) (Pythagorean)",
          "[7][integrate][reciprocal][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(coth(x), integer(2));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫coth(2x+1) dx (affine argument scaling)",
          "[7][integrate][reciprocal][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = coth(integer(2) * x + integer(1));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

// ----- ∫sechⁿ / ∫cschⁿ via by-parts reduction (regression, INT-28) -----------
// Hyperbolic analogue of INT-27. The Pythagorean sign differs (coth²−csch²=1),
// so the csch rest term is subtracted:
// ∫sechⁿ = sech^(n-2)·tanh/((n-1)a) + (n-2)/(n-1)·∫sech^(n-2);
// ∫cschⁿ = −csch^(n-2)·coth/((n-1)a) − (n-2)/(n-1)·∫csch^(n-2).
// Verified by differentiation against the oracle.
TEST_CASE("integrate: ∫sech(x)^3 dx (by-parts reduction)",
          "[7][integrate][reciprocal][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(sech(x), integer(3));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫sech(x)^4 dx (reduction to ∫sech²)",
          "[7][integrate][reciprocal][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(sech(x), integer(4));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫csch(x)^3 dx (by-parts reduction, sign-flipped rest)",
          "[7][integrate][reciprocal][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(csch(x), integer(3));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫csch(x)^4 dx (reduction to ∫csch²)",
          "[7][integrate][reciprocal][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(csch(x), integer(4));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫sech(2x)^3 dx (affine argument scaling)",
          "[7][integrate][reciprocal][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(sech(integer(2) * x), integer(3));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

// ----- tan² via the Pythagorean identity (regression, INT-8) -----------------
// ∫tan²(u) du fell through (only sin²/cos² had a trig-reduction rewrite). Now
// tan²(u) → 1/cos²(u) − 1, so ∫tan²(u) = tan(u)/a − u for an affine u. SymPP
// spells tan(u) as sin(u)/cos(u); the oracle confirms equivalence.
TEST_CASE("integrate: ∫tan(x)^2 dx = tan(x) - x",
          "[7][integrate][trig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(pow(tan(x), integer(2)), x);
    REQUIRE(r->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(r->str(), "tan(x) - x"));
}

TEST_CASE("integrate: ∫tan(2x+1)^2 dx = tan(2x+1)/2 - x",
          "[7][integrate][trig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(pow(tan(integer(2) * x + integer(1)), integer(2)), x);
    REQUIRE(oracle.equivalent(r->str(), "tan(2*x + 1)/2 - x"));
}

// ----- Repeated-root quadratic denominator (regression, INT-9) ---------------
// ∫1/(a·x²+b·x+c) with discriminant 0 = a perfect square a·(x−r)². try_rational
// (apart) didn't decompose the repeated root and the arctan branch needs D>0,
// so it fell through. Now returns −2/(2ax+b). Distinct-real-root and
// irreducible cases are unaffected.
TEST_CASE("integrate: ∫1/(x^2+2x+1) dx = -1/(x+1)",
          "[7][integrate][arctan][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto den = pow(x, integer(2)) + integer(2) * x + integer(1);
    auto r = integrate(pow(den, integer(-1)), x);
    REQUIRE(r->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(r->str(), "-1/(x + 1)"));
}

TEST_CASE("integrate: ∫1/(4x^2+4x+1) dx = -1/(2(2x+1))",
          "[7][integrate][arctan][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto den = integer(4) * pow(x, integer(2)) + integer(4) * x + integer(1);
    auto r = integrate(pow(den, integer(-1)), x);
    REQUIRE(oracle.equivalent(r->str(), "-1/(4*x + 2)"));
}

// ----- Heurisch u-substitution with a nested inner function (INT-10) ---------
// ∫1/(x·log(x)) = log(log(x)) via u = log(x). Heurisch only collected function
// *args* and Pow *bases* as candidates, so a function buried as a Mul factor
// (here log(x) inside the (x·log(x))⁻¹ base) was never tried. Now the function
// application itself is a candidate.
TEST_CASE("integrate: ∫1/(x*log(x)) dx = log(log(x))",
          "[7][integrate][heurisch][regression]") {
    auto x = symbol("x");
    auto f = pow(x * log(x), integer(-1));
    auto F = integrate(f, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(simplify(diff(F, x) - f) == S::Zero());
}

TEST_CASE("integrate: ∫1/(x*log(x)^2) dx = -1/log(x)",
          "[7][integrate][heurisch][regression]") {
    auto x = symbol("x");
    auto f = pow(x * pow(log(x), integer(2)), integer(-1));
    auto F = integrate(f, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(simplify(diff(F, x) - f) == S::Zero());
}

// ----- Gaussian integral → erf (regression, INT-11) -------------------------
// ∫exp(−a·x²) dx = √π·erf(√a·x)/(2√a) — the error-function antiderivative.
// Verified by differentiation (erf' was fixed in DIFF-2).
TEST_CASE("integrate: ∫exp(-x^2) dx = sqrt(pi)*erf(x)/2",
          "[7][integrate][erf][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto f = exp(mul(S::NegativeOne(), pow(x, integer(2))));
    auto F = integrate(f, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(F->str(), "sqrt(pi)*erf(x)/2"));
}

TEST_CASE("integrate: ∫2*exp(-x^2)/sqrt(pi) dx = erf(x)",
          "[7][integrate][erf][regression]") {
    auto x = symbol("x");
    auto f = exp(mul(S::NegativeOne(), pow(x, integer(2))));
    auto F = integrate(f, x);
    // d/dx of the result is the integrand.
    REQUIRE(simplify(diff(F, x) - f) == S::Zero());
}

// Positive-coefficient Gaussian → erfi (regression, ERFI): ∫exp(a·x²) dx =
// √π·erfi(√a·x)/(2√a). Verified by differentiation against the oracle.
TEST_CASE("integrate: ∫exp(x^2) dx = sqrt(pi)*erfi(x)/2",
          "[7][integrate][erf][erfi][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = exp(pow(x, integer(2)));
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(F->str(), "sqrt(pi)*erfi(x)/2"));
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

// ----- Hyperbolic table integrals (regression, INT-12) -----------------------
// Hyperbolic analogues of INT-3: ∫tanh = log(cosh), ∫1/cosh² = tanh,
// ∫1/sinh² = −cosh/sinh (= −coth), of an affine argument. Oracle-checked
// (SymPP spells −coth as −cosh/sinh).
TEST_CASE("integrate: ∫tanh(x) dx = log(cosh(x))",
          "[7][integrate][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(tanh(x), x);
    REQUIRE(r->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(r->str(), "log(cosh(x))"));
}

TEST_CASE("integrate: ∫1/cosh(x)^2 dx = tanh(x)",
          "[7][integrate][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(pow(cosh(x), integer(-2)), x);
    REQUIRE(oracle.equivalent(r->str(), "tanh(x)"));
}

TEST_CASE("integrate: ∫1/sinh(x)^2 dx = -cosh(x)/sinh(x)",
          "[7][integrate][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(pow(sinh(x), integer(-2)), x);
    REQUIRE(oracle.equivalent(r->str(), "-cosh(x)/sinh(x)"));
}

TEST_CASE("integrate: ∫tanh(2x+1) dx = log(cosh(2x+1))/2",
          "[7][integrate][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto r = integrate(tanh(integer(2) * x + integer(1)), x);
    REQUIRE(oracle.equivalent(r->str(), "log(cosh(2*x + 1))/2"));
}

// ----- Cyclic integration by parts (regression, issue #7 / INT-1) ------------
// ∫exp(x)*sin(x) dx used to recurse exp·sin → exp·cos → exp·sin … without
// bound and SEGFAULT the process. It now returns a closed form (verified by
// differentiation) and never crashes.
TEST_CASE("integrate: ∫exp(x)*sin(x) dx (cyclic by parts, was segfault)",
          "[7][integrate][byparts][regression]") {
    auto x = symbol("x");
    auto f = exp(x) * sin(x);
    auto F = integrate(f, x);
    // Must be actually evaluated, not the unevaluated Integral(...) marker.
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(simplify(diff(F, x) - f) == S::Zero());
}

TEST_CASE("integrate: ∫exp(x)*cos(x) dx (cyclic by parts)",
          "[7][integrate][byparts][regression]") {
    auto x = symbol("x");
    auto f = exp(x) * cos(x);
    auto F = integrate(f, x);
    REQUIRE(simplify(diff(F, x) - f) == S::Zero());
}

TEST_CASE("integrate: ∫exp(2x)*sin(3x) dx (distinct linear coefficients)",
          "[7][integrate][byparts][regression]") {
    auto x = symbol("x");
    auto f = exp(integer(2) * x) * sin(integer(3) * x);
    auto F = integrate(f, x);
    REQUIRE(simplify(diff(F, x) - f) == S::Zero());
}

TEST_CASE("integrate: ∫x^2*exp(x) dx still works (polynomial by parts intact)",
          "[7][integrate][byparts][regression]") {
    auto x = symbol("x");
    auto f = pow(x, integer(2)) * exp(x);
    auto F = integrate(f, x);
    REQUIRE(simplify(diff(F, x) - f) == S::Zero());
}

// ----- Linear u-substitution & rational integration --------------------------

TEST_CASE("integrate: ∫sin(2x+1) dx = -cos(2x+1)/2",
          "[7][integrate][linear_sub][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = sin(integer(2) * x + integer(1));
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(r->str(), "-cos(2*x + 1)/2"));
}

TEST_CASE("integrate: ∫1/(3x + 2) dx = log(3x + 2)/3",
          "[7][integrate][linear_sub][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(integer(3) * x + integer(2), integer(-1));
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(r->str(), "log(3*x + 2)/3"));
}

TEST_CASE("integrate: ∫(2x+5)^4 dx = (2x+5)^5/10",
          "[7][integrate][linear_sub][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(integer(2) * x + integer(5), integer(4));
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(r->str(), "(2*x + 5)**5/10"));
}

TEST_CASE("integrate: ∫exp(5x - 3) dx = exp(5x - 3)/5",
          "[7][integrate][linear_sub][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = exp(integer(5) * x - integer(3));
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(r->str(), "exp(5*x - 3)/5"));
}

TEST_CASE("integrate: ∫(3x+5)/((x-1)(x+2)) via apart",
          "[7][integrate][rational][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto num = integer(3) * x + integer(5);
    auto den = (x - integer(1)) * (x + integer(2));
    auto e = num / den;
    auto r = integrate(e, x);
    // Expected: (8/3)*log(x - 1) + (1/3)*log(x + 2). Compare via SymPy.
    auto resp = oracle.send({{"op", "integrate"},
                             {"expr", "(3*x + 5)/((x-1)*(x+2))"},
                             {"var", "x"}});
    REQUIRE(resp.ok);
    REQUIRE(oracle.equivalent(r->str(),
                              resp.raw.at("result").get<std::string>()));
}

TEST_CASE("integrate: improper rational gets polynomial part",
          "[7][integrate][rational][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // x^3/(x^2 - 1) — improper. Verify via diff-back: d/dx(integral) == integrand.
    auto e = pow(x, integer(3)) / (pow(x, integer(2)) - integer(1));
    auto r = integrate(e, x);
    auto check = oracle.equivalent(diff(r, x)->str(), e->str());
    REQUIRE(check);
}

// ----- Trig reductions -------------------------------------------------------

TEST_CASE("integrate: ∫sin²(x) dx via reduction (diff back)",
          "[7][integrate][trig_reduction][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(sin(x), integer(2));
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(diff(r, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫cos²(x) dx via reduction (diff back)",
          "[7][integrate][trig_reduction][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(cos(x), integer(2));
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(diff(r, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫sin(x)*cos(x) dx (diff back)",
          "[7][integrate][trig_reduction][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = sin(x) * cos(x);
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(diff(r, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫sin(2x)*cos(3x) dx product-to-sum",
          "[7][integrate][trig_reduction][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = sin(integer(2) * x) * cos(integer(3) * x);
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(diff(r, x)->str(), e->str()));
}

// ----- Trig powers ∫sinᵐcosⁿ (INT-18) ----------------------------------------
// Odd powers via u = sin/cos substitution; both-even via half-angle reduction.
TEST_CASE("integrate: odd-power trig (sin³, cos³, cos⁵, sin³cos²)",
          "[7][integrate][trig_power][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto db = [&](Expr e) {
        return oracle.equivalent(diff(integrate(e, x), x)->str(), e->str());
    };
    REQUIRE(db(pow(sin(x), integer(3))));
    REQUIRE(db(pow(cos(x), integer(3))));
    REQUIRE(db(pow(cos(x), integer(5))));
    REQUIRE(db(pow(sin(x), integer(3)) * pow(cos(x), integer(2))));
    // Affine argument.
    REQUIRE(db(pow(sin(integer(2) * x), integer(3))));
    // Closed-form check for sin³.
    REQUIRE(oracle.equivalent(integrate(pow(sin(x), integer(3)), x)->str(),
                              "cos(x)**3/3 - cos(x)"));
}

TEST_CASE("integrate: even-power trig (sin²cos², sin⁴) via half-angle",
          "[7][integrate][trig_power][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto db = [&](Expr e) {
        return oracle.equivalent(diff(integrate(e, x), x)->str(), e->str());
    };
    REQUIRE(db(pow(sin(x), integer(2)) * pow(cos(x), integer(2))));
    REQUIRE(db(pow(sin(x), integer(4))));
    REQUIRE(db(pow(cos(x), integer(4))));
}

// ----- tan powers and hyperbolic powers (INT-19) -----------------------------
TEST_CASE("integrate: tan(x)^n via reduction",
          "[7][integrate][tan_power][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto db = [&](Expr e) {
        return oracle.equivalent(diff(integrate(e, x), x)->str(), e->str());
    };
    REQUIRE(db(pow(tan(x), integer(3))));
    REQUIRE(db(pow(tan(x), integer(4))));
    REQUIRE(db(pow(tan(x), integer(5))));
    REQUIRE(db(pow(tan(integer(2) * x), integer(3))));  // affine argument
}

TEST_CASE("integrate: sinh^m * cosh^n powers",
          "[7][integrate][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto db = [&](Expr e) {
        return oracle.equivalent(diff(integrate(e, x), x)->str(), e->str());
    };
    REQUIRE(db(pow(sinh(x), integer(2))));
    REQUIRE(db(pow(cosh(x), integer(2))));
    REQUIRE(db(pow(sinh(x), integer(3))));
    REQUIRE(db(pow(cosh(x), integer(3))));
    REQUIRE(db(pow(sinh(x), integer(2)) * pow(cosh(x), integer(2))));
    REQUIRE(db(pow(sinh(integer(2) * x), integer(2))));  // affine argument
}

TEST_CASE("integrate: tanh(x)^n / coth(x)^n via reduction (INT-30)",
          "[7][integrate][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        // No unevaluated Integral marker leaks through.
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        return oracle.equivalent(diff(F, x)->str(), e->str());
    };
    REQUIRE(db(pow(tanh(x), integer(2))));
    REQUIRE(db(pow(tanh(x), integer(3))));
    REQUIRE(db(pow(tanh(x), integer(4))));
    REQUIRE(db(pow(coth(x), integer(2))));
    REQUIRE(db(pow(coth(x), integer(3))));
    REQUIRE(db(pow(coth(x), integer(4))));
    REQUIRE(db(pow(tanh(integer(2) * x), integer(3))));  // affine argument
}

TEST_CASE("integrate: ∫sin²(2x+1) dx — affine inside reduction",
          "[7][integrate][trig_reduction]") {
    // SymPy's plain simplify can't always reduce cos(2*(2x+1)) ↔ cos(4x+2),
    // so use a numeric-substitution sanity check at x = 1/2.
    auto x = symbol("x");
    auto e = pow(sin(integer(2) * x + integer(1)), integer(2));
    auto r = integrate(e, x);
    auto db = diff(r, x);
    // Verify the integration completed (didn't return Integral(...)).
    REQUIRE(r->type_id() != TypeId::Function);
    // Substitute a concrete value and compare numerics through evalf.
    auto e_at = evalf(subs(e, x, S::Half()), 20);
    auto db_at = evalf(subs(db, x, S::Half()), 20);
    auto& oracle = Oracle::instance();
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", e_at->str() + " - (" + db_at->str() + ")"},
                             {"prec", 30}, {"tol", 12}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

// ----- Integration by parts --------------------------------------------------

TEST_CASE("integrate: ∫log(x) dx = x*log(x) - x",
          "[7][integrate][parts][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = log(x);
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(r->str(), "x*log(x) - x"));
}

TEST_CASE("integrate: ∫log(2x+1) dx",
          "[7][integrate][parts][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = log(integer(2) * x + integer(1));
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(diff(r, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫x*exp(x) dx = (x-1)*exp(x)",
          "[7][integrate][parts][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = x * exp(x);
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(diff(r, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫x*sin(x) dx = -x*cos(x) + sin(x)",
          "[7][integrate][parts][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = x * sin(x);
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(r->str(), "-x*cos(x) + sin(x)"));
}

TEST_CASE("integrate: ∫x*cos(2x+1) dx (parts with affine)",
          "[7][integrate][parts][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = x * cos(integer(2) * x + integer(1));
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(diff(r, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫x²*exp(x) dx (recursive parts)",
          "[7][integrate][parts][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(x, integer(2)) * exp(x);
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(diff(r, x)->str(), e->str()));
}

// INT-13: ∫poly·cosh / poly·sinh fell through to the unevaluated Integral
// because the by-parts target set was {exp,sin,cos}. sinh/cosh now integrate
// by parts too (they don't cycle like exp·sin/cos).
TEST_CASE("integrate: ∫x*cosh(x) dx = x*sinh(x) - cosh(x)",
          "[7][integrate][parts][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = x * cosh(x);
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(r->str(), "x*sinh(x) - cosh(x)"));
}

TEST_CASE("integrate: ∫x*sinh(x) dx = x*cosh(x) - sinh(x)",
          "[7][integrate][parts][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = x * sinh(x);
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(r->str(), "x*cosh(x) - sinh(x)"));
}

TEST_CASE("integrate: ∫x²*cosh(x) dx (recursive parts)",
          "[7][integrate][parts][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(x, integer(2)) * cosh(x);
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(diff(r, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫x*cosh(2x+1) dx (parts with affine)",
          "[7][integrate][parts][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = x * cosh(integer(2) * x + integer(1));
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(diff(r, x)->str(), e->str()));
}

// INT-14: ∫log(x)^n and ∫poly·log(x)^n fell through to the unevaluated
// Integral — by parts only handled a single (power-1) log factor. log^n now
// integrates by repeated parts (u = log^n, dv = rest dx), recursing down to
// the ∫log case.
TEST_CASE("integrate: ∫log(x)² dx = x*log(x)² - 2x*log(x) + 2x",
          "[7][integrate][parts][log][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(log(x), integer(2));
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(r->str(), "x*log(x)**2 - 2*x*log(x) + 2*x"));
}

TEST_CASE("integrate: ∫log(x)³ dx (repeated parts)",
          "[7][integrate][parts][log][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(log(x), integer(3));
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(diff(r, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫x*log(x)² dx (poly times log power)",
          "[7][integrate][parts][log][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = x * pow(log(x), integer(2));
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(diff(r, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫log(2x)² dx (affine argument, b = 0)",
          "[7][integrate][parts][log][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(log(integer(2) * x), integer(2));
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(diff(r, x)->str(), e->str()));
}

// INT-15 (regression): ∫exp(x)/x once recursed exp(x)/x → exp(x)/x² → … forever
// in by-parts (u = x^(-1) is not a polynomial, so du grows). The by-parts guard
// made it terminate; it now closes to the exponential integral Ei(x) (EXPINT-1).
// The test completing at all still proves termination.
TEST_CASE("integrate: ∫exp(x)/x = Ei(x), terminating (INT-15)",
          "[7][integrate][expint][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = exp(x) / x;
    auto r = integrate(e, x);
    REQUIRE(r->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(r, x)->str(), e->str()));
}

// ----- Exponential/sine/cosine integral functions (regression, EXPINT-1) -----
// ∫sin(c·x)/x = Si(c·x), ∫cos(c·x)/x = Ci(c·x), ∫exp(c·x)/x = Ei(c·x) — the
// non-elementary special-integral functions, for a monomial argument c·x.
// Verified by differentiation against the oracle.
TEST_CASE("integrate: ∫sin(x)/x = Si(x)",
          "[7][integrate][expint][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = sin(x) / x;
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫cos(x)/x = Ci(x)",
          "[7][integrate][expint][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = cos(x) / x;
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫sinh(x)/x = Shi(x)",
          "[7][integrate][expint][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = sinh(x) / x;
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫cosh(x)/x = Chi(x)",
          "[7][integrate][expint][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = cosh(x) / x;
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫sin(3x)/x = Si(3x) (scaled argument)",
          "[7][integrate][expint][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = sin(integer(3) * x) / x;
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

// ----- ∫ of a special-integral function by parts (regression, EXPINT-BYPARTS)
// ∫f(x) dx = x·f(x) − ∫x·f'(x), closing because x·f' is elementary: ∫erf =
// x·erf+e^(-x²)/√π, ∫Si = x·Si+cos, ∫Ei = x·Ei−e^x, etc. Verified by oracle.
TEST_CASE("integrate: ∫erf(x) dx = x*erf(x) + exp(-x^2)/sqrt(pi)",
          "[7][integrate][expint][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = erf(x);
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫Si(x) dx = x*Si(x) + cos(x)",
          "[7][integrate][expint][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = sinint(x);
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫Ei(x) dx = x*Ei(x) - exp(x)",
          "[7][integrate][expint][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = expint_ei(x);
    auto F = integrate(e, x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
}

TEST_CASE("integrate: ∫erfi(x) and ∫Chi(x) by parts",
          "[7][integrate][expint][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    for (auto e : {erfi(x), coshint(x)}) {
        auto F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
    }
}

// ----- ∫ inverse trig/hyperbolic by parts (regression, INVTRIG-BYPARTS) ------
// ∫f(x) dx = x·f(x) − ∫x·f'(x); the x·f' term is a rational or x/√(quadratic),
// the latter handled by try_x_over_sqrt_quadratic. ∫asin = x·asin+√(1-x²),
// ∫atan = x·atan−log(1+x²)/2, ∫asinh = x·asinh−√(x²+1), etc. Oracle-verified.
TEST_CASE("integrate: ∫asin/acos/atan/acot by parts",
          "[7][integrate][invtrig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    for (auto e : {asin(x), acos(x), atan(x), acot(x)}) {
        auto F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
    }
}

TEST_CASE("integrate: ∫asinh/acosh/atanh by parts",
          "[7][integrate][invtrig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    for (auto e : {asinh(x), acosh(x), atanh(x)}) {
        auto F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
    }
}

TEST_CASE("integrate: ∫x/sqrt(a*x^2+c) = sqrt(a*x^2+c)/a",
          "[7][integrate][invtrig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // ∫x/√(1-x²) = -√(1-x²); ∫x/√(x²+1) = √(x²+1).
    auto e1 = x * pow(integer(1) - pow(x, integer(2)), rational(-1, 2));
    auto F1 = integrate(e1, x);
    REQUIRE(F1->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(F1, x)->str(), e1->str()));
    auto e2 = x * pow(pow(x, integer(2)) + integer(1), rational(-1, 2));
    auto F2 = integrate(e2, x);
    REQUIRE(oracle.equivalent(diff(F2, x)->str(), e2->str()));
}

TEST_CASE("integrate: completing the square for √-quadratics (INT-31)",
          "[7][integrate][invtrig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        return oracle.equivalent(diff(F, x)->str(), e->str());
    };
    Expr q1 = pow(x, integer(2)) + x + integer(1);        // x²+x+1
    Expr q2 = integer(2) * x - pow(x, integer(2));        // 2x-x² (a<0)
    Expr q3 = pow(x, integer(2)) + integer(2) * x + integer(5);  // x²+2x+5
    REQUIRE(db(pow(q1, rational(-1, 2))));                 // ∫1/√(x²+x+1)
    REQUIRE(db(pow(q2, rational(-1, 2))));                 // ∫1/√(2x-x²)
    REQUIRE(db(pow(q3, rational(1, 2))));                  // ∫√(x²+2x+5)
    // Linear numerator over a completed-square radical.
    REQUIRE(db((integer(2) * x + integer(3)) * pow(q1, rational(-1, 2))));
    REQUIRE(db((x - integer(1)) * pow(q3, rational(-1, 2))));
}

TEST_CASE("integrate: improper rational over an irreducible quadratic (INT-32)",
          "[7][integrate][rational][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        return oracle.equivalent(diff(F, x)->str(), e->str());
    };
    Expr d = pow(x, integer(2)) + integer(1);   // irreducible x²+1
    // Polynomial division leaves a proper part apart() can't split further.
    REQUIRE(db(pow(x, integer(2)) / d));         // ∫x²/(x²+1) = x − atan(x)
    REQUIRE(db(pow(x, integer(3)) / d));         // ∫x³/(x²+1)
    REQUIRE(db(pow(x, integer(4)) / d));
    REQUIRE(db((pow(x, integer(2)) + integer(2)) / d));
}

TEST_CASE("integrate: polynomial × inverse-trig by parts (INT-32)",
          "[7][integrate][invtrig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        return oracle.equivalent(diff(F, x)->str(), e->str());
    };
    REQUIRE(db(x * atan(x)));                     // x·atan
    REQUIRE(db(pow(x, integer(2)) * atan(x)));    // x²·atan
    REQUIRE(db(x * atanh(x)));                     // x·atanh
    REQUIRE(db(x * acot(x)));                      // x·acot
    REQUIRE(db((x + integer(1)) * atan(x)));      // (x+1)·atan
}

TEST_CASE("integrate: trig × hyperbolic and exp × hyperbolic products (INT-34)",
          "[7][integrate][hyperbolic][regression]") {
    auto x = symbol("x");
    // The antiderivatives come out in exponential form while the integrand is in
    // sinh/cosh form; SymPy's simplify can't bridge the two, and its numeric
    // `.equals` sampling is flaky here. Verify deterministically: diff(F) − e
    // must evaluate to ~0 at fixed rational points.
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        Expr resid = diff(F, x) - e;
        for (Expr pt : {rational(1, 2), rational(7, 5), rational(-3, 4)}) {
            double d = std::stod(evalf(subs(resid, x, pt))->str());
            if (std::abs(d) > 1e-9) return false;
        }
        return true;
    };
    REQUIRE(db(sin(x) * sinh(x)));
    REQUIRE(db(cos(x) * cosh(x)));
    REQUIRE(db(sin(x) * cosh(x)));
    REQUIRE(db(cos(x) * sinh(x)));
    REQUIRE(db(exp(x) * sinh(x)));                // collapses to pure exponential
    REQUIRE(db(exp(integer(2) * x) * cosh(x)));
    REQUIRE(db(sin(integer(2) * x) * sinh(integer(3) * x)));  // affine arguments
}

TEST_CASE("integrate: polynomial × trig/hyperbolic power by parts (INT-35)",
          "[7][integrate][parts][regression]") {
    auto x = symbol("x");
    // The antiderivatives keep sinⁿ/cosⁿ or multiple-angle terms that SymPy's
    // simplify won't always reduce against the integrand, so verify the
    // derivative numerically at fixed rational points (deterministic).
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        Expr resid = diff(F, x) - e;
        for (Expr pt : {rational(1, 2), rational(7, 5), rational(-3, 4)}) {
            double d = std::stod(evalf(subs(resid, x, pt))->str());
            if (std::abs(d) > 1e-9) return false;
        }
        return true;
    };
    REQUIRE(db(x * pow(cos(x), integer(2))));               // ∫x·cos²x
    REQUIRE(db(x * pow(sin(x), integer(2))));               // ∫x·sin²x
    REQUIRE(db(x * pow(sin(x), integer(3))));               // odd power
    REQUIRE(db(pow(x, integer(2)) * pow(cos(x), integer(2))));  // x²·cos²x
    REQUIRE(db(x * pow(cosh(x), integer(2))));              // hyperbolic power
    REQUIRE(db(x * pow(cos(integer(2) * x), integer(2))));  // affine argument
}

TEST_CASE("integrate: 1/(quadratic)^n reduction formula (INT-37)",
          "[7][integrate][rational][regression]") {
    auto x = symbol("x");
    // The reduction emits atan/log mixed with rational terms that SymPy's
    // simplify won't reduce against the integrand, so verify the derivative
    // numerically at fixed rational points (deterministic).
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        Expr resid = diff(F, x) - e;
        for (Expr pt : {rational(1, 2), rational(7, 5), rational(-3, 4)}) {
            double d = std::stod(evalf(subs(resid, x, pt))->str());
            if (std::abs(d) > 1e-9) return false;
        }
        return true;
    };
    Expr q = pow(x, integer(2)) + integer(1);
    REQUIRE(db(pow(q, integer(-2))));                              // 1/(x²+1)²
    REQUIRE(db(pow(q, integer(-3))));                              // 1/(x²+1)³
    REQUIRE(db(pow(pow(x, integer(2)) + integer(4), integer(-2))));  // 1/(x²+4)²
    REQUIRE(db(pow(integer(2) * pow(x, integer(2)) + integer(3),
                   integer(-2))));                                 // 1/(2x²+3)²
    // Completing the square: 1/(x²+2x+5)².
    REQUIRE(db(pow(pow(x, integer(2)) + integer(2) * x + integer(5),
                   integer(-2))));
    // Linear numerator over a power of an irreducible quadratic (INT-39).
    Expr q2 = pow(x, integer(2)) + integer(1);
    REQUIRE(db((x + integer(1)) * pow(q2, integer(-2))));          // (x+1)/(x²+1)²
    REQUIRE(db((integer(2) * x + integer(3)) * pow(q2, integer(-2))));
    REQUIRE(db((integer(3) * x + integer(2)) * pow(q2, integer(-3))));
    // Linear numerator over a completed square.
    REQUIRE(db((x - integer(1))
               * pow(pow(x, integer(2)) + integer(2) * x + integer(5),
                     integer(-2))));
}

TEST_CASE("integrate: rational with repeated factors via partial fractions (INT-38)",
          "[7][integrate][rational][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        return oracle.equivalent(diff(F, x)->str(), e->str());
    };
    // apart now decomposes repeated factors; try_rational accepts (poly)^(-n)
    // denominators. ∫1/((x-1)²(x+1)), ∫1/(x²(x+1)), ∫x³/(x²+1)².
    REQUIRE(db(pow(pow(x - integer(1), integer(2)) * (x + integer(1)),
                   integer(-1))));
    REQUIRE(db(pow(pow(x, integer(2)) * (x + integer(1)), integer(-1))));
    REQUIRE(db(pow(x, integer(3))
               * pow(pow(x, integer(2)) + integer(1), integer(-2))));
    // Regression: distinct-factor and irreducible-quadratic cases still close.
    REQUIRE(db(pow((x - integer(1)) * (x + integer(1)), integer(-1))));
    REQUIRE(db(pow(pow(x, integer(3)) + integer(1), integer(-1))));
}

TEST_CASE("integrate: polynomial over √(quadratic) via reduction (INT-40)",
          "[7][integrate][invtrig][regression]") {
    auto x = symbol("x");
    // ∫xⁿ/√(quadratic) reduction; results mix √ and asin/asinh that simplify
    // can't reduce against the integrand, so check the derivative numerically.
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        Expr resid = diff(F, x) - e;
        for (Expr pt : {rational(1, 3), rational(-2, 5), rational(3, 7)}) {
            double d = std::stod(evalf(subs(resid, x, pt))->str());
            if (std::abs(d) > 1e-9) return false;
        }
        return true;
    };
    REQUIRE(db(pow(x, integer(2))
               * pow(integer(1) - pow(x, integer(2)), rational(-1, 2))));   // x²/√(1−x²)
    REQUIRE(db(pow(x, integer(3))
               * pow(pow(x, integer(2)) + integer(1), rational(-1, 2))));   // x³/√(x²+1)
    REQUIRE(db(pow(x, integer(2))
               * pow(pow(x, integer(2)) + integer(2) * x + integer(5),
                     rational(-1, 2))));                                     // completed square
    // The INT-32-deferred polynomial × inverse-function cases now close.
    REQUIRE(db(pow(x, integer(2)) * asin(x)));
    REQUIRE(db(pow(x, integer(2)) * asinh(x)));
}

TEST_CASE("integrate: log of a polynomial by parts (INT-41)",
          "[7][integrate][parts][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        return oracle.equivalent(diff(F, x)->str(), e->str());
    };
    // ∫log(g) = x·log(g) − ∫x·g'/g, closing when the remainder is rational.
    REQUIRE(db(log(pow(x, integer(2)) + integer(1))));        // x·log + atan
    REQUIRE(db(log(pow(x, integer(2)) - integer(1))));        // reducible → logs
    REQUIRE(db(log(pow(x, integer(2)) + x + integer(1))));    // irreducible → log+atan
    // Regression: log(affine) closed form still works.
    REQUIRE(db(log(x)));
    REQUIRE(db(log(integer(2) * x + integer(3))));
}

TEST_CASE("integrate: heurisch substitution into an irreducible quadratic (INT-36)",
          "[7][integrate][heurisch][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        return oracle.equivalent(diff(F, x)->str(), e->str());
    };
    // After the u = g substitution the integrand becomes c/(1+u²); the table and
    // try_rational don't close a bare irreducible quadratic, so the dedicated
    // arctan helper now backs them up. ∫g'/(1+g²) = atan(g).
    REQUIRE(db(cos(x) / (integer(1) + pow(sin(x), integer(2)))));   // atan(sin x)
    REQUIRE(db(sin(x) / (integer(1) + pow(cos(x), integer(2)))));   // -atan(cos x)
    REQUIRE(db(exp(x) / (integer(1) + pow(exp(x), integer(2)))));   // atan(exp x)
    REQUIRE(db(pow(x * (integer(1) + pow(log(x), integer(2))), integer(-1))));
}

TEST_CASE("integrate: Weierstrass substitution for rational trig (INT-33)",
          "[7][integrate][weierstrass][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        return oracle.equivalent(diff(F, x)->str(), e->str());
    };
    // ∫1/(a + b·cos x), ∫1/(a + b·sin x), and a mixed denominator. The
    // half-angle output stresses SymPy's simplify, so these are the subset the
    // oracle reliably confirms (others integrate correctly — verified by direct
    // numeric sampling — but SymPy can't reduce their tan(x/2) derivative form).
    REQUIRE(db(pow(integer(2) + cos(x), integer(-1))));               // a²>b²: atan
    REQUIRE(db(pow(integer(3) + integer(5) * cos(x), integer(-1))));  // a²<b²: log
    REQUIRE(db(pow(integer(2) - cos(x), integer(-1))));               // negated b
    REQUIRE(db(pow(integer(1) + sin(x), integer(-1))));               // a=b: rational
    REQUIRE(db(pow(integer(4) + integer(5) * sin(x), integer(-1))));  // a²<b²
    REQUIRE(db(pow(integer(2) + cos(x) + sin(x), integer(-1))));      // sin+cos
    // The dedicated integrators still win for simple trig (not Weierstrass).
    REQUIRE(integrate(sin(x), x)->str() == "-cos(x)");
    // Regression: a trig integrand that is NOT rational in tan(x/2) — √(tan x)
    // substitutes to √(2t/(1−t²)), a non-elementary algebraic integral. The
    // is_rational_in guard must bail to the marker instead of looping forever.
    REQUIRE(integrate(pow(tan(x), rational(1, 2)), x)->str().find("Integral(")
            != std::string::npos);
    REQUIRE(integrate(pow(sin(x), rational(1, 2)), x)->str().find("Integral(")
            != std::string::npos);
    // Regression: a trig function raised to a power (tan²/sec²) substitutes to a
    // high-degree nested rational in t whose normalisation/integration runs
    // away. The has_trig_power_of guard must bail to the marker, not hang.
    REQUIRE(integrate(pow(integer(1) + pow(tan(x), integer(2)), integer(-1)), x)
                ->str().find("Integral(") != std::string::npos);
    // sec²·(1 + tan²)⁻¹ reduces to 1 (sec² = 1 + tan²), so trigsimp resolves it
    // before the Weierstrass path is reached: it now closes (∫1 = x, returned
    // here as the equivalent atan(tan x)) instead of bailing — must not hang.
    {
        Expr F = integrate(pow(sec(x), integer(2))
                               * pow(integer(1) + pow(tan(x), integer(2)),
                                     integer(-1)),
                           x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        REQUIRE(oracle.equivalent(simplify(diff(F, x))->str(), "1"));
    }
    // But trig to the FIRST power inside a polynomial denominator still works.
    REQUIRE(integrate(pow(integer(1) + tan(x), integer(-1)), x)->str()
                .find("Integral(") == std::string::npos);
}

// ----- manualintegrate orchestrator ------------------------------------------

TEST_CASE("manualintegrate: returns Some on tractable integrand",
          "[7][manualintegrate][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = x * sin(x);
    auto r = manualintegrate(e, x);
    REQUIRE(r.has_value());
    REQUIRE(oracle.equivalent((*r)->str(), "-x*cos(x) + sin(x)"));
}

TEST_CASE("manualintegrate: returns None on intractable integrand",
          "[7][manualintegrate]") {
    auto x = symbol("x");
    // exp(x³) has no elementary antiderivative and no special-function table
    // entry (the Gaussian exp(c·x²) → erf/erfi rule needs a degree-2 exponent).
    auto e = exp(pow(x, integer(3)));
    auto r = manualintegrate(e, x);
    REQUIRE(!r.has_value());
}

TEST_CASE("manualintegrate: chain reaches log(x) via parts",
          "[7][manualintegrate][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = log(x);
    auto r = manualintegrate(e, x);
    REQUIRE(r.has_value());
    REQUIRE(oracle.equivalent((*r)->str(), "x*log(x) - x"));
}

// ----- Heurisch (chain-rule-reverse u-substitution) --------------------------

TEST_CASE("integrate: ∫2x*exp(x²) dx = exp(x²) (heurisch)",
          "[7][integrate][heurisch][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = integer(2) * x * exp(pow(x, integer(2)));
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(r->str(), "exp(x**2)"));
}

TEST_CASE("integrate: ∫x*exp(x²) dx = exp(x²)/2 (heurisch with constant)",
          "[7][integrate][heurisch][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = x * exp(pow(x, integer(2)));
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(r->str(), "exp(x**2)/2"));
}

TEST_CASE("integrate: ∫cos(x²)*2x dx = sin(x²)",
          "[7][integrate][heurisch][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = cos(pow(x, integer(2))) * integer(2) * x;
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(r->str(), "sin(x**2)"));
}

TEST_CASE("integrate: ∫sin(x³)*3x² dx = -cos(x³)",
          "[7][integrate][heurisch][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = sin(pow(x, integer(3))) * integer(3) * pow(x, integer(2));
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(r->str(), "-cos(x**3)"));
}

TEST_CASE("integrate: ∫(2x)/(x²+1) dx = log(x²+1) (logarithmic derivative)",
          "[7][integrate][heurisch][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = (integer(2) * x) / (pow(x, integer(2)) + integer(1));
    auto r = integrate(e, x);
    REQUIRE(oracle.equivalent(diff(r, x)->str(), e->str()));
}

#include <sympp/integrals/quadrature.hpp>

// ----- Numeric quadrature (vpaintegral) --------------------------------------

TEST_CASE("vpaintegral: ∫_0^1 x dx = 1/2",
          "[7][vpaintegral][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto v = vpaintegral(x, x, S::Zero(), S::One(), 15);
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", v->str() + " - 0.5"},
                             {"prec", 30}, {"tol", 10}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

TEST_CASE("vpaintegral: ∫_0^1 x² dx = 1/3",
          "[7][vpaintegral][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto v = vpaintegral(pow(x, integer(2)), x, S::Zero(), S::One(), 15);
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", v->str() + " - 1/3"},
                             {"prec", 30}, {"tol", 10}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

TEST_CASE("vpaintegral: ∫_0^pi sin(x) dx = 2",
          "[7][vpaintegral][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto v = vpaintegral(sin(x), x, S::Zero(), S::Pi(), 15);
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", v->str() + " - 2"},
                             {"prec", 30}, {"tol", 10}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

TEST_CASE("vpaintegral: ∫_0^1 exp(x) dx = e - 1",
          "[7][vpaintegral][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto v = vpaintegral(exp(x), x, S::Zero(), S::One(), 15);
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", v->str() + " - (E - 1)"},
                             {"prec", 30}, {"tol", 10}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

TEST_CASE("vpaintegral: ∫_0^1 exp(x²) dx (intractable symbolically)",
          "[7][vpaintegral][oracle]") {
    auto& oracle = Oracle::instance();
    // Known value: ≈ 1.4626517459071816 — not expressible in elementary
    // closed form. SymPP's symbolic integrator returns an Integral marker;
    // vpaintegral computes it numerically.
    auto x = symbol("x");
    auto v = vpaintegral(exp(pow(x, integer(2))), x, S::Zero(), S::One(), 15);
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", v->str() + " - 1.4626517459071816"},
                             {"prec", 30}, {"tol", 10}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

// INT-DEF-1: definite integration with an infinite bound now evaluates the
// antiderivative as a LIMIT, not by literal substitution of oo (which gave nan
// from the ∞·0 boundary term). Recovers the Gamma integrals ∫₀^∞ x^n·e^(-x) = n!.
TEST_CASE("integrate: improper integrals over [0, oo) (INT-DEF-1)",
          "[7][integrate][definite][regression]") {
    auto x = symbol("x");
    auto oo = S::Infinity();
    auto negx = mul(S::NegativeOne(), x);
    // ∫₀^∞ x^n·e^(-x) dx = n!.
    REQUIRE(integrate(x * exp(negx), x, S::Zero(), oo) == integer(1));
    REQUIRE(integrate(pow(x, integer(2)) * exp(negx), x, S::Zero(), oo)
            == integer(2));
    REQUIRE(integrate(pow(x, integer(3)) * exp(negx), x, S::Zero(), oo)
            == integer(6));
    REQUIRE(integrate(pow(x, integer(4)) * exp(negx), x, S::Zero(), oo)
            == integer(24));
    // Scaled exponent: ∫₀^∞ x·e^(-2x) dx = 1/4.
    REQUIRE(integrate(x * exp(mul(integer(-2), x)), x, S::Zero(), oo)
            == rational(1, 4));
    // Finite-bound integrals are unchanged.
    REQUIRE(integrate(pow(x, integer(2)), x, S::Zero(), integer(1))
            == rational(1, 3));
}

// INT-IMPROPER-1: improper rational functions (deg numerator ≥ deg denominator)
// over a LINEAR denominator used to come back unevaluated. try_rational does the
// polynomial division, but when apart left the proper remainder as a single
// c/(x+a) term the code only closed a degree-2 denominator and dropped the rest
// to the Integral marker. It now integrates the remainder through the general
// integrator (affine-log rule), so ∫x/(x+1) = x − log(x+1), etc.
TEST_CASE("integrate: improper rational over a linear denominator (INT-IMPROPER-1)",
          "[7][integrate][rational][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    for (const Expr& e : {x / (x + integer(1)),
                          pow(x, integer(2)) / (x + integer(1)),
                          (pow(x, integer(2)) + integer(1)) / (x - integer(1)),
                          (x + integer(1)) / x}) {
        auto F = integrate(e, x);
        INFO("integrand: " << e->str() << "  antiderivative: " << F->str());
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
    }
    // Quadratic-denominator improper cases still close (no regression).
    auto g = pow(x, integer(3)) / (pow(x, integer(2)) - integer(1));
    auto Fg = integrate(g, x);
    REQUIRE(Fg->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(diff(Fg, x)->str(), g->str()));
}

// INT-RECIP-1: 1/cos(x) and 1/sin(x) written as Pow(cos/sin, -1) used to fall
// through to the marker, even though the Sec/Csc functions and cos(x)**(-2)
// already integrated. The reciprocal first power now routes to the same sec/csc
// antiderivatives (including affine arguments).
TEST_CASE("integrate: reciprocal trig as a Pow (INT-RECIP-1)",
          "[7][integrate][trig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    for (const Expr& e : {pow(cos(x), integer(-1)),
                          pow(sin(x), integer(-1)),
                          pow(cos(integer(2) * x + integer(1)), integer(-1)),
                          pow(sin(integer(3) * x), integer(-1))}) {
        auto F = integrate(e, x);
        INFO("integrand: " << e->str() << "  antiderivative: " << F->str());
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
    }
    // 1/cos(x) and sec(x) now give the same antiderivative.
    REQUIRE(integrate(pow(cos(x), integer(-1)), x)
            == integrate(sec(x), x));
}

// INT-RECIP-2: hyperbolic analogue of INT-RECIP-1 — 1/cosh(x) and 1/sinh(x)
// written as Pow(cosh/sinh, -1) route to the sech/csch antiderivatives.
TEST_CASE("integrate: reciprocal hyperbolic as a Pow (INT-RECIP-2)",
          "[7][integrate][hyperbolic][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    for (const Expr& e : {pow(cosh(x), integer(-1)),
                          pow(sinh(x), integer(-1)),
                          pow(cosh(integer(2) * x), integer(-1)),
                          pow(sinh(integer(3) * x + integer(1)), integer(-1))}) {
        auto F = integrate(e, x);
        INFO("integrand: " << e->str() << "  antiderivative: " << F->str());
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
    }
    // 1/cosh(x) and sech(x) now give the same antiderivative.
    REQUIRE(integrate(pow(cosh(x), integer(-1)), x)
            == integrate(sech(x), x));
}
