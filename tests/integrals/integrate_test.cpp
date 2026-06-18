// integrate tests.

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <string>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/combinatorial.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/functions/miscellaneous.hpp>
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

// INT-ABS-DEF-1: a definite integral of an abs/sign integrand has no elementary
// antiderivative, but is piecewise-smooth — split at the sign changes (roots of
// the argument) and sum. Previously the unevaluated marker was substituted at the
// bounds, producing garbage like −Integral(1, −1) + Integral(1, 1).
TEST_CASE("integrate: definite integral of abs/sign by sign-splitting (INT-ABS-DEF-1)",
          "[7][integrate][regression]") {
    auto x = symbol("x");
    auto I = [&](const Expr& e, int a, int b) {
        return integrate(e, x, integer(a), integer(b));
    };
    REQUIRE(I(abs(x), -1, 1) == integer(1));               // two triangles
    REQUIRE(I(abs(x), -2, 3) == rational(13, 2));          // asymmetric
    REQUIRE(I(abs(x - integer(1)), 0, 2) == integer(1));   // shifted root
    REQUIRE(I(x * abs(x), -1, 1) == integer(0));           // odd integrand
    REQUIRE(I(sign(x), -1, 2) == integer(1));              // ∫sign = −1+2
    REQUIRE(I(abs(x) + pow(x, integer(2)), -1, 1)
            == rational(5, 3));  // marker buried in a sum: 1 + 2/3
    // No interior root: the integrand is smooth over [0,1], plain ∫x = 1/2.
    REQUIRE(I(abs(x), 0, 1) == rational(1, 2));
    // |cos x| over [0, π] = 2 (root at π/2).
    auto two = integrate(abs(cos(x)), x, S::Zero(), S::Pi());
    REQUIRE(two == integer(2));
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

// INT-LOGSUB-1: integrands built from log(x) close under u = log(x)
// (x = eᵘ, dx = eᵘ du): ∫cos(log x), ∫sin(log x), ∫cos(2·log x) become cyclic
// exp·trig integrals in u, and ∫log(log x)/x becomes ∫log(u). Verified by
// diff-back against the oracle.
TEST_CASE("integrate: substitution u = log(x) for log-composite integrands (INT-LOGSUB-1)",
          "[7][integrate][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        return oracle.equivalent(diff(F, x)->str(), e->str());
    };
    REQUIRE(db(cos(log(x))));                          // ∫cos(log x)
    REQUIRE(db(sin(log(x))));                          // ∫sin(log x)
    REQUIRE(db(cos(integer(2) * log(x))));             // affine inner 2·log x
    REQUIRE(db(log(log(x)) / x));                      // ∫log(log x)/x → ∫log u
    REQUIRE(db(sin(log(x)) * cos(log(x))));            // product of logs
    // Known closed form for ∫cos(log x).
    REQUIRE(oracle.equivalent(integrate(cos(log(x)), x)->str(),
                              "x*sin(log(x))/2 + x*cos(log(x))/2"));
    // Integrands without log(x) are untouched: ∫1/x = log x, ∫x·log x by parts.
    REQUIRE(oracle.equivalent(integrate(pow(x, integer(-1)), x)->str(),
                              "log(x)"));
    REQUIRE(oracle.equivalent(integrate(x * log(x), x)->str(),
                              "x**2*log(x)/2 - x**2/4"));
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

// INT-TANSEC-1: ∫tan^m·sec^n (and cot^m·csc^n). m odd → u = sec; m even → rewrite
// tan²=sec²−1 to pure sec powers. ∫tan³·sec = sec³/3 − sec, etc. Verified by
// diff-back at sample points (the closed forms keep sec/tan that SymPy's
// simplify won't reduce against the integrand symbolically).
TEST_CASE("integrate: tan^m·sec^n products (INT-TANSEC-1)",
          "[7][integrate][trig][regression]") {
    auto x = symbol("x");
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        Expr resid = diff(F, x) - e;
        for (Expr pt : {rational(2, 5), rational(3, 7), rational(-1, 4)}) {
            double d = std::stod(evalf(subs(resid, x, pt))->str());
            if (std::fabs(d) > 1e-7) return false;
        }
        return true;
    };
    REQUIRE(db(pow(tan(x), integer(3)) * sec(x)));         // m odd
    REQUIRE(db(pow(tan(x), integer(2)) * sec(x)));         // m even → sec powers
    REQUIRE(db(tan(x) * pow(sec(x), integer(3))));         // m=1
    REQUIRE(db(pow(tan(x), integer(3)) * pow(sec(x), integer(3))));
    REQUIRE(db(pow(tan(x), integer(2)) * pow(sec(x), integer(2))));
    REQUIRE(db(pow(cot(x), integer(3)) * csc(x)));         // cot/csc analogue
    REQUIRE(db(pow(cot(x), integer(2)) * csc(x)));
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

// INT-ARCTAN-PARAM-1: ∫1/(quadratic) with SYMBOLIC positive coefficients —
// ∫1/(x²+a²) = atan(x/a)/a (a > 0), ∫1/(ax²+b) = atan(x√(a/b))/√(ab). Fires only
// when the discriminant is provably positive, matching SymPy under positivity
// assumptions. Relies on the Mul-positivity fix (is_positive(4·a²) = true).
TEST_CASE("integrate: arctan over a symbolic-coefficient quadratic (INT-ARCTAN-PARAM-1)",
          "[7][integrate][arctan][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto a = symbol("a", AssumptionMask{}.set_positive(true));
    auto b = symbol("b", AssumptionMask{}.set_positive(true));
    // Symbolic antiderivatives carry an (a > 0)-branch (√(4a²) = 2a), which the
    // oracle can't reduce without assumptions — so verify after substituting
    // concrete positive values for the parameters.
    auto chk = [&](const Expr& integrand, const ExprMap<Expr>& sub) {
        auto F = integrate(integrand, x);
        INFO("integrand: " << integrand->str() << "  F: " << F->str());
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        Expr Fc = xreplace(F, sub);
        Expr ic = xreplace(integrand, sub);
        REQUIRE(oracle.equivalent(diff(Fc, x)->str(), ic->str()));
    };
    // 1/(x²+a²), 1/(x²+a), 1/(a·x²+b).
    chk(pow(pow(x, integer(2)) + pow(a, integer(2)), integer(-1)),
        {{a, integer(2)}});
    chk(pow(pow(x, integer(2)) + a, integer(-1)), {{a, integer(3)}});
    chk(pow(a * pow(x, integer(2)) + b, integer(-1)),
        {{a, integer(2)}, {b, integer(5)}});
    // Negative discriminant (Δ = 4a² > 0) → the log form: 1/(a²−x²), 1/(x²−a²).
    chk(pow(pow(a, integer(2)) - pow(x, integer(2)), integer(-1)),
        {{a, integer(3)}});
    chk(pow(pow(x, integer(2)) - pow(a, integer(2)), integer(-1)),
        {{a, integer(3)}});
    // is_positive propagates through an even power and a positive coefficient.
    REQUIRE(is_positive(integer(4) * pow(a, integer(2))) == true);
}

// INT-SQRTQUAD-PARAM-1: ∫1/√(quadratic) with symbolic positive coefficients —
// ∫1/√(x²+a²) = asinh(x/a), ∫1/√(a²−x²) = asin(x/a). The sign-gated branches in
// try_sqrt_quadratic already handled symbolic coefficients; only a rational-only
// gate blocked them (now relaxed). These match SymPy's forms exactly.
TEST_CASE("integrate: sqrt-quadratic with symbolic coefficients (INT-SQRTQUAD-PARAM-1)",
          "[7][integrate][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto a = symbol("a", AssumptionMask{}.set_positive(true));
    // ∫1/√(x²+a²) = asinh(x/a), exactly SymPy's form.
    REQUIRE(oracle.equivalent(
        integrate(pow(pow(x, integer(2)) + pow(a, integer(2)), rational(-1, 2)),
                  x)
            ->str(),
        "asinh(x/a)"));
    // ∫1/√(a²−x²) = asin(x/a).
    REQUIRE(oracle.equivalent(
        integrate(pow(pow(a, integer(2)) - pow(x, integer(2)), rational(-1, 2)),
                  x)
            ->str(),
        "asin(x/a)"));
    // ∫1/√(x²+a) = asinh(x/√a).
    REQUIRE(oracle.equivalent(
        integrate(pow(pow(x, integer(2)) + a, rational(-1, 2)), x)->str(),
        "asinh(x/sqrt(a))"));
    // ∫√(a²−x²) = (a²·asin(x/a) + x·√(a²−x²))/2 (verified by diff-back at a=2).
    auto F = integrate(
        pow(pow(a, integer(2)) - pow(x, integer(2)), rational(1, 2)), x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    Expr Fc = xreplace(F, ExprMap<Expr>{{a, integer(2)}});
    Expr ic = pow(integer(4) - pow(x, integer(2)), rational(1, 2));
    REQUIRE(oracle.equivalent(diff(Fc, x)->str(), ic->str()));
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

// INT-XSQRTQUAD-NUM-1: a non-polynomial numerator over √(quadratic) must NOT be
// treated as a constant coefficient. Poly(asin(x), x) sees asin(x) as a degree-0
// "coefficient", so ∫asin(x)/√(1−x²) was wrongly pulled out as asin(x)·∫1/√Q =
// asin(x)² (a factor-of-2 error; correct is asin²/2). Guard the numerator
// coefficients var-free; the residual is then verified by heurisch's diff-back.
TEST_CASE("integrate: non-polynomial numerator over sqrt-quadratic (INT-XSQRTQUAD-NUM-1)",
          "[7][integrate][invtrig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto isqrt = pow(integer(1) - pow(x, integer(2)), rational(-1, 2));  // 1/√(1−x²)
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        return oracle.equivalent(diff(F, x)->str(), e->str());
    };
    // ∫asin/√(1−x²) = asin²/2 (was the wrong asin²); ∫asin²/√ = asin³/3.
    REQUIRE(db(asin(x) * isqrt));
    REQUIRE(db(pow(asin(x), integer(2)) * isqrt));
    REQUIRE(db(acos(x) * isqrt));  // −acos²/2 (was the wrong acos·asin)
    // Genuine numerator forms (bare x, linear) are unaffected.
    REQUIRE(db(x * isqrt));                                   // −√(1−x²)
    REQUIRE(db((integer(2) * x + integer(1)) * isqrt));       // −2√(1−x²)+asin
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

// INT-POLYSQRTQUAD-1: ∫P(x)·√(quadratic) for an even-power P. The u = Q
// substitution closes odd powers (∫x·√(1−x²)), but even powers fell through;
// rewriting P·√Q = (P·Q)/√Q routes them to the polynomial-over-√(quadratic)
// reduction. Verified by diff-back.
TEST_CASE("integrate: polynomial × sqrt(quadratic) (INT-POLYSQRTQUAD-1)",
          "[7][integrate][invtrig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        return oracle.equivalent(diff(F, x)->str(), e->str());
    };
    auto sq = [&](const Expr& q) { return pow(q, rational(1, 2)); };
    auto omx2 = integer(1) - pow(x, integer(2));
    REQUIRE(db(pow(x, integer(2)) * sq(omx2)));           // ∫x²·√(1−x²)
    REQUIRE(db(pow(x, integer(4)) * sq(omx2)));           // ∫x⁴·√(1−x²)
    REQUIRE(db(pow(x, integer(2))
               * sq(integer(4) - pow(x, integer(2)))));    // non-unit radicand
    // Higher half-power: ∫x²·(1−x²)^(3/2).
    REQUIRE(db(pow(x, integer(2)) * pow(omx2, rational(3, 2))));
    // Odd-power and bare-radical cases keep their existing (cleaner) forms.
    REQUIRE(oracle.equivalent(integrate(x * sq(omx2), x)->str(),
                              "-(1-x**2)**(3/2)/3"));
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
// INT-CONSTBASEEXP-1: constant-base exponential ∫P(x)·a^(b·x+c). SymPP integrated
// the natural base eˣ but left aˣ unevaluated (exp(x·ln a) even canonicalizes back
// to a^x). Reduces by parts to ∫a^g = a^g/(b·ln a). Verified by diff-back.
TEST_CASE("integrate: constant-base exponential ∫P(x)·a^x (INT-CONSTBASEEXP-1)",
          "[7][integrate][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        return oracle.equivalent(diff(F, x)->str(), e->str());
    };
    REQUIRE(db(pow(integer(2), x)));                       // ∫2ˣ
    REQUIRE(db(x * pow(integer(2), x)));                   // ∫x·2ˣ
    REQUIRE(db(pow(x, integer(2)) * pow(integer(2), x)));  // ∫x²·2ˣ
    REQUIRE(db(x * pow(integer(3), x)));                   // different base
    REQUIRE(db((x + integer(1)) * pow(integer(2), x)));    // mixed polynomial
    REQUIRE(db(x * pow(integer(2), mul(S::NegativeOne(), x))));  // decaying 2^(−x)
    REQUIRE(db(pow(integer(2), integer(3) * x)));          // affine exponent 3x
    // ∫2ˣ closed form.
    REQUIRE(oracle.equivalent(integrate(pow(integer(2), x), x)->str(),
                              "2**x/log(2)"));
    // Natural base is still handled by the elementary path, unchanged.
    REQUIRE(oracle.equivalent(integrate(x * exp(x), x)->str(),
                              "x*exp(x) - exp(x)"));
}

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

// INT-GAUSSMOMENT-1: polynomial × Gaussian ∫P(x)·exp(c·x²) dx — the Gaussian
// moments. Reduces ∫xⁿ·exp(c·x²) by parts to ∫exp(c·x²) (erf/erfi, even n) and
// ∫x·exp(c·x²) = exp(c·x²)/(2c) (odd n). Without it ∫x²·exp(−x²) was unevaluated,
// and the improper ∫₀^∞ x²·exp(−x²) produced garbage. Verified by diff-back.
TEST_CASE("integrate: polynomial times Gaussian → erf moments (INT-GAUSSMOMENT-1)",
          "[7][integrate][erf][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        return oracle.equivalent(diff(F, x)->str(), e->str());
    };
    auto g = exp(mul(S::NegativeOne(), pow(x, integer(2))));  // exp(−x²)
    REQUIRE(db(pow(x, integer(2)) * g));               // even: erf moment
    REQUIRE(db(pow(x, integer(3)) * g));               // odd: elementary
    REQUIRE(db(pow(x, integer(4)) * g));
    REQUIRE(db((pow(x, integer(2)) + integer(1)) * g));  // mixed polynomial
    REQUIRE(db(pow(x, integer(2)) * exp(pow(x, integer(2)))));  // positive c → erfi
    // The known closed form for the canonical second moment.
    REQUIRE(oracle.equivalent(integrate(pow(x, integer(2)) * g, x)->str(),
                              "-x*exp(-x**2)/2 + sqrt(pi)*erf(x)/4"));
    // Improper second moment: ∫₀^∞ x²·exp(−x²) = √π/4.
    REQUIRE(oracle.equivalent(
        integrate(pow(x, integer(2)) * g, x, S::Zero(), S::Infinity())->str(),
        "sqrt(pi)/4"));
}

// INT-GAUSSSHIFT-1: a Gaussian with a linear term ∫P(x)·exp(a·x²+b·x+c), b ≠ 0.
// Completing the square via u = x + b/(2a) reduces it to a pure Gaussian/moment in
// u. Without it ∫exp(−(x−1)²), ∫exp(−x²+x), ∫x·exp(−(x−1)²) were unevaluated.
TEST_CASE("integrate: Gaussian with linear term via completing the square (INT-GAUSSSHIFT-1)",
          "[7][integrate][erf][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        return oracle.equivalent(diff(F, x)->str(), e->str());
    };
    // exp(−(x−1)²) = exp(−x²+2x−1).
    REQUIRE(db(exp(mul(S::NegativeOne(), pow(x - integer(1), integer(2))))));
    // exp(−x²+x): a non-monic shift x₀ = 1/2, constant e^(1/4).
    REQUIRE(db(exp(mul(S::NegativeOne(), pow(x, integer(2))) + x)));
    // x·exp(−(x−1)²): polynomial residual under the shift → erf moment.
    REQUIRE(db(x * exp(mul(S::NegativeOne(), pow(x - integer(1), integer(2))))));
    // Positive leading coefficient → erfi.
    REQUIRE(db(exp(pow(x, integer(2)) + x)));
    // General a≠1: exp(−2x²+3x−1).
    REQUIRE(db(exp(mul(integer(-2), pow(x, integer(2))) + integer(3) * x
                   - integer(1))));
}

// INT-GAUSS-PARAM-1: the Gaussian with a symbolic positive coefficient —
// ∫exp(−a·x²) = √π·erf(√a·x)/(2√a), ∫exp(a·x²) = √π·erfi(√a·x)/(2√a). try_gaussian
// already branched on is_negative/is_positive; a leftover rational-only gate
// blocked symbolic coefficients (now relaxed). Matches SymPy exactly.
TEST_CASE("integrate: parametric Gaussian → erf/erfi (INT-GAUSS-PARAM-1)",
          "[7][integrate][erf][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto a = symbol("a", AssumptionMask{}.set_positive(true));
    REQUIRE(oracle.equivalent(
        integrate(exp(mul(mul(S::NegativeOne(), a), pow(x, integer(2)))), x)
            ->str(),
        "sqrt(pi)*erf(sqrt(a)*x)/(2*sqrt(a))"));
    REQUIRE(oracle.equivalent(
        integrate(exp(a * pow(x, integer(2))), x)->str(),
        "sqrt(pi)*erfi(sqrt(a)*x)/(2*sqrt(a))"));
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

// INT-PROD2SUM-1: product-to-sum for sin·sin and cos·cos (not just sin·cos), and
// for three-way products (reduced pair by pair). ∫sin(2x)·sin(3x), ∫cos·cos·cos,
// ∫x·sin(2x)·cos(3x) close via the identities + linearity. Verified by diff-back.
TEST_CASE("integrate: product-to-sum sin·sin, cos·cos, triple (INT-PROD2SUM-1)",
          "[7][integrate][trig_reduction][regression]") {
    auto x = symbol("x");
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        // diff-back numerically (some forms keep cos² that SymPy's simplify
        // won't reduce against the integrand symbolically).
        Expr resid = diff(F, x) - e;
        for (Expr pt : {rational(2, 5), rational(7, 4), rational(11, 3)}) {
            double d = std::stod(evalf(subs(resid, x, pt))->str());
            if (std::fabs(d) > 1e-7) return false;
        }
        return true;
    };
    REQUIRE(db(sin(integer(2) * x) * sin(integer(3) * x)));   // sin·sin
    REQUIRE(db(cos(integer(2) * x) * cos(integer(5) * x)));   // cos·cos
    REQUIRE(db(cos(x) * cos(integer(2) * x) * cos(integer(3) * x)));  // triple cos
    REQUIRE(db(sin(x) * sin(integer(2) * x) * sin(integer(3) * x)));  // triple sin
    REQUIRE(db(x * sin(integer(2) * x) * cos(integer(3) * x)));  // with poly factor
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

// INT-RECIPTRIG-PARTS-1: ∫x·sec²(x) = ∫x/cos²(x) and the reciprocal-square trig /
// hyperbolic family, by parts with v = ∫target tabulated (∫1/cos² = tan, etc.).
// The by-parts target whitelist was widened to negative integer powers; a
// recursive marker check bails on odd reciprocal powers (∫x/cos x is
// non-elementary) instead of returning a partial Integral(...) garbage.
TEST_CASE("integrate: polynomial × reciprocal-square trig by parts (INT-RECIPTRIG-PARTS-1)",
          "[7][integrate][parts][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        return oracle.equivalent(diff(F, x)->str(), e->str());
    };
    REQUIRE(db(x / pow(cos(x), integer(2))));   // ∫x·sec²x = x·tan x + log cos x
    REQUIRE(db(x / pow(sin(x), integer(2))));   // ∫x·csc²x = −x·cot x + log sin x
    REQUIRE(db(x / pow(cosh(x), integer(2))));  // ∫x·sech²x = x·tanh x − log cosh x
    REQUIRE(db(x / pow(sinh(x), integer(2))));  // ∫x·csch²x
    // Non-elementary odd reciprocal power (∫x/cos x) must stay an unevaluated
    // marker, not a partial garbage answer with embedded Integral(...) terms.
    REQUIRE(integrate(x / cos(x), x)->str().find("Integral(") != std::string::npos);
    // Positive-power and bare cases are unchanged.
    REQUIRE(oracle.equivalent(diff(integrate(x * pow(cos(x), integer(2)), x),
                                   x)->str(),
                              (x * pow(cos(x), integer(2)))->str()));
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

// INT-EXPINT-POWER-1: ∫f(c·x)/xⁿ for f ∈ {sin,cos,exp,sinh,cosh}, n ≥ 2 — by
// parts down to the n = 1 Si/Ci/Ei base case. Previously only n = 1 was handled.
TEST_CASE("integrate: ∫f(c·x)/x^n reduces to Si/Ci/Ei (INT-EXPINT-POWER-1)",
          "[7][integrate][expint][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto db = [&](Expr e) {
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        return oracle.equivalent(diff(F, x)->str(), e->str());
    };
    REQUIRE(db(sin(x) / pow(x, integer(2))));            // = Ci(x) − sin(x)/x
    REQUIRE(db(cos(x) / pow(x, integer(2))));            // = −Si(x) − cos(x)/x
    REQUIRE(db(exp(x) / pow(x, integer(2))));            // = Ei(x) − exp(x)/x
    REQUIRE(db(sin(x) / pow(x, integer(3))));            // two-step reduction
    REQUIRE(db(sinh(x) / pow(x, integer(2))));           // hyperbolic → Chi
    REQUIRE(db(sin(integer(2) * x) / pow(x, integer(2))));  // scaled argument
    // n = 1 base cases still close to the bare special-integral functions.
    REQUIRE(oracle.equivalent(integrate(sin(x) / x, x)->str(), "Si(x)"));
    REQUIRE(oracle.equivalent(integrate(cos(x) / x, x)->str(), "Ci(x)"));
    // Squared trig over a power: the half-angle identity reduces sin²/cos² to a
    // (1∓cos 2x)/2 form, then the linear + Si/Ci paths close it.
    REQUIRE(db(pow(sin(x), integer(2)) / pow(x, integer(2))));  // Si(2x)+…
    REQUIRE(db(pow(cos(x), integer(2)) / pow(x, integer(2))));
    REQUIRE(db(pow(sin(x), integer(2)) / x));                  // Ci-form
    REQUIRE(db(pow(cos(x), integer(2)) / x));
    REQUIRE(db(pow(sin(x), integer(2)) / pow(x, integer(3))));
    // A pure trig × trig product is NOT hijacked by the half-angle rewrite — it
    // keeps the clean sin^m·cos^n closed form.
    REQUIRE(oracle.equivalent(
        integrate(pow(sin(x), integer(3)) * pow(cos(x), integer(2)), x)->str(),
        "cos(x)**5/5 - cos(x)**3/3"));
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
    // INT-INVTRIG-RECIP-1: a reciprocal-power dv = 1/xⁿ for the rational-derivative
    // inverse functions (atan/acot/atanh) — the residual stays rational, so it
    // closes where the polynomial-only gate previously bailed.
    REQUIRE(db(atan(x) / pow(x, integer(2))));    // ∫atan(x)/x²
    REQUIRE(db(atan(x) / pow(x, integer(3))));    // ∫atan(x)/x³
    REQUIRE(db(acot(x) / pow(x, integer(2))));    // ∫acot(x)/x²
    REQUIRE(db(atanh(x) / pow(x, integer(2))));   // ∫atanh(x)/x²
    // ∫atan(x)/x is genuinely non-elementary (residual log(x)/(x²+1)) — must stay
    // an unevaluated marker, not a wrong closed form.
    REQUIRE(integrate(atan(x) / x, x)->str().find("Integral(")
            != std::string::npos);
    // The √-derivative functions over 1/xⁿ are NOT admitted (the residual is
    // non-rational and the rational path mis-handles it) — they bail to a marker.
    REQUIRE(integrate(asin(x) / pow(x, integer(2)), x)->str().find("Integral(")
            != std::string::npos);
    // INT-INVTRIG-SQ-1: x·f(x)² for a rational-derivative inverse function. The
    // by-parts target now admits a positive integer power f^k (recursing down a
    // power each step), and a rational dv keeps the residual rational so it closes.
    // ∫x·atan(x)² = x²·atan²/2 − x·atan + atan²/2 + log(x²+1)/2.
    REQUIRE(db(x * pow(atan(x), integer(2))));
    REQUIRE(db(x * pow(acot(x), integer(2))));
    // A non-elementary case (no polynomial factor) stays an unevaluated marker.
    REQUIRE(integrate(pow(atan(x), integer(2)), x)->str().find("Integral(")
            != std::string::npos);
    // INT-INVTRIG-SQRT-SQ-1: the √-derivative inverse functions f² (asin/acos/
    // asinh/acosh). A dv = P(x)/√(quadratic) is admitted, and a numeric diff-back
    // guards the broadened recursion (it once produced a wrong factor before the
    // try_x_over_sqrt_quadratic coefficient bug was fixed). ∫asin² =
    // x·asin² − 2x + 2√(1−x²)·asin.
    REQUIRE(db(pow(asin(x), integer(2))));        // bare asin²
    REQUIRE(db(pow(acos(x), integer(2))));
    REQUIRE(db(pow(asinh(x), integer(2))));
    REQUIRE(db(x * pow(asin(x), integer(2))));     // x·asin²
    REQUIRE(db(pow(asin(x), integer(3))));         // asin³
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
    // Degenerate a=b for cosine: ∫1/(1+cos x) = tan(x/2). together() leaves the
    // half-angle denominator un-cancelled here (it collapses to a constant only
    // after full cancellation), which previously made try_rational emit zoo·log 2;
    // cancel() reduces it first.
    REQUIRE(db(pow(integer(1) + cos(x), integer(-1))));               // a=b: tan(x/2)
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
    // Numerator-bearing rational trig: cos(x)/(1+cos x) = x − tan(x/2), returned
    // as −tan(x/2) + 2·atan(tan x/2). The half-angle substitution leaves a nested
    // fraction in the denominator (1 + (1−t²)/(1+t²)); flatten_ratio must reduce
    // it so the integral closes instead of bailing. SymPy's simplify can't reduce
    // the derivative form symbolically, so diff-back is verified numerically.
    {
        Expr e = cos(x) * pow(integer(1) + cos(x), integer(-1));
        Expr F = integrate(e, x);
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        Expr residual = subs(diff(F, x) - e, x, rational(1, 2));
        REQUIRE(std::abs(std::stod(oracle.evalf(residual->str()))) < 1e-9);
    }
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

// INT-POLYPROD-1: a product or power of polynomials — x³·(1−x)², x²·(1−x) — was
// returned unevaluated (no handler expands a product of polynomial factors), so
// the definite integral garbled. integrate() now expands such a product into a
// sum and integrates term-wise as a last resort. Matches SymPy.
TEST_CASE("integrate: products of polynomials (INT-POLYPROD-1)",
          "[7][integrate][regression]") {
    auto x = symbol("x");
    auto one = integer(1);

    // Indefinite: x³·(1−x)² now integrates instead of returning the marker; its
    // derivative recovers the (expanded) integrand.
    {
        auto F = integrate(pow(x, integer(3)) * pow(one - x, integer(2)), x);
        REQUIRE(F->str().find("Integral") == std::string::npos);
        REQUIRE(simplify(diff(F, x)
                         - (pow(x, integer(3)) * pow(one - x, integer(2))))
                == S::Zero());
    }

    // Definite Beta-like integrals: ∫₀¹ x³(1−x)² = ∫₀¹ x²(1−x)³ = 1/60.
    REQUIRE(integrate(pow(x, integer(3)) * pow(one - x, integer(2)), x,
                      S::Zero(), one)
            == rational(1, 60));
    REQUIRE(integrate(pow(x, integer(2)) * pow(one - x, integer(3)), x,
                      S::Zero(), one)
            == rational(1, 60));
    // ∫₀¹ x²(1−x) = 1/12, ∫₀¹ x(x+1)(x+2) = 9/4, ∫₀¹ (2x−1)³ = 0.
    REQUIRE(integrate(pow(x, integer(2)) * (one - x), x, S::Zero(), one)
            == rational(1, 12));
    REQUIRE(integrate(x * (x + one) * (x + integer(2)), x, S::Zero(), one)
            == rational(9, 4));
    REQUIRE(integrate(pow(integer(2) * x - one, integer(3)), x, S::Zero(), one)
            == S::Zero());

    // By-parts and rational integrands are unaffected by the expand fallback.
    REQUIRE(integrate(x * sin(x), x) == sin(x) - x * cos(x));
    REQUIRE(integrate(pow(x * x + one, integer(-1)), x) == atan(x));
}

// INT-DIRICHLET-1 / INT-CLEANMARKER-1: ∫₀^∞ (1−cos x)/x² = π/2 (a Dirichlet
// integral). Its antiderivative Si(x) + cos(x)/x − 1/x is in a factored form
// −1·(−Si(x) − cos(x)/x) − 1/x whose c·(sum) limit folded to a wrong 0; with the
// constant-factor pull-out in the limit engine it resolves. Separately, a
// genuinely unintegrable definite integral now returns a clean unevaluated
// Integral instead of garbled −Integral(nan,a)+Integral(…,b). Matches SymPy.
TEST_CASE("integrate: Dirichlet (1-cos x)/x² and clean markers (INT-DIRICHLET-1)",
          "[7][integrate][definite][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto oo = S::Infinity();

    // ∫₀^∞ (1−cos x)/x² = π/2.
    REQUIRE(oracle.equivalent(
        integrate((integer(1) - cos(x)) * pow(x, integer(-2)), x, S::Zero(), oo)
            ->str(),
        "pi/2"));

    // A definite integral with no closed form returns a clean unevaluated
    // Integral marker — not garbage with embedded nan. (x/(e^x+1) was the old
    // example here, but it is now evaluated in closed form by the Bose–Einstein /
    // Fermi–Dirac rule, see INT-BOSE-EINSTEIN-1; x/(x+e^x) stays non-elementary.)
    auto uneval =
        integrate(x * pow(x + exp(x), integer(-1)), x, S::Zero(), oo);
    REQUIRE(uneval->str().find("nan") == std::string::npos);
    REQUIRE(uneval->str().rfind("Integral(", 0) == 0);

    // Previously-working improper integrals are unchanged.
    REQUIRE(oracle.equivalent(
        integrate(pow(sin(x), integer(2)) * pow(x, integer(-2)), x, S::Zero(),
                  oo)
            ->str(),
        "pi/2"));
}

// INT-PARAMSIGN-1: parametric improper integrals ∫₀^∞ f(x; a) where the
// boundary limit depends on the sign of a symbolic coefficient. With a declared
// positive, exp(-a·x)→0 at +∞ and the closed forms resolve (1/a, 1/a², …). With
// a *generic* a the boundary limit cannot decide the sign, so Newton-Leibniz
// previously leaked garbage like (−exp(a·−oo)+1)/a — a raw ±oo buried inside
// exp/Si. integrate() now detects an infinity stranded below the root and falls
// back to a clean unevaluated Integral, matching SymPy (which returns a
// Piecewise gated on sign(a); the unevaluated form is its "otherwise" branch).
TEST_CASE("integrate: parametric improper integrals by parameter sign (INT-PARAMSIGN-1)",
          "[7][integrate][definite][regression]") {
    auto x = symbol("x");
    auto ap = symbol("a", AssumptionMask{}.set_positive(true));
    auto ag = symbol("a");  // sign unknown
    auto oo = S::Infinity();
    auto z = S::Zero();
    auto negax = [&](const Expr& a) { return mul(mul(S::NegativeOne(), a), x); };

    // Positive a: closed forms resolve and contain no stray infinity.
    REQUIRE(simplify(integrate(exp(negax(ap)), x, z, oo))
            == pow(ap, integer(-1)));
    REQUIRE(simplify(integrate(x * exp(negax(ap)), x, z, oo))
            == pow(ap, integer(-2)));
    REQUIRE(simplify(integrate(pow(x, integer(2)) * exp(negax(ap)), x, z, oo))
            == integer(2) * pow(ap, integer(-3)));

    // Generic a: clean unevaluated marker, never garbage with an embedded oo.
    for (const Expr& f : {exp(negax(ag)), x * exp(negax(ag)),
                          sin(mul(ag, x)) * pow(x, integer(-1))}) {
        auto r = integrate(f, x, z, oo);
        REQUIRE(r->str().rfind("Integral(", 0) == 0);
        // The only infinity is the bound "oo"; no stray "-oo" leaked from a
        // boundary substitution like exp(a·-oo).
        REQUIRE(r->str().find("-oo") == std::string::npos);
    }
}

// INT-GAUSSFOURIER-1: the Fourier integral of a real Gaussian. The integrand
// exp(-a x²)·cos(b x) has no elementary antiderivative, so Newton-Leibniz
// garbled it; the closed form is ∫_{-∞}^{∞} = √(π/a)·exp(-b²/(4a)), half that
// over [0,∞), and the sine version is 0 over the symmetric line. Matches SymPy
// (which actually crashes on the sine case — SymPP returns the correct 0).
TEST_CASE("integrate: Gaussian Fourier integrals (INT-GAUSSFOURIER-1)",
          "[7][integrate][definite][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto oo = S::Infinity();
    auto noo = S::NegativeInfinity();
    auto gauss = [&](const Expr& a) { return exp(mul(a, pow(x, integer(2)))); };

    // ∫_{-∞}^{∞} exp(-x²)·cos(x) = √π·exp(-1/4).
    REQUIRE(oracle.equivalent(
        integrate(gauss(integer(-1)) * cos(x), x, noo, oo)->str(),
        "sqrt(pi)*exp(-1/4)"));
    // Half line: ½·√π·exp(-1/4).
    REQUIRE(oracle.equivalent(
        integrate(gauss(integer(-1)) * cos(x), x, S::Zero(), oo)->str(),
        "sqrt(pi)*exp(-1/4)/2"));
    // Scaled: ∫_{-∞}^{∞} exp(-2x²)·cos(3x) = √(π/2)·exp(-9/8).
    REQUIRE(oracle.equivalent(
        integrate(gauss(integer(-2)) * cos(integer(3) * x), x, noo, oo)->str(),
        "sqrt(2)*sqrt(pi)*exp(-Rational(9,8))/2"));
    // ∫_{-∞}^{∞} exp(-x²)·cos(2x) = √π·exp(-1).
    REQUIRE(oracle.equivalent(
        integrate(gauss(integer(-1)) * cos(integer(2) * x), x, noo, oo)->str(),
        "sqrt(pi)*exp(-1)"));
    // Constant prefactor carries through.
    REQUIRE(oracle.equivalent(
        integrate(integer(5) * gauss(integer(-1)) * cos(x), x, noo, oo)->str(),
        "5*sqrt(pi)*exp(-1/4)"));
    // Odd integrand over the symmetric line integrates to 0.
    REQUIRE(integrate(gauss(integer(-1)) * sin(x), x, noo, oo) == S::Zero());
    // The pure Gaussian (no trig factor) is unaffected.
    REQUIRE(oracle.equivalent(integrate(gauss(integer(-1)), x, noo, oo)->str(),
                              "sqrt(pi)"));
}

// INT-POLY-EXP-TRIG-1: ∫₀^∞ xⁿ·e^(−ax)·sin/cos(bx) — the antiderivative
// (poly × exp × trig) was correct, but the upper-bound limit returned nan
// because the engine could not see that a bounded oscillation times a decaying
// envelope vanishes (x·cos x·e^(−x) → 0). With that limit rule in place these
// close. Matches SymPy.
TEST_CASE("integrate: polynomial × decaying exp × trig (INT-POLY-EXP-TRIG-1)",
          "[7][integrate][definite][regression]") {
    auto x = symbol("x");
    auto oo = S::Infinity();
    auto edecay = exp(mul(S::NegativeOne(), x));        // e^(−x)
    auto edecay2 = exp(mul(integer(-2), x));            // e^(−2x)

    // ∫₀^∞ x·e^(−x)·sin x = 1/2, x²·e^(−x)·sin x = 1/2, x·e^(−x)·cos x = 0.
    REQUIRE(integrate(x * edecay * sin(x), x, S::Zero(), oo) == rational(1, 2));
    REQUIRE(integrate(pow(x, integer(2)) * edecay * sin(x), x, S::Zero(), oo)
            == rational(1, 2));
    REQUIRE(integrate(x * edecay * cos(x), x, S::Zero(), oo) == S::Zero());
    // Scaled exponent / frequency: ∫₀^∞ x·e^(−2x)·sin(3x) = 12/169.
    REQUIRE(integrate(x * edecay2 * sin(integer(3) * x), x, S::Zero(), oo)
            == rational(12, 169));
    // The n = 0 cases (no polynomial factor) are unchanged.
    REQUIRE(integrate(edecay * sin(x), x, S::Zero(), oo) == rational(1, 2));
}

// INT-SINSQ-1: ∫₀^∞ sin²(bx)/x² = π|b|/2 (a Dirichlet-family integral). The
// antiderivative is correct but factored — −½·(−2·Si(2x) − cos(2x)/x) − 1/(2x) —
// which hid the bounded Si inside a product so the boundary limits folded to
// wrong values. Expanding the antiderivative before Newton–Leibniz lets the
// per-term limit rules resolve each piece. Matches SymPy.
TEST_CASE("integrate: sin²(bx)/x² over [0, oo) (INT-SINSQ-1)",
          "[7][integrate][definite][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto oo = S::Infinity();

    // ∫₀^∞ sin²x/x² = π/2.
    REQUIRE(oracle.equivalent(
        integrate(pow(sin(x), integer(2)) * pow(x, integer(-2)), x, S::Zero(),
                  oo)
            ->str(),
        "pi/2"));
    // ∫₀^∞ sin²(2x)/x² = π.
    REQUIRE(oracle.equivalent(
        integrate(pow(sin(integer(2) * x), integer(2)) * pow(x, integer(-2)), x,
                  S::Zero(), oo)
            ->str(),
        "pi"));
    // A finite-bound version still evaluates (regression guard on the expand
    // retry): ∫₁² sin²x/x² matches SymPy.
    REQUIRE(oracle.equivalent(
        integrate(pow(sin(x), integer(2)) * pow(x, integer(-2)), x, integer(1),
                  integer(2))
            ->str(),
        "Si(4) - Si(2) + cos(4)/4 - cos(2)/2 + Rational(1,4)"));
}

// INT-DEF-2: improper integrals whose antiderivative carries log and atan terms
// that individually diverge at ∞ but combine to a finite limit. The upper-bound
// evaluation (a limit at ∞) now resolves the ∞ − ∞ between logs and the
// atan(arg→∞) = π/2 terms (supported by the oo + √2 = oo arithmetic fix).
// ∫₀^∞ 1/(1+x⁴) = π√2/4 used to come back as nan.
TEST_CASE("integrate: improper integrals with log/atan antiderivatives (INT-DEF-2)",
          "[7][integrate][definite][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto oo = S::Infinity();
    // ∫₀^∞ 1/(1+x⁴) = π√2/4.
    REQUIRE(oracle.equivalent(
        integrate(pow(integer(1) + pow(x, integer(4)), integer(-1)), x,
                  S::Zero(), oo)
            ->str(),
        "sqrt(2)*pi/4"));
    // ∫₀^∞ 1/(x⁴+x²+1) = π√3/6.
    REQUIRE(oracle.equivalent(
        integrate(pow(pow(x, integer(4)) + pow(x, integer(2)) + integer(1),
                      integer(-1)),
                  x, S::Zero(), oo)
            ->str(),
        "sqrt(3)*pi/6"));
    // ∫₁^∞ 1/(x(x+1)) = log 2 (a pure ∞ − ∞ between logs).
    REQUIRE(oracle.equivalent(
        integrate(pow(x * (x + integer(1)), integer(-1)), x, integer(1), oo)
            ->str(),
        "log(2)"));
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

// INT-MONOMIAL-SUB-1: monomial substitution u = x^d closes integrands of the
// form x^(d-1)·f(x^d) — e.g. x/(x⁴+1) → ½atan(x²), x³/(x⁸+1) → ¼atan(x⁴),
// x/√(1−x⁴) → ½asin(x²). The structural u-substitution (try_heurisch) misses
// these because x⁴ does not contain x² as a subexpression. Verified by
// differentiating the antiderivative back to the integrand.
TEST_CASE("integrate: monomial substitution u = x^d (INT-MONOMIAL-SUB-1)",
          "[7][integrate][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    const std::vector<Expr> integrands = {
        x * pow(pow(x, integer(4)) + integer(1), integer(-1)),       // x/(x⁴+1)
        pow(x, integer(3))
            * pow(pow(x, integer(8)) + integer(1), integer(-1)),     // x³/(x⁸+1)
        x * pow(integer(1) - pow(x, integer(4)), rational(-1, 2)),   // x/√(1−x⁴)
    };
    for (const Expr& e : integrands) {
        auto F = integrate(e, x);
        INFO("integrand: " << e->str() << "  antiderivative: " << F->str());
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
    }
    // Explicit closed forms.
    REQUIRE(oracle.equivalent(
        integrate(x * pow(pow(x, integer(4)) + integer(1), integer(-1)), x)
            ->str(),
        "atan(x**2)/2"));
    // A plain x/(x²+1) is unaffected (handled before, log form, not atan(x)).
    REQUIRE(oracle.equivalent(
        integrate(x * pow(pow(x, integer(2)) + integer(1), integer(-1)), x)
            ->str(),
        "log(x**2 + 1)/2"));
}

// INT-RATIONAL-NOPARTIAL-1: try_rational must not return a half-answer with a
// leaked Integral marker. When apart() can't fully decompose the integrand it
// now bails, so a cleaner strategy can take over (∫x²/(x⁶+1) closes to
// ⅓atan(x³) via monomial substitution) or the whole integral is returned
// honestly unevaluated (∫1/(x⁶+1)) — never `…atan(x) + Integral(…)`.
TEST_CASE("integrate: no partial rational results (INT-RATIONAL-NOPARTIAL-1)",
          "[7][integrate][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // x²/(x⁶+1) now closes cleanly.
    auto F = integrate(
        pow(x, integer(2)) * pow(pow(x, integer(6)) + integer(1), integer(-1)),
        x);
    REQUIRE(F->str().find("Integral(") == std::string::npos);
    REQUIRE(oracle.equivalent(F->str(), "atan(x**3)/3"));
    // No leaked partials: any remaining Integral must be the whole result, not a
    // summand alongside closed-form terms.
    auto no_partial = [](const Expr& r) {
        const std::string s = r->str();
        const auto pos = s.find("Integral(");
        return pos == std::string::npos || pos == 0;
    };
    for (const Expr& e : {pow(pow(x, integer(6)) + integer(1), integer(-1)),
                          pow(pow(x, integer(5)) + integer(1), integer(-1))}) {
        auto r = integrate(e, x);
        INFO("integrand: " << e->str() << "  result: " << r->str());
        REQUIRE(no_partial(r));
    }
    // Fully-solvable rationals are unaffected.
    REQUIRE(oracle.equivalent(
        diff(integrate(pow(pow(x, integer(4)) - integer(1), integer(-1)), x), x)
            ->str(),
        "1/(x**4 - 1)"));
}

// INT-RADICAL-SUB-1: integrands that are functions of x^(1/d) close via the
// substitution u = x^(1/d) (x = u^d, dx = d·u^(d-1) du): ∫exp(√x), ∫1/(1+√x),
// ∫1/(1+x^(1/3)). Previously unevaluated. Verified by differentiating back.
TEST_CASE("integrate: radical substitution u = x^(1/d) (INT-RADICAL-SUB-1)",
          "[7][integrate][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    const std::vector<Expr> integrands = {
        exp(pow(x, rational(1, 2))),                              // exp(√x)
        sin(pow(x, rational(1, 2))),                              // sin(√x)
        pow(integer(1) + pow(x, rational(1, 2)), integer(-1)),   // 1/(1+√x)
        pow(integer(1) + pow(x, rational(1, 3)), integer(-1)),   // 1/(1+x^(1/3))
        pow(pow(x, rational(1, 2)) + x, integer(-1)),            // 1/(√x + x)
    };
    for (const Expr& e : integrands) {
        auto F = integrate(e, x);
        INFO("integrand: " << e->str() << "  antiderivative: " << F->str());
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
    }
    // Explicit closed form: ∫exp(√x) = 2√x·exp(√x) − 2exp(√x).
    REQUIRE(oracle.equivalent(
        integrate(exp(pow(x, rational(1, 2))), x)->str(),
        "2*sqrt(x)*exp(sqrt(x)) - 2*exp(sqrt(x))"));
    // A plain power is still handled by the power rule (not this path).
    REQUIRE(oracle.equivalent(integrate(pow(x, rational(1, 2)), x)->str(),
                              "2*x**(3/2)/3"));
}

// INT-LINRADICAL-SUB-1: integrands containing √(a·x+b) (non-trivial linear
// inner) close via u = √(a·x+b): ∫1/(x√(x+1)), ∫√(x+1)/x. Previously
// unevaluated. Verified by differentiating back to the integrand.
TEST_CASE("integrate: linear-radical substitution √(ax+b) (INT-LINRADICAL-SUB-1)",
          "[7][integrate][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto sx1 = pow(x + integer(1), rational(1, 2));  // √(x+1)
    const std::vector<Expr> integrands = {
        pow(x * sx1, integer(-1)),                       // 1/(x·√(x+1))
        sx1 * pow(x, integer(-1)),                       // √(x+1)/x
        pow(sx1 + integer(1), integer(-1)),              // 1/(√(x+1)+1)
        x * sx1,                                         // x·√(x+1)
    };
    for (const Expr& e : integrands) {
        auto F = integrate(e, x);
        INFO("integrand: " << e->str() << "  antiderivative: " << F->str());
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
    }
    // A bare √(x+1) (no other var dependence) is the power rule's job.
    REQUIRE(oracle.equivalent(integrate(sx1, x)->str(),
                              "2*(x + 1)**(3/2)/3"));
}

// INT-EXP-SUB-1: an integrand rational in eˣ closes via u = eˣ (each
// exp(k·x+d) → e^d·uᵏ, dx = du/u), e.g. ∫1/(eˣ+e⁻ˣ) = atan(eˣ). Previously
// exp(2x)/exp(−x) terms (distinct nodes from exp(x)) blocked the substitution.
// Verified by differentiating back.
TEST_CASE("integrate: exponential substitution u=eˣ (INT-EXP-SUB-1)",
          "[7][integrate][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto ex = [&](const Expr& a) { return exp(a); };
    const std::vector<Expr> integrands = {
        pow(ex(x) + ex(-x), integer(-1)),                       // 1/(eˣ+e⁻ˣ)
        ex(x) * pow(ex(integer(2) * x) + integer(1), integer(-1)),  // eˣ/(e²ˣ+1)
        ex(integer(2) * x) * pow(integer(1) + ex(x), integer(-1)),  // e²ˣ/(1+eˣ)
        pow(ex(x) + ex(integer(2) * x), integer(-1)),           // 1/(eˣ+e²ˣ)
    };
    for (const Expr& e : integrands) {
        auto F = integrate(e, x);
        INFO("integrand: " << e->str() << "  antiderivative: " << F->str());
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
    }
    // The headline closed form matches SymPy exactly.
    REQUIRE(oracle.equivalent(
        integrate(pow(ex(x) + ex(-x), integer(-1)), x)->str(), "atan(exp(x))"));
}

// INT-RECIP-SUB-1: ∫dx/(xⁿ·√(a·x²+c)) closes via the reciprocal substitution
// x = 1/u: substituting clears the x in the denominator, and pulling u out of
// the radical ((a·u⁻²+c)^e = u^(−2e)·(a+c·u²)^e) leaves an ordinary
// √(quadratic) integral. The antiderivative matches SymPy's (it carries the
// same x>0 principal-branch convention, so diff-back is checked on x>0).
TEST_CASE("integrate: reciprocal substitution x=1/u (INT-RECIP-SUB-1)",
          "[7][integrate][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto recip_sqrt = [&](const Expr& xn, const Expr& a, const Expr& c) {
        // 1 / (xn · √(a·x² + c))
        return pow(xn * pow(a * pow(x, integer(2)) + c, rational(1, 2)),
                   integer(-1));
    };
    // ∫1/(x√(x²+1)) = −asinh(1/x), ∫1/(x√(x²+4)) = −asinh(2/x)/2.
    REQUIRE(oracle.equivalent(
        integrate(recip_sqrt(x, integer(1), integer(1)), x)->str(),
        "-asinh(1/x)"));
    REQUIRE(oracle.equivalent(
        integrate(recip_sqrt(x, integer(1), integer(4)), x)->str(),
        "-asinh(2/x)/2"));
    REQUIRE(oracle.equivalent(
        integrate(recip_sqrt(x, integer(9), integer(1)), x)->str(),
        "-asinh(1/(3*x))"));
    // Higher denominator powers also close (form differs from SymPy by the
    // x>0 |x| branch, so just require a closed form — no Integral marker).
    for (int n : {2, 3}) {
        auto F = integrate(recip_sqrt(pow(x, integer(n)), integer(1),
                                      integer(1)),
                           x);
        INFO("n=" << n << " antiderivative: " << F->str());
        REQUIRE(F->str().find("Integral(") == std::string::npos);
    }
}

// INT-QUAD-IRRATIONAL-1: ∫1/(a·x²+b·x+c) where the discriminant is positive but
// the roots are irrational (no rational factorization) — 1/(x²−3), 1/(2x²−3),
// 1/(x²+x−1). The partial-fraction logs carry √Δ. Previously unevaluated (the
// rational path only factors over ℚ; the arctan path only handled Δ<0).
TEST_CASE("integrate: 1/quadratic with irrational roots (INT-QUAD-IRRATIONAL-1)",
          "[7][integrate][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto recip_quad = [&](const Expr& q) { return pow(q, integer(-1)); };
    const std::vector<Expr> integrands = {
        recip_quad(pow(x, integer(2)) - integer(3)),                // 1/(x²−3)
        recip_quad(integer(3) - pow(x, integer(2))),                // 1/(3−x²)
        recip_quad(integer(2) * pow(x, integer(2)) - integer(3)),   // 1/(2x²−3)
        recip_quad(pow(x, integer(2)) + x - integer(1)),            // 1/(x²+x−1)
    };
    for (const Expr& e : integrands) {
        auto F = integrate(e, x);
        INFO("integrand: " << e->str() << "  antiderivative: " << F->str());
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
    }
    // Rational-root and negative-discriminant cases are unchanged.
    REQUIRE(oracle.equivalent(
        integrate(pow(pow(x, integer(2)) - integer(1), integer(-1)), x)->str(),
        "log(x - 1)/2 - log(x + 1)/2"));
    REQUIRE(oracle.equivalent(
        integrate(pow(pow(x, integer(2)) + integer(1), integer(-1)), x)->str(),
        "atan(x)"));
}

// INT-BIQUADRATIC-1: ∫1/(a₄x⁴+a₂x²+a₀) for a biquadratic that is irreducible
// over ℚ but factors into two real quadratics over ℚ(√q) — the canonical
// ∫1/(x⁴+1). Closed via that factorization and partial fractions; previously
// unevaluated. Verified by differentiating back.
TEST_CASE("integrate: irreducible biquadratic 1/(x⁴+q) (INT-BIQUADRATIC-1)",
          "[7][integrate][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto recip4 = [&](const Expr& num4, const Expr& c) {
        return pow(num4 * pow(x, integer(4)) + c, integer(-1));
    };
    const std::vector<Expr> integrands = {
        recip4(integer(1), integer(1)),   // 1/(x⁴+1)
        recip4(integer(1), integer(4)),   // 1/(x⁴+4)
        recip4(integer(2), integer(2)),   // 1/(2x⁴+2)
        recip4(integer(1), integer(9)),   // 1/(x⁴+9)
    };
    for (const Expr& e : integrands) {
        auto F = integrate(e, x);
        INFO("integrand: " << e->str() << "  antiderivative: " << F->str());
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
    }
    // Biquadratics that factor over ℚ are still handled by the rational path.
    REQUIRE(oracle.equivalent(
        diff(integrate(pow(pow(x, integer(4)) + pow(x, integer(2)) + integer(1),
                           integer(-1)),
                       x),
             x)
            ->str(),
        "1/(x**4 + x**2 + 1)"));
}

// INT-BIQUAD-NUM-1: ∫(n₂x²+n₀)/(a₄x⁴+a₂x²+a₀) — an even polynomial numerator
// over an irreducible biquadratic, e.g. ∫x²/(x⁴+1). Closed via the same
// ℚ(√q) factorization as INT-BIQUADRATIC-1 with the numerator distributed
// across the two real-quadratic partial fractions; previously unevaluated.
// Verified by differentiating back.
TEST_CASE("integrate: even numerator over biquadratic (INT-BIQUAD-NUM-1)",
          "[7][integrate][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto over4 = [&](const Expr& numer, const Expr& a2, const Expr& a0) {
        return numer
               * pow(pow(x, integer(4)) + a2 * pow(x, integer(2)) + a0,
                     integer(-1));
    };
    const std::vector<Expr> integrands = {
        over4(pow(x, integer(2)), integer(0), integer(1)),                 // x²/(x⁴+1)
        over4(pow(x, integer(2)) + integer(1), integer(0), integer(1)),    // (x²+1)/(x⁴+1)
        over4(pow(x, integer(2)), integer(1), integer(1)),                 // x²/(x⁴+x²+1)
        over4(integer(3) * pow(x, integer(2)) + integer(2), integer(0),
              integer(4)),                                    // (3x²+2)/(x⁴+4)
    };
    for (const Expr& e : integrands) {
        auto F = integrate(e, x);
        INFO("integrand: " << e->str() << "  antiderivative: " << F->str());
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
    }
}

// INT-SINCOS-QUOT-1: ∫ sinᵐ(x)cosⁿ(x) dx for a sin/cos quotient (one negative
// power) with at least one ODD exponent, e.g. ∫cos²/sin, ∫sin³/cos², ∫cos²/sin³.
// The substitution u=sin(x) (cos odd) or u=cos(x) (sin odd) makes the integrand
// rational in u, which the rational path closes; previously all unevaluated.
// Verified by differentiating back.
TEST_CASE("integrate: sin/cos quotients with an odd power (INT-SINCOS-QUOT-1)",
          "[7][integrate][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto sc = [&](long m, long n) {
        return pow(sin(x), integer(m)) * pow(cos(x), integer(n));
    };
    const std::vector<Expr> integrands = {
        sc(-1, 2),  // cos²/sin
        sc(2, -1),  // sin²/cos
        sc(3, -1),  // sin³/cos
        sc(-1, 3),  // cos³/sin
        sc(3, -2),  // sin³/cos²
        sc(-3, 2),  // cos²/sin³
    };
    for (const Expr& e : integrands) {
        auto F = integrate(e, x);
        INFO("integrand: " << e->str() << "  antiderivative: " << F->str());
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
    }
}

// INT-SINCOS-QUOT-2: ∫ sinᵐ(x)cosⁿ(x) dx for a sin/cos quotient with BOTH exponents
// even, e.g. ∫cos⁴/sin², ∫cos²/sin² (=cot²), ∫sin⁴/cos². The substitution t=tan(x)
// makes the integrand rational in t, which the rational path closes; previously all
// unevaluated. Verified by differentiating back.
TEST_CASE("integrate: even/even sin/cos quotients (INT-SINCOS-QUOT-2)",
          "[7][integrate][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto sc = [&](long m, long n) {
        return pow(sin(x), integer(m)) * pow(cos(x), integer(n));
    };
    const std::vector<Expr> integrands = {
        sc(-2, 4),  // cos⁴/sin²
        sc(-2, 2),  // cos²/sin² = cot²
        sc(4, -2),  // sin⁴/cos²
        sc(-4, 2),  // cos²/sin⁴
        sc(2, -4),  // sin²/cos⁴
        sc(-2, 6),  // cos⁶/sin²
    };
    for (const Expr& e : integrands) {
        auto F = integrate(e, x);
        INFO("integrand: " << e->str() << "  antiderivative: " << F->str());
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        REQUIRE(oracle.equivalent(diff(F, x)->str(), e->str()));
    }
}

// INT-TRIGPROD-1: ∫ ∏ sin/cos(affineᵢ)^kᵢ dx — a product of sin/cos powers whose
// arguments are NOT all equal, e.g. ∫sin²(x)cos(2x), ∫sin³(x)cos(2x),
// ∫sin²(2x)cos(x). Product-to-sum linearization reduces it to a sum of single
// sin/cos terms the table integrates; previously all unevaluated. Verified by
// high-precision numeric diff-back (SymPy's symbolic simplify can't reliably
// reduce a trig product, so the standard oracle.equivalent check is unreliable).
TEST_CASE("integrate: trig products with mixed arguments (INT-TRIGPROD-1)",
          "[7][integrate][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    const std::vector<Expr> integrands = {
        pow(sin(x), integer(2)) * cos(integer(2) * x),       // sin²(x)cos(2x)
        pow(cos(x), integer(2)) * cos(integer(2) * x),       // cos²(x)cos(2x)
        pow(sin(x), integer(2)) * sin(integer(2) * x),       // sin²(x)sin(2x)
        pow(sin(x), integer(3)) * cos(integer(2) * x),       // sin³(x)cos(2x)
        pow(sin(integer(2) * x), integer(2)) * cos(x),       // sin²(2x)cos(x)
        sin(x) * cos(integer(2) * x),                        // sin(x)cos(2x)
    };
    for (const Expr& e : integrands) {
        auto F = integrate(e, x);
        INFO("integrand: " << e->str() << "  antiderivative: " << F->str());
        REQUIRE(F->str().find("Integral(") == std::string::npos);
        Expr resid = diff(F, x) - e;
        for (long num : {3, 7, 13, 19}) {
            Expr at = subs(resid, x, rational(num, 29));
            auto resp = oracle.send({{"op", "evalf_is_zero"},
                                     {"expr", at->str()},
                                     {"prec", 40},
                                     {"tol", 20}});
            REQUIRE(resp.ok);
            REQUIRE(resp.raw.at("result").get<bool>());
        }
    }
}

// INT-LOGPOW-BOUNDARY-1: ∫₀¹ (log x)ⁿ — the antiderivative x·(log x)ⁿ − … is
// known, but the lower boundary value needs lim_{x→0} x·(log x)ⁿ = 0, a 0·∞ form
// that is finite only from the RIGHT (the left side is complex). The boundary
// evaluation used a two-sided limit and got nan; it now approaches each bound
// from inside the interval (lower from the right, upper from the left). Matches
// SymPy: ∫₀¹ (log x)² = 2, (log x)³ = −6, (log x)⁴ = 24.
TEST_CASE("integrate: log-power definite integrals via one-sided boundaries (INT-LOGPOW-BOUNDARY-1)",
          "[7][integrate][definite][regression]") {
    auto x = symbol("x");
    auto zero = S::Zero();
    auto one = S::One();
    auto L = [&](const Expr& e) { return log(e); };

    REQUIRE(integrate(pow(L(x), integer(2)), x, zero, one) == integer(2));
    REQUIRE(integrate(pow(L(x), integer(3)), x, zero, one) == integer(-6));
    REQUIRE(integrate(pow(L(x), integer(4)), x, zero, one) == integer(24));
    // With a polynomial weight: ∫₀¹ x·(log x)² = 1/4, ∫₀¹ √x·log x = −4/9.
    REQUIRE(integrate(x * pow(L(x), integer(2)), x, zero, one)
            == rational(1, 4));
    REQUIRE(integrate(pow(x, rational(1, 2)) * L(x), x, zero, one)
            == rational(-4, 9));
    // ∫₀¹ log x = −1 (the n = 1 case kept working) is unchanged.
    REQUIRE(integrate(L(x), x, zero, one) == S::NegativeOne());
    // A non-singular boundary is unaffected: ∫₀¹ x² = 1/3.
    REQUIRE(integrate(pow(x, integer(2)), x, zero, one) == rational(1, 3));
}

// INT-BOSE-EINSTEIN-1: ∫₀^∞ x^p/(e^{cx}∓1) dx in closed form. Bose–Einstein
// (e^{cx}−1): Γ(p+1)·ζ(p+1)/c^{p+1} (Debye / Stefan–Boltzmann). Fermi–Dirac
// (e^{cx}+1): × (1−2^{−p}) (Dirichlet eta). The antiderivative is non-elementary,
// so SymPP previously returned an unevaluated marker (as SymPy still does). The
// closed forms were verified numerically against SymPy quadrature to 1e-9; the
// symbolic equalities below (π²/6, π⁴/15, …) are confirmed via the oracle.
TEST_CASE("integrate: Bose-Einstein / Fermi-Dirac integrals (INT-BOSE-EINSTEIN-1)",
          "[7][integrate][definite][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto zero = S::Zero();
    auto oo = S::Infinity();
    auto kernel = [&](const Expr& c, int s) {
        return pow(exp(c * x) + integer(s), integer(-1));
    };

    // Bose–Einstein, c = 1: Γ(p+1)ζ(p+1).
    REQUIRE(oracle.equivalent(
        integrate(x * kernel(S::One(), -1), x, zero, oo)->str(), "pi**2/6"));
    REQUIRE(oracle.equivalent(
        integrate(pow(x, integer(3)) * kernel(S::One(), -1), x, zero, oo)->str(),
        "pi**4/15"));
    REQUIRE(oracle.equivalent(
        integrate(pow(x, integer(2)) * kernel(S::One(), -1), x, zero, oo)->str(),
        "2*zeta(3)"));
    REQUIRE(oracle.equivalent(
        integrate(pow(x, integer(5)) * kernel(S::One(), -1), x, zero, oo)->str(),
        "8*pi**6/63"));
    // Scaled rate c = 2: Γ(3)ζ(3)/2³ = ζ(3)/4.
    REQUIRE(oracle.equivalent(
        integrate(pow(x, integer(2)) * kernel(integer(2), -1), x, zero, oo)
            ->str(),
        "zeta(3)/4"));

    // Fermi–Dirac (e^{x}+1): × (1 − 2^{−p}).
    REQUIRE(oracle.equivalent(
        integrate(x * kernel(S::One(), 1), x, zero, oo)->str(), "pi**2/12"));
    REQUIRE(oracle.equivalent(
        integrate(pow(x, integer(3)) * kernel(S::One(), 1), x, zero, oo)->str(),
        "7*pi**4/120"));

    // Guards: a finite upper bound or a non-monomial numerator must NOT fire.
    REQUIRE(integrate(x * kernel(S::One(), -1), x, zero, S::One())
                ->str()
                .rfind("Integral(", 0) == 0);
    // Plain Γ integral (no kernel) is unaffected: ∫₀^∞ x²e^{−x} = 2.
    REQUIRE(integrate(pow(x, integer(2)) * exp(mul(S::NegativeOne(), x)), x,
                      zero, oo)
            == integer(2));
}

// INT-GAMMA-LOG-1: ∫₀^∞ x^p·e^{−cx}·log x = Γ(p+1)(ψ(p+1) − log c)/c^{p+1}, the
// s-derivative of the Γ-integral. Non-elementary antiderivative, so SymPP
// previously returned nan. ψ = digamma evaluates at integer p, giving closed
// forms in EulerGamma. Matches SymPy: ∫₀^∞ e^{−x} log x = −γ, x³ variant = 11−6γ.
TEST_CASE("integrate: gamma-log integrals ∫x^p e^{-cx} log x (INT-GAMMA-LOG-1)",
          "[7][integrate][definite][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto zero = S::Zero();
    auto oo = S::Infinity();
    auto E = [&](const Expr& c) { return exp(mul(c, x)); };  // e^{c·x}

    // ∫₀^∞ e^{−x} log x = −γ.
    REQUIRE(oracle.equivalent(
        integrate(E(S::NegativeOne()) * log(x), x, zero, oo)->str(),
        "-EulerGamma"));
    // ∫₀^∞ x·e^{−x} log x = 1 − γ.
    REQUIRE(oracle.equivalent(
        integrate(x * E(S::NegativeOne()) * log(x), x, zero, oo)->str(),
        "1 - EulerGamma"));
    // ∫₀^∞ x³·e^{−x} log x = 11 − 6γ.
    REQUIRE(oracle.equivalent(
        integrate(pow(x, integer(3)) * E(S::NegativeOne()) * log(x), x, zero, oo)
            ->str(),
        "11 - 6*EulerGamma"));
    // Scaled rate c = 2: ∫₀^∞ e^{−2x} log x = −γ/2 − log(2)/2.
    REQUIRE(oracle.equivalent(
        integrate(E(integer(-2)) * log(x), x, zero, oo)->str(),
        "-EulerGamma/2 - log(2)/2"));

    // Guards: no log factor (plain Γ-integral) and a finite bound must NOT fire.
    REQUIRE(integrate(pow(x, integer(3)) * E(S::NegativeOne()), x, zero, oo)
            == integer(6));
    REQUIRE(integrate(x * E(S::NegativeOne()), x, zero, oo) == integer(1));
}

// INT-LOG-QUADRATIC-1: ∫₀^∞ log(x)/(x²+A) dx = π·log(A)/(4√A) for a positive
// constant A (equivalently π·log(a)/(2a) with A = a²). Non-elementary
// antiderivative, so SymPP returned an unevaluated marker. Gives the classic
// ∫₀^∞ log(x)/(x²+1) = 0 and ∫₀^∞ log(x)/(x²+4) = π·log(4)/8. Matches SymPy.
TEST_CASE("integrate: ∫log(x)/(x²+A) over [0,∞) (INT-LOG-QUADRATIC-1)",
          "[7][integrate][definite][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto zero = S::Zero();
    auto oo = S::Infinity();
    auto denom = [&](const Expr& A) {
        return pow(pow(x, integer(2)) + A, integer(-1));
    };

    REQUIRE(integrate(log(x) * denom(S::One()), x, zero, oo) == S::Zero());
    REQUIRE(oracle.equivalent(
        integrate(log(x) * denom(integer(4)), x, zero, oo)->str(),
        "pi*log(4)/8"));
    REQUIRE(oracle.equivalent(
        integrate(log(x) * denom(integer(9)), x, zero, oo)->str(),
        "pi*log(9)/12"));
    REQUIRE(oracle.equivalent(
        integrate(log(x) * denom(integer(2)), x, zero, oo)->str(),
        "pi*log(2)/(4*sqrt(2))"));
    // A constant multiplier carries through; log(1)·anything is still 0.
    REQUIRE(integrate(integer(3) * log(x) * denom(S::One()), x, zero, oo)
            == S::Zero());

    // Guards: no log factor (plain arctangent integral) and a finite bound.
    REQUIRE(oracle.equivalent(integrate(denom(S::One()), x, zero, oo)->str(),
                              "pi/2"));
    REQUIRE(integrate(log(x) * denom(S::One()), x, zero, S::One())
                ->str()
                .rfind("Integral(", 0) == 0);
}

// INT-CALLBUDGET-1: a non-elementary integrand whose partial-fraction × by-parts
// expansion branches without terminating — log(x)/(x²−1) splits into log(x)/(x±1)
// pieces that by-parts ping-pongs between log(x)/(x−a) and log(x−a)/x forever.
// The depth guard bounds linear recursion but not this exponential blow-up, so
// the call used to hang. A per-top-level call-count budget now bounds total work
// and bails to a clean unevaluated marker. Legitimate integrals (well under the
// budget) are unaffected — checked here alongside.
TEST_CASE("integrate: call-budget backstop on branching non-elementary integrands (INT-CALLBUDGET-1)",
          "[7][integrate][regression]") {
    auto x = symbol("x");

    // Used to hang; now returns a clean unevaluated Integral marker, terminating.
    auto r = integrate(log(x) * pow(pow(x, integer(2)) - integer(1),
                                    integer(-1)),
                       x);
    REQUIRE(r->str().rfind("Integral(", 0) == 0);

    // The backstop must not affect ordinary closed-form integrals.
    REQUIRE(integrate(pow(x, integer(3)) * exp(x), x)->str().find("Integral")
            == std::string::npos);
    REQUIRE(integrate(pow(pow(x, integer(4)) - integer(1), integer(-1)), x)
                ->str()
                .find("Integral") == std::string::npos);
    REQUIRE(integrate(pow(x, integer(2)) * atan(x), x)->str().find("Integral")
            == std::string::npos);
}

// INT-CSCH-1: ∫₀^∞ x^p/sinh(cx) = 2·Γ(p+1)·(1−2^{−(p+1)})·ζ(p+1)/c^{p+1} for
// p>0, c>0 (from 1/sinh = 2Σe^{−(2k+1)cx} and the odd-denominator zeta). Gives
// ∫₀^∞ x/sinh x = π²/4, x³/sinh x = π⁴/8, x²/sinh x = 7ζ(3)/2. Non-elementary
// antiderivative; verified numerically against SymPy quadrature.
TEST_CASE("integrate: ∫x^p/sinh(cx) over [0,oo) (INT-CSCH-1)",
          "[7][integrate][definite][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto zero = S::Zero();
    auto oo = S::Infinity();
    auto csch = [&](const Expr& c) {
        return pow(sinh(c * x), integer(-1));
    };

    REQUIRE(oracle.equivalent(integrate(x * csch(S::One()), x, zero, oo)->str(),
                              "pi**2/4"));
    REQUIRE(oracle.equivalent(
        integrate(pow(x, integer(3)) * csch(S::One()), x, zero, oo)->str(),
        "pi**4/8"));
    REQUIRE(oracle.equivalent(
        integrate(pow(x, integer(2)) * csch(S::One()), x, zero, oo)->str(),
        "7*zeta(3)/2"));
    // Scaled rate c = 2: π²/16.
    REQUIRE(oracle.equivalent(integrate(x * csch(integer(2)), x, zero, oo)->str(),
                              "pi**2/16"));
    // Constant multiplier carries through.
    REQUIRE(oracle.equivalent(
        integrate(integer(2) * x * csch(S::One()), x, zero, oo)->str(),
        "pi**2/2"));

    // Guards: cosh kernel (a different closed form) and the e^x−1 Bose kernel are
    // left to their own rules; the latter still evaluates.
    REQUIRE(integrate(x * pow(cosh(x), integer(-1)), x, zero, oo)
                ->str()
                .rfind("Integral(", 0) == 0);
    REQUIRE(oracle.equivalent(
        integrate(x * pow(exp(x) - integer(1), integer(-1)), x, zero, oo)->str(),
        "pi**2/6"));
}

// INT-LOG1PX2-1: ∫₀^∞ log(1+c·x²)/x² dx = π·√c for a positive constant c (by
// differentiating under the integral). The constant inside the log must be 1 so
// the integrand ~ c at 0 (convergent). Non-elementary antiderivative, so SymPP
// returned a marker. Gives ∫₀^∞ log(1+x²)/x² = π, log(1+4x²)/x² = 2π. Matches SymPy.
TEST_CASE("integrate: ∫log(1+cx²)/x² over [0,oo) (INT-LOG1PX2-1)",
          "[7][integrate][definite][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto zero = S::Zero();
    auto oo = S::Infinity();
    auto logterm = [&](const Expr& c) {
        return log(integer(1) + c * pow(x, integer(2)));
    };
    auto invx2 = pow(x, integer(-2));

    REQUIRE(oracle.equivalent(integrate(logterm(S::One()) * invx2, x, zero, oo)
                                  ->str(),
                              "pi"));
    REQUIRE(oracle.equivalent(
        integrate(logterm(integer(4)) * invx2, x, zero, oo)->str(), "2*pi"));
    REQUIRE(oracle.equivalent(
        integrate(logterm(integer(9)) * invx2, x, zero, oo)->str(), "3*pi"));
    REQUIRE(oracle.equivalent(
        integrate(logterm(integer(2)) * invx2, x, zero, oo)->str(),
        "sqrt(2)*pi"));
    // Leading constant carries through.
    REQUIRE(oracle.equivalent(
        integrate(integer(3) * logterm(S::One()) * invx2, x, zero, oo)->str(),
        "3*pi"));

    // Guard: a non-unit constant inside the log diverges at 0, so the rule must
    // not fire — left as an unevaluated marker.
    REQUIRE(integrate(log(integer(2) + pow(x, integer(2))) * invx2, x, zero, oo)
                ->str()
                .rfind("Integral(", 0) == 0);
}

// INT-INCGAMMA-1: ∫ xᵖ·e⁻ˣ dx = γ(p+1, x) for a symbolic exponent (the lower
// incomplete gamma is the antiderivative). Integer powers keep the elementary
// by-parts result; the definite [0,a] form folds to γ(s, a).
TEST_CASE("integrate: ∫x^(s-1)·e^{-x} = lowergamma(s, x) (INT-INCGAMMA-1)",
          "[7][integrate][incompletegamma][oracle][regression]") {
    auto x = symbol("x");
    auto s = symbol("s");
    auto a = symbol("a");
    auto emx = exp(mul(S::NegativeOne(), x));

    // Indefinite: symbolic exponent → lower incomplete gamma.
    REQUIRE(integrate(pow(x, add(s, integer(-1))) * emx, x)
            == lowergamma(s, x));
    REQUIRE(integrate(pow(x, s) * emx, x) == lowergamma(add(s, integer(1)), x));
    // Constant coefficient factors through.
    REQUIRE(integrate(integer(3) * pow(x, a) * emx, x)
            == mul(integer(3), lowergamma(add(a, integer(1)), x)));

    // The antiderivative differentiates back to the integrand.
    Expr F = integrate(pow(x, add(s, integer(-1))) * emx, x);
    REQUIRE(simplify(diff(F, x))
            == simplify(pow(x, add(s, integer(-1))) * emx));

    // Definite [0, a] folds to γ(s, a) (γ(s, 0) = 0).
    REQUIRE(integrate(pow(x, add(s, integer(-1))) * emx, x, S::Zero(), a)
            == lowergamma(s, a));

    // No regression: a non-negative integer power stays elementary (no gamma).
    Expr elem = integrate(pow(x, integer(2)) * emx, x);
    REQUIRE(elem->str().find("lowergamma") == std::string::npos);
    REQUIRE(simplify(diff(elem, x)) == simplify(pow(x, integer(2)) * emx));
}
