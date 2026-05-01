// integrate tests.

#include <catch2/catch_test_macros.hpp>

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
    // exp(x²) has no elementary antiderivative — and we don't have erf
    // recognition in the integral table.
    auto e = exp(pow(x, integer(2)));
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
