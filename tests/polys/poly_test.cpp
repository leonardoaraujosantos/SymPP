// Poly tests.
//
// References:
//   sympy/polys/tests/test_polytools.py — Poly.degree, .coeffs, .as_expr
//   roots: sympy/polys/polyroots.py

#include <algorithm>

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/expand.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/polys/poly.hpp>
#include <sympp/polys/rootof.hpp>
#include <sympp/core/float.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- Construction & accessors ----------------------------------------------

TEST_CASE("Poly: degree from expression", "[4][poly]") {
    auto x = symbol("x");
    auto e = pow(x, integer(3)) + integer(2) * pow(x, integer(2)) - integer(7) * x + integer(5);
    Poly p(e, x);
    REQUIRE(p.degree() == 3);
}

TEST_CASE("Poly: coefficients in increasing degree", "[4][poly]") {
    auto x = symbol("x");
    auto e = pow(x, integer(3)) + integer(2) * pow(x, integer(2))
             - integer(7) * x + integer(5);
    Poly p(e, x);
    REQUIRE(p.coeffs().size() == 4);
    REQUIRE(p.coeffs()[0] == integer(5));
    REQUIRE(p.coeffs()[1] == integer(-7));
    REQUIRE(p.coeffs()[2] == integer(2));
    REQUIRE(p.coeffs()[3] == integer(1));
}

TEST_CASE("Poly: leading_coeff and constant_term", "[4][poly]") {
    auto x = symbol("x");
    auto e = integer(3) * pow(x, integer(2)) + integer(7);
    Poly p(e, x);
    REQUIRE(p.leading_coeff() == integer(3));
    REQUIRE(p.constant_term() == integer(7));
}

// ----- Arithmetic ------------------------------------------------------------

TEST_CASE("Poly: addition", "[4][poly]") {
    auto x = symbol("x");
    Poly a(pow(x, integer(2)) + integer(1), x);
    Poly b(integer(2) * x, x);
    Poly c = a + b;
    REQUIRE(c.degree() == 2);
    REQUIRE(c.coeffs()[0] == integer(1));
    REQUIRE(c.coeffs()[1] == integer(2));
    REQUIRE(c.coeffs()[2] == integer(1));
}

TEST_CASE("Poly: multiplication", "[4][poly]") {
    auto x = symbol("x");
    // (x + 1) * (x - 1) = x^2 - 1
    Poly a(x + integer(1), x);
    Poly b(x - integer(1), x);
    Poly c = a * b;
    REQUIRE(c.degree() == 2);
    REQUIRE(c.coeffs()[0] == integer(-1));
    REQUIRE(c.coeffs()[1] == integer(0));
    REQUIRE(c.coeffs()[2] == integer(1));
}

// ----- eval ------------------------------------------------------------------

TEST_CASE("Poly: eval at integer", "[4][poly]") {
    auto x = symbol("x");
    Poly p(pow(x, integer(2)) + integer(2) * x + integer(3), x);
    // p(2) = 4 + 4 + 3 = 11
    REQUIRE(p.eval(integer(2)) == integer(11));
}

// ----- as_expr roundtrip -----------------------------------------------------

TEST_CASE("Poly: as_expr round-trips through SymPy", "[4][poly][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(x, integer(3)) - integer(7) * x + integer(2);
    Poly p(e, x);
    REQUIRE(oracle.equivalent(p.as_expr()->str(), "x**3 - 7*x + 2"));
}

// ----- Roots -----------------------------------------------------------------

TEST_CASE("Poly: linear root", "[4][poly][roots]") {
    auto x = symbol("x");
    Poly p(integer(3) * x - integer(7), x);
    auto roots = p.roots();
    REQUIRE(roots.size() == 1);
    REQUIRE(roots[0] == rational(7, 3));
}

TEST_CASE("Poly: quadratic roots", "[4][poly][roots]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // x^2 - 5x + 6 = (x - 2)(x - 3) -> roots 2, 3
    Poly p(pow(x, integer(2)) - integer(5) * x + integer(6), x);
    auto roots = p.roots();
    REQUIRE(roots.size() == 2);
    // The roots come in (-b ± sqrt(disc))/(2a) order; check both numeric values
    // collapse via the oracle to the expected set.
    REQUIRE(oracle.equivalent(roots[0]->str(), "3"));
    REQUIRE(oracle.equivalent(roots[1]->str(), "2"));
}

TEST_CASE("Poly: quadratic with irrational roots", "[4][poly][roots][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // x^2 - 2 -> roots ±sqrt(2)
    Poly p(pow(x, integer(2)) - integer(2), x);
    auto roots = p.roots();
    REQUIRE(roots.size() == 2);
    REQUIRE(oracle.equivalent(roots[0]->str(), "sqrt(2)"));
    REQUIRE(oracle.equivalent(roots[1]->str(), "-sqrt(2)"));
}

// ----- Polynomial division ---------------------------------------------------

namespace {
struct DivOracle { std::string quotient, remainder; };
[[nodiscard]] DivOracle oracle_div(std::string_view p, std::string_view q,
                                   std::string_view var) {
    auto& oracle = Oracle::instance();
    auto resp = oracle.send({{"op", "div"}, {"p", p}, {"q", q}, {"var", var}});
    REQUIRE(resp.ok);
    return DivOracle{resp.raw.at("quotient").get<std::string>(),
                     resp.raw.at("remainder").get<std::string>()};
}
}  // namespace

TEST_CASE("Poly: divmod exact (x^2 - 1) / (x - 1)", "[4][poly][divmod][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    Poly p(pow(x, integer(2)) - integer(1), x);
    Poly d(x - integer(1), x);
    auto [q, r] = p.divmod(d);
    auto ref = oracle_div("x**2 - 1", "x - 1", "x");
    REQUIRE(oracle.equivalent(q.as_expr()->str(), ref.quotient));
    REQUIRE(oracle.equivalent(r.as_expr()->str(), ref.remainder));
    REQUIRE(r.is_zero());
}

TEST_CASE("Poly: divmod with remainder (x^2 + 1) / (x + 1)", "[4][poly][divmod][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    Poly p(pow(x, integer(2)) + integer(1), x);
    Poly d(x + integer(1), x);
    auto [q, r] = p.divmod(d);
    auto ref = oracle_div("x**2 + 1", "x + 1", "x");
    REQUIRE(oracle.equivalent(q.as_expr()->str(), ref.quotient));
    REQUIRE(oracle.equivalent(r.as_expr()->str(), ref.remainder));
}

TEST_CASE("Poly: divmod cubic by quadratic", "[4][poly][divmod][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // (x^3 - 6x^2 + 11x - 6) / (x^2 - 3x + 2) = x - 3, remainder 0
    Poly p(pow(x, integer(3)) - integer(6) * pow(x, integer(2))
               + integer(11) * x - integer(6), x);
    Poly d(pow(x, integer(2)) - integer(3) * x + integer(2), x);
    auto [q, r] = p.divmod(d);
    auto ref = oracle_div("x**3 - 6*x**2 + 11*x - 6",
                          "x**2 - 3*x + 2", "x");
    REQUIRE(oracle.equivalent(q.as_expr()->str(), ref.quotient));
    REQUIRE(oracle.equivalent(r.as_expr()->str(), ref.remainder));
}

TEST_CASE("Poly: divmod with rational coefficients", "[4][poly][divmod][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // (2x^3 + 3x^2 - x + 5) / (2x + 1)
    Poly p(integer(2) * pow(x, integer(3)) + integer(3) * pow(x, integer(2))
               - x + integer(5), x);
    Poly d(integer(2) * x + integer(1), x);
    auto [q, r] = p.divmod(d);
    auto ref = oracle_div("2*x**3 + 3*x**2 - x + 5", "2*x + 1", "x");
    REQUIRE(oracle.equivalent(q.as_expr()->str(), ref.quotient));
    REQUIRE(oracle.equivalent(r.as_expr()->str(), ref.remainder));
}

TEST_CASE("Poly: divmod degree(divisor) > degree(dividend)", "[4][poly][divmod]") {
    auto x = symbol("x");
    // (x + 1) / (x^2 + 1)  -> q=0, r=x+1
    Poly p(x + integer(1), x);
    Poly d(pow(x, integer(2)) + integer(1), x);
    auto [q, r] = p.divmod(d);
    REQUIRE(q.is_zero());
    REQUIRE(r.as_expr() == p.as_expr());
}

TEST_CASE("Poly: divmod by zero throws", "[4][poly][divmod]") {
    auto x = symbol("x");
    Poly p(pow(x, integer(2)) + x, x);
    Poly z(std::vector<Expr>{S::Zero()}, x);
    REQUIRE_THROWS_AS(p.divmod(z), std::invalid_argument);
}

TEST_CASE("Poly: monic normalization", "[4][poly][divmod]") {
    auto x = symbol("x");
    Poly p(integer(3) * pow(x, integer(2)) + integer(6) * x + integer(9), x);
    auto m = p.monic();
    // Now leading coeff is 1, others divided by 3
    REQUIRE(m.leading_coeff() == S::One());
    REQUIRE(m.coeffs()[1] == integer(2));
    REQUIRE(m.coeffs()[0] == integer(3));
}

// ----- GCD -------------------------------------------------------------------

TEST_CASE("Poly: gcd of (x^2 - 1) and (x - 1)", "[4][poly][gcd][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    Poly a(pow(x, integer(2)) - integer(1), x);
    Poly b(x - integer(1), x);
    auto g = gcd(a, b);
    // Expect x - 1 (monic).
    auto resp = oracle.send({{"op", "gcd"},
                             {"a", "x**2 - 1"}, {"b", "x - 1"}});
    REQUIRE(resp.ok);
    REQUIRE(oracle.equivalent(g.as_expr()->str(),
                              resp.raw.at("result").get<std::string>()));
}

TEST_CASE("Poly: gcd is monic (rational scaling)", "[4][poly][gcd][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // (2x^2 - 2) and (3x - 3) — common factor (x - 1)
    Poly a(integer(2) * pow(x, integer(2)) - integer(2), x);
    Poly b(integer(3) * x - integer(3), x);
    auto g = gcd(a, b);
    auto resp = oracle.send({{"op", "gcd"},
                             {"a", "2*x**2 - 2"}, {"b", "3*x - 3"}});
    REQUIRE(resp.ok);
    REQUIRE(oracle.equivalent(g.as_expr()->str(),
                              resp.raw.at("result").get<std::string>()));
    // Sanity: leading coeff is 1.
    REQUIRE(g.leading_coeff() == S::One());
}

TEST_CASE("Poly: gcd of coprime polynomials is 1", "[4][poly][gcd][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    Poly a(pow(x, integer(2)) + integer(1), x);
    Poly b(x + integer(2), x);
    auto g = gcd(a, b);
    REQUIRE(oracle.equivalent(g.as_expr()->str(), "1"));
}

TEST_CASE("Poly: gcd shared cubic factor", "[4][poly][gcd][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // (x^3 - 1) = (x-1)(x^2+x+1) and (x^2 - 1) = (x-1)(x+1) — gcd is x - 1
    Poly a(pow(x, integer(3)) - integer(1), x);
    Poly b(pow(x, integer(2)) - integer(1), x);
    auto g = gcd(a, b);
    auto resp = oracle.send({{"op", "gcd"},
                             {"a", "x**3 - 1"}, {"b", "x**2 - 1"}});
    REQUIRE(resp.ok);
    REQUIRE(oracle.equivalent(g.as_expr()->str(),
                              resp.raw.at("result").get<std::string>()));
}

TEST_CASE("Poly: gcd with zero", "[4][poly][gcd]") {
    auto x = symbol("x");
    Poly a(integer(2) * pow(x, integer(2)) + integer(4), x);
    Poly z(std::vector<Expr>{S::Zero()}, x);
    auto g1 = gcd(a, z);
    auto g2 = gcd(z, a);
    // Both should equal the monic form of a: x^2 + 2
    REQUIRE(g1.leading_coeff() == S::One());
    REQUIRE(g2.leading_coeff() == S::One());
    REQUIRE(g1.coeffs()[0] == integer(2));
    REQUIRE(g2.coeffs()[0] == integer(2));
}

// ----- Square-free factorization (Yun) ---------------------------------------

TEST_CASE("Poly: sqf_list of (x+1)^2", "[4][poly][sqf][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // x^2 + 2x + 1 = (x + 1)^2
    Poly f(pow(x, integer(2)) + integer(2) * x + integer(1), x);
    auto sl = sqf_list(f);
    REQUIRE(sl.factors.size() == 1);
    REQUIRE(sl.factors[0].second == 2);
    REQUIRE(oracle.equivalent(sl.factors[0].first.as_expr()->str(), "x + 1"));
    REQUIRE(oracle.equivalent(sl.content->str(), "1"));
}

TEST_CASE("Poly: sqf_list strips content from 2(x+1)^2", "[4][poly][sqf][oracle]") {
    auto x = symbol("x");
    // 2x^2 + 4x + 2 = 2 (x+1)^2
    Poly f(integer(2) * pow(x, integer(2)) + integer(4) * x + integer(2), x);
    auto sl = sqf_list(f);
    REQUIRE(sl.content == integer(2));
    REQUIRE(sl.factors.size() == 1);
    REQUIRE(sl.factors[0].second == 2);
}

TEST_CASE("Poly: sqf_list of square-free polynomial", "[4][poly][sqf]") {
    auto x = symbol("x");
    // x^2 - 1 = (x-1)(x+1) — both multiplicity 1, so reported as one factor
    // x^2 - 1 (the whole thing since it's square-free).
    Poly f(pow(x, integer(2)) - integer(1), x);
    auto sl = sqf_list(f);
    REQUIRE(sl.factors.size() == 1);
    REQUIRE(sl.factors[0].second == 1);
    // The factor is the whole polynomial (already square-free).
    REQUIRE(sl.factors[0].first.degree() == 2);
}

TEST_CASE("Poly: sqf_list of (x-1)^3 (x+2)^2", "[4][poly][sqf][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // Build (x-1)^3 (x+2)^2 by expanding: ... or just construct via Poly mul.
    Poly a(x - integer(1), x);
    Poly b(x + integer(2), x);
    Poly f = a * a * a * b * b;
    auto sl = sqf_list(f);
    REQUIRE(sl.factors.size() == 2);
    // Sort by multiplicity to get a deterministic order.
    auto factors = sl.factors;
    std::sort(factors.begin(), factors.end(),
              [](auto& l, auto& r) { return l.second < r.second; });
    REQUIRE(factors[0].second == 2);
    REQUIRE(factors[1].second == 3);
    REQUIRE(oracle.equivalent(factors[0].first.as_expr()->str(), "x + 2"));
    REQUIRE(oracle.equivalent(factors[1].first.as_expr()->str(), "x - 1"));
}

TEST_CASE("Poly: sqf_list of zero", "[4][poly][sqf]") {
    auto x = symbol("x");
    Poly z(std::vector<Expr>{S::Zero()}, x);
    auto sl = sqf_list(z);
    REQUIRE(sl.factors.empty());
    REQUIRE(sl.content == S::Zero());
}

TEST_CASE("Poly: sqf_list of constant", "[4][poly][sqf]") {
    auto x = symbol("x");
    Poly c(std::vector<Expr>{integer(7)}, x);
    auto sl = sqf_list(c);
    REQUIRE(sl.factors.empty());
    REQUIRE(sl.content == integer(7));
}

// ----- Rational roots --------------------------------------------------------

TEST_CASE("Poly: rational_roots of x^2 - 4", "[4][poly][rational_roots]") {
    auto x = symbol("x");
    Poly f(pow(x, integer(2)) - integer(4), x);
    auto roots = rational_roots(f);
    REQUIRE(roots.size() == 2);
    // Order: synthetic division proceeds by smallest p first; ±2 in some order.
    std::sort(roots.begin(), roots.end(),
              [](const Expr& a, const Expr& b) { return a->str() < b->str(); });
    REQUIRE(roots[0] == integer(-2));
    REQUIRE(roots[1] == integer(2));
}

TEST_CASE("Poly: rational_roots with multiplicity", "[4][poly][rational_roots]") {
    auto x = symbol("x");
    // (x-3)^2 = x^2 - 6x + 9 — root 3 with multiplicity 2
    Poly f(pow(x, integer(2)) - integer(6) * x + integer(9), x);
    auto roots = rational_roots(f);
    REQUIRE(roots.size() == 2);
    REQUIRE(roots[0] == integer(3));
    REQUIRE(roots[1] == integer(3));
}

TEST_CASE("Poly: rational_roots p/q form", "[4][poly][rational_roots]") {
    auto x = symbol("x");
    // 6x^2 - 5x + 1 = (2x-1)(3x-1) — roots 1/2 and 1/3
    Poly f(integer(6) * pow(x, integer(2)) - integer(5) * x + integer(1), x);
    auto roots = rational_roots(f);
    REQUIRE(roots.size() == 2);
    std::sort(roots.begin(), roots.end(),
              [](const Expr& a, const Expr& b) { return a->str() < b->str(); });
    // Lex-sorted "1/2" < "1/3" (3 > 2 at the last char).
    REQUIRE(roots[0] == rational(1, 2));
    REQUIRE(roots[1] == rational(1, 3));
}

TEST_CASE("Poly: rational_roots zero root", "[4][poly][rational_roots]") {
    auto x = symbol("x");
    // x^3 - x — roots 0, 1, -1
    Poly f(pow(x, integer(3)) - x, x);
    auto roots = rational_roots(f);
    REQUIRE(roots.size() == 3);
}

TEST_CASE("Poly: rational_roots no rationals (irreducible quadratic)",
          "[4][poly][rational_roots]") {
    auto x = symbol("x");
    // x^2 + 1 has no rational roots
    Poly f(pow(x, integer(2)) + integer(1), x);
    auto roots = rational_roots(f);
    REQUIRE(roots.empty());
}

// ----- Cubic / quartic closed-form roots -------------------------------------

namespace {
[[nodiscard]] bool oracle_is_zero(std::string_view expr) {
    auto& oracle = Oracle::instance();
    if (oracle.equivalent(expr, "0")) return true;
    // Symbolic simplify can't always prove nested-radical Cardano outputs
    // are zero. Fall back to a high-precision numeric check.
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", expr}, {"prec", 50}, {"tol", 25}});
    return resp.ok && resp.raw.at("result").get<bool>();
}
}  // namespace

TEST_CASE("Poly: cubic with rational roots (deflation)", "[4][poly][roots][cubic][oracle]") {
    auto x = symbol("x");
    // (x-1)(x-2)(x-3) = x^3 - 6x^2 + 11x - 6
    Poly f(pow(x, integer(3)) - integer(6) * pow(x, integer(2))
               + integer(11) * x - integer(6), x);
    auto roots = f.roots();
    REQUIRE(roots.size() == 3);
    for (const auto& r : roots) {
        auto val = f.eval(r);
        REQUIRE(oracle_is_zero(val->str()));
    }
}

TEST_CASE("Poly: cubic x^3 - 8 (Cardano via deflation)", "[4][poly][roots][cubic][oracle]") {
    auto x = symbol("x");
    Poly f(pow(x, integer(3)) - integer(8), x);
    auto roots = f.roots();
    REQUIRE(roots.size() == 3);
    for (const auto& r : roots) {
        auto val = f.eval(r);
        REQUIRE(oracle_is_zero(val->str()));
    }
}

TEST_CASE("Poly: cubic x^3 + x + 1 (full Cardano, no rationals)",
          "[4][poly][roots][cubic][oracle]") {
    auto x = symbol("x");
    Poly f(pow(x, integer(3)) + x + integer(1), x);
    auto roots = f.roots();
    REQUIRE(roots.size() == 3);
    for (const auto& r : roots) {
        auto val = f.eval(r);
        REQUIRE(oracle_is_zero(val->str()));
    }
}

TEST_CASE("Poly: depressed cubic q = 0", "[4][poly][roots][cubic][oracle]") {
    auto x = symbol("x");
    // x^3 - 4x = x(x-2)(x+2). Has rational root 0; deflates.
    Poly f(pow(x, integer(3)) - integer(4) * x, x);
    auto roots = f.roots();
    REQUIRE(roots.size() == 3);
    for (const auto& r : roots) {
        auto val = f.eval(r);
        REQUIRE(oracle_is_zero(val->str()));
    }
}

TEST_CASE("Poly: quartic with all rational roots", "[4][poly][roots][quartic][oracle]") {
    auto x = symbol("x");
    // (x-1)(x+1)(x-2)(x+2) = (x^2-1)(x^2-4) = x^4 - 5x^2 + 4
    Poly f(pow(x, integer(4)) - integer(5) * pow(x, integer(2)) + integer(4), x);
    auto roots = f.roots();
    REQUIRE(roots.size() == 4);
    for (const auto& r : roots) {
        auto val = f.eval(r);
        REQUIRE(oracle_is_zero(val->str()));
    }
}

TEST_CASE("Poly: biquadratic x^4 + 5x^2 + 6", "[4][poly][roots][quartic][oracle]") {
    auto x = symbol("x");
    // No rational roots; biquadratic in x²: u² + 5u + 6 = (u+2)(u+3)
    // x² = -2 or x² = -3 — all four roots are complex.
    Poly f(pow(x, integer(4)) + integer(5) * pow(x, integer(2)) + integer(6), x);
    auto roots = f.roots();
    REQUIRE(roots.size() == 4);
    for (const auto& r : roots) {
        auto val = f.eval(r);
        REQUIRE(oracle_is_zero(val->str()));
    }
}

TEST_CASE("Poly: quintic returns rational roots only", "[4][poly][roots]") {
    auto x = symbol("x");
    // (x-1) * (x^4 + 1) — rational root 1; the rest is degree 4 quartic
    // with no rational roots (we'd recurse).
    // To keep this test focused, use (x-1)(x-2)(x^3 + 1) which has 5 roots
    // total, 2 rational (1, 2) and 3 from cubic via Cardano.
    // Build: (x-1)(x-2)(x³+1) — multiply manually:
    //   = (x²-3x+2)(x³+1) = x^5 + x² - 3x^4 - 3x + 2x³ + 2
    //   = x^5 - 3x^4 + 2x³ + x² - 3x + 2
    Poly f(pow(x, integer(5)) - integer(3) * pow(x, integer(4))
               + integer(2) * pow(x, integer(3)) + pow(x, integer(2))
               - integer(3) * x + integer(2), x);
    auto roots = f.roots();
    // Expect 5 roots: 2 rational + 3 cubic via Cardano (extracted recursively).
    REQUIRE(roots.size() == 5);
}

// ----- factor over Z (Kronecker) ---------------------------------------------

TEST_CASE("Poly: factor x^2 - 1", "[4][poly][factor][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(x, integer(2)) - integer(1);
    auto f = factor(e, x);
    REQUIRE(oracle.equivalent(f->str(), "(x - 1)*(x + 1)"));
}

TEST_CASE("Poly: factor x^2 + 2x + 1 = (x+1)^2", "[4][poly][factor][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(x, integer(2)) + integer(2) * x + integer(1);
    auto f = factor(e, x);
    REQUIRE(oracle.equivalent(f->str(), "(x + 1)**2"));
}

TEST_CASE("Poly: factor x^3 - 1", "[4][poly][factor][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(x, integer(3)) - integer(1);
    auto f = factor(e, x);
    // Over Z: (x - 1)(x^2 + x + 1)
    REQUIRE(oracle.equivalent(f->str(), "(x - 1)*(x**2 + x + 1)"));
}

TEST_CASE("Poly: factor x^4 - 1", "[4][poly][factor][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(x, integer(4)) - integer(1);
    auto f = factor(e, x);
    // (x-1)(x+1)(x²+1)
    REQUIRE(oracle.equivalent(f->str(), "(x - 1)*(x + 1)*(x**2 + 1)"));
}

TEST_CASE("Poly: factor 6x^2 - 5x + 1", "[4][poly][factor][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = integer(6) * pow(x, integer(2)) - integer(5) * x + integer(1);
    auto f = factor(e, x);
    REQUIRE(oracle.equivalent(f->str(), "(2*x - 1)*(3*x - 1)"));
}

TEST_CASE("Poly: factor irreducible quadratic stays put",
          "[4][poly][factor][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(x, integer(2)) + integer(1);
    auto f = factor(e, x);
    REQUIRE(oracle.equivalent(f->str(), "x**2 + 1"));
}

TEST_CASE("Poly: factor x^4 + 4 (Sophie Germain)",
          "[4][poly][factor][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // x^4 + 4 = (x²-2x+2)(x²+2x+2)
    auto e = pow(x, integer(4)) + integer(4);
    auto f = factor(e, x);
    REQUIRE(oracle.equivalent(
        f->str(), "(x**2 - 2*x + 2)*(x**2 + 2*x + 2)"));
}

TEST_CASE("Poly: factor with content 2(x-1)(x+1)",
          "[4][poly][factor][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // 2x^2 - 2 = 2(x-1)(x+1)
    auto e = integer(2) * pow(x, integer(2)) - integer(2);
    auto f = factor(e, x);
    REQUIRE(oracle.equivalent(f->str(), "2*(x - 1)*(x + 1)"));
}

TEST_CASE("Poly: factor degree-5 with mixed roots",
          "[4][poly][factor][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // (x-1)(x-2)(x³+1) = (x-1)(x-2)(x+1)(x²-x+1)
    // Build polynomial directly:
    auto e = pow(x, integer(5)) - integer(3) * pow(x, integer(4))
             + integer(2) * pow(x, integer(3)) + pow(x, integer(2))
             - integer(3) * x + integer(2);
    auto f = factor(e, x);
    REQUIRE(oracle.equivalent(
        f->str(),
        "(x - 2)*(x - 1)*(x + 1)*(x**2 - x + 1)"));
}

TEST_CASE("Poly: factor_list returns multiplicities",
          "[4][poly][factor]") {
    auto x = symbol("x");
    // (x-1)^3 (x+2)^2
    Poly a(x - integer(1), x);
    Poly b(x + integer(2), x);
    Poly f = a * a * a * b * b;
    auto fl = factor_list(f);
    REQUIRE(fl.factors.size() == 2);
    std::size_t total_mult = 0;
    for (auto& [p, m] : fl.factors) total_mult += m;
    REQUIRE(total_mult == 5);
}

// ----- together / cancel / apart / horner ------------------------------------

TEST_CASE("together: 1/x + 1/y", "[4][together][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto e = pow(x, integer(-1)) + pow(y, integer(-1));
    auto t = together(e);
    auto resp = oracle.send({{"op", "together"}, {"expr", "1/x + 1/y"}});
    REQUIRE(resp.ok);
    REQUIRE(oracle.equivalent(t->str(),
                              resp.raw.at("result").get<std::string>()));
}

TEST_CASE("together: 1/(x-1) - 1/(x+1)", "[4][together][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(x - integer(1), integer(-1)) - pow(x + integer(1), integer(-1));
    auto t = together(e);
    REQUIRE(oracle.equivalent(t->str(), "2/((x - 1)*(x + 1))"));
}

TEST_CASE("cancel: (x^2 - 1)/(x - 1) = x + 1", "[4][cancel][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = (pow(x, integer(2)) - integer(1)) / (x - integer(1));
    auto c = cancel(e, x);
    REQUIRE(oracle.equivalent(c->str(), "x + 1"));
}

TEST_CASE("cancel: (x^2 - 4)/(x^2 - 5x + 6) = (x+2)/(x-3)",
          "[4][cancel][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto num = pow(x, integer(2)) - integer(4);
    auto den = pow(x, integer(2)) - integer(5) * x + integer(6);
    auto e = num / den;
    auto c = cancel(e, x);
    REQUIRE(oracle.equivalent(c->str(), "(x + 2)/(x - 3)"));
}

TEST_CASE("cancel: nothing in common returns input", "[4][cancel][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = (x + integer(1)) / (x - integer(1));
    auto c = cancel(e, x);
    REQUIRE(oracle.equivalent(c->str(), "(x + 1)/(x - 1)"));
}

TEST_CASE("apart: (3x + 5)/((x-1)(x+2))", "[4][apart][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto num = integer(3) * x + integer(5);
    auto den = (x - integer(1)) * (x + integer(2));
    auto e = num / den;
    auto a = apart(e, x);
    auto resp = oracle.send({{"op", "apart"},
                             {"expr", "(3*x + 5)/((x-1)*(x+2))"},
                             {"var", "x"}});
    REQUIRE(resp.ok);
    REQUIRE(oracle.equivalent(a->str(),
                              resp.raw.at("result").get<std::string>()));
}

TEST_CASE("apart: 1/((x-1)(x-2)(x-3)) — three distinct linear",
          "[4][apart][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto den = (x - integer(1)) * (x - integer(2)) * (x - integer(3));
    auto e = pow(den, integer(-1));
    auto a = apart(e, x);
    auto resp = oracle.send({{"op", "apart"},
                             {"expr", "1/((x-1)*(x-2)*(x-3))"},
                             {"var", "x"}});
    REQUIRE(resp.ok);
    REQUIRE(oracle.equivalent(a->str(),
                              resp.raw.at("result").get<std::string>()));
}

TEST_CASE("apart: improper fraction extracts polynomial part",
          "[4][apart][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // (x^3 + 1)/(x^2 - 1) — improper. After division: x + (x+1)/(x²-1).
    // (x+1)/((x-1)(x+1)) cancels to 1/(x-1).
    // So result is x + 1/(x-1).
    auto e = (pow(x, integer(3)) + integer(1)) / (pow(x, integer(2)) - integer(1));
    auto a = apart(e, x);
    REQUIRE(oracle.equivalent(a->str(), "x + 1/(x - 1)"));
}

TEST_CASE("horner: 2x^3 + 3x^2 + x + 5 = x*(x*(2*x + 3) + 1) + 5",
          "[4][horner][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = integer(2) * pow(x, integer(3)) + integer(3) * pow(x, integer(2))
             + x + integer(5);
    auto h = horner(e, x);
    // Horner form is algebraically equivalent — let oracle confirm.
    REQUIRE(oracle.equivalent(h->str(),
                              "2*x**3 + 3*x**2 + x + 5"));
}

TEST_CASE("horner: degree 4 nests four deep", "[4][horner]") {
    auto x = symbol("x");
    auto e = pow(x, integer(4)) + integer(2) * pow(x, integer(3))
             + integer(3) * pow(x, integer(2)) + integer(4) * x + integer(5);
    auto h = horner(e, x);
    // Nested form should not contain x^k for k > 1 anywhere; structurally it's
    // just adds and muls of x. Verify via expansion equality.
    REQUIRE(expand(h) == expand(e));
}

// ----- Resultant / Discriminant ----------------------------------------------

TEST_CASE("resultant: x^2 - 1, x - 1 share root → 0",
          "[4][resultant][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = resultant(pow(x, integer(2)) - integer(1), x - integer(1), x);
    REQUIRE(oracle.equivalent(e->str(), "0"));
}

TEST_CASE("resultant: coprime polys → nonzero",
          "[4][resultant][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // x² + 1 and x + 2 — no common root over ℚ. Resultant = (-2)² + 1 = 5.
    auto e = resultant(pow(x, integer(2)) + integer(1), x + integer(2), x);
    auto resp = oracle.send({{"op", "resultant"},
                             {"p", "x**2 + 1"}, {"q", "x + 2"},
                             {"var", "x"}});
    REQUIRE(resp.ok);
    REQUIRE(oracle.equivalent(e->str(),
                              resp.raw.at("result").get<std::string>()));
}

TEST_CASE("resultant: arbitrary cubic vs quadratic",
          "[4][resultant][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto p = pow(x, integer(3)) - integer(2) * x + integer(1);
    auto q = pow(x, integer(2)) + integer(3) * x - integer(4);
    auto e = resultant(p, q, x);
    auto resp = oracle.send({{"op", "resultant"},
                             {"p", "x**3 - 2*x + 1"},
                             {"q", "x**2 + 3*x - 4"},
                             {"var", "x"}});
    REQUIRE(resp.ok);
    REQUIRE(oracle.equivalent(e->str(),
                              resp.raw.at("result").get<std::string>()));
}

TEST_CASE("discriminant: quadratic ax^2 + bx + c → b² - 4ac",
          "[4][discriminant][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // 2x² + 5x + 3 — disc = 25 - 24 = 1
    auto e = discriminant(integer(2) * pow(x, integer(2)) + integer(5) * x
                          + integer(3), x);
    REQUIRE(oracle.equivalent(e->str(), "1"));
}

TEST_CASE("discriminant: x^3 - 1 has discriminant -27",
          "[4][discriminant][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = discriminant(pow(x, integer(3)) - integer(1), x);
    REQUIRE(oracle.equivalent(e->str(), "-27"));
}

TEST_CASE("discriminant: cubic with double root vanishes",
          "[4][discriminant][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // (x-1)²(x-2) = x³ - 4x² + 5x - 2 — has a double root, disc = 0.
    auto e = discriminant(pow(x, integer(3)) - integer(4) * pow(x, integer(2))
                          + integer(5) * x - integer(2), x);
    REQUIRE(oracle.equivalent(e->str(), "0"));
}

// ----- RootOf ----------------------------------------------------------------

TEST_CASE("RootOf: construction and printing", "[4][rootof]") {
    auto x = symbol("x");
    auto p = pow(x, integer(5)) - x - integer(1);
    auto r = root_of(p, x, 0);
    auto s = r->str();
    REQUIRE(s.starts_with("CRootOf("));
    REQUIRE(s.ends_with(", 0)"));
}

TEST_CASE("RootOf: equality and hash-cons", "[4][rootof]") {
    auto x = symbol("x");
    auto p = pow(x, integer(3)) + x + integer(1);
    auto r1 = root_of(p, x, 0);
    auto r2 = root_of(p, x, 0);
    REQUIRE(r1 == r2);
    REQUIRE(r1.get() == r2.get());  // hash-cons identity
}

TEST_CASE("RootOf: distinct indices not equal", "[4][rootof]") {
    auto x = symbol("x");
    auto p = pow(x, integer(3)) + x + integer(1);
    auto r0 = root_of(p, x, 0);
    auto r1 = root_of(p, x, 1);
    REQUIRE(r0 != r1);
}

TEST_CASE("RootOf: evalf of x^2 - 2 root 0 ≈ -sqrt(2)",
          "[4][rootof][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto p = pow(x, integer(2)) - integer(2);
    auto r = root_of(p, x, 0);
    auto v = evalf(r, 25);
    // Real roots in ascending order: -sqrt(2), +sqrt(2). Index 0 = -sqrt(2).
    auto resp = oracle.send({{"op", "evalf"}, {"expr", "-sqrt(2)"}, {"prec", 25}});
    REQUIRE(resp.ok);
    // Compare evalf strings — both should agree to ~25 digits.
    REQUIRE(oracle.equivalent(v->str() + " - ("
                              + resp.raw.at("result").get<std::string>() + ")",
                              "0"));
}

TEST_CASE("RootOf: evalf x^3 - 2 (real cube root of 2)",
          "[4][rootof][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto p = pow(x, integer(3)) - integer(2);
    auto r = root_of(p, x, 0);  // unique real root: 2^(1/3)
    auto v = evalf(r, 20);
    // Numerical compare to high precision: should be ≈ 1.2599210498948732...
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", v->str() + " - 2**(1/3)"},
                             {"prec", 30}, {"tol", 18}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

TEST_CASE("RootOf: evalf irreducible quintic x^5 - x - 1",
          "[4][rootof][oracle]") {
    auto x = symbol("x");
    // Bring factor to attention: famous irreducible quintic with one real
    // root ~1.1673... — Abel-Ruffini case; closed-form impossible.
    auto p = pow(x, integer(5)) - x - integer(1);
    auto r = root_of(p, x, 0);
    auto v = evalf(r, 15);
    // Substitute back: poly evaluated at v should be near zero.
    // We can't easily plug v into Poly::eval here without numeric eval of
    // the resulting Add — but we can compare to the known value 1.1673...
    auto& oracle = Oracle::instance();
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", v->str() + " - 1.1673039782614187"},
                             {"prec", 30}, {"tol", 10}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

TEST_CASE("RootOf: evalf out-of-range index returns self",
          "[4][rootof]") {
    auto x = symbol("x");
    // x^2 + 1 — no real roots. RootOf(_, _, 0) should not evalf.
    auto p = pow(x, integer(2)) + integer(1);
    auto r = root_of(p, x, 0);
    auto v = evalf(r, 15);
    // Should return RootOf unchanged (no real roots).
    REQUIRE(v == r);
}

TEST_CASE("Poly: diff lowers degree", "[4][poly][diff]") {
    auto x = symbol("x");
    // (x^3 + 2x^2 - 7x + 5)' = 3x^2 + 4x - 7
    Poly p(pow(x, integer(3)) + integer(2) * pow(x, integer(2))
               - integer(7) * x + integer(5), x);
    auto dp = p.diff();
    REQUIRE(dp.degree() == 2);
    REQUIRE(dp.coeffs()[2] == integer(3));
    REQUIRE(dp.coeffs()[1] == integer(4));
    REQUIRE(dp.coeffs()[0] == integer(-7));
}
