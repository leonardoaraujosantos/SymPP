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
#include <sympp/core/type_id.hpp>
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

// POLY-FACTOR-ROOTS-1: a quartic with no rational root that nevertheless
// factors over ℚ into two quadratics (x⁴+x²+1 = (x²+x+1)(x²−x+1)) must be
// solved through that factorization — yielding clean ±1/2 ± √3·i/2 roots —
// rather than via Ferrari's resolvent, which returns nested radicals like
// sqrt((I*sqrt(3) - 1)/2). The factoring path also makes higher-degree
// products solvable.
TEST_CASE("Poly: factors reducible quartic into clean roots (POLY-FACTOR-ROOTS-1)",
          "[4][poly][roots][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto set_equal = [&](const std::vector<Expr>& got,
                         const std::vector<std::string>& want) {
        if (got.size() != want.size()) return false;
        std::vector<bool> used(got.size(), false);
        for (const auto& w : want) {
            bool hit = false;
            for (std::size_t i = 0; i < got.size(); ++i) {
                if (!used[i] && oracle.equivalent(got[i]->str(), w)) {
                    used[i] = hit = true;
                    break;
                }
            }
            if (!hit) return false;
        }
        return true;
    };
    // x⁴ + x² + 1 = (x²+x+1)(x²−x+1): the four primitive 6th/3rd roots.
    Poly q(pow(x, integer(4)) + pow(x, integer(2)) + integer(1), x);
    auto qr = q.roots();
    REQUIRE(set_equal(qr, {"-1/2 + sqrt(3)*I/2", "-1/2 - sqrt(3)*I/2",
                           "1/2 + sqrt(3)*I/2", "1/2 - sqrt(3)*I/2"}));
    // None of the roots should be a nested radical (no "**(1/2)" wrapping a
    // complex subexpression — i.e. no Ferrari fallback).
    for (const auto& r : qr) {
        REQUIRE(r->str().find(")**(1/2)") == std::string::npos);
    }
    // x⁶ − 1: ±1 plus the four roots above.
    Poly s(pow(x, integer(6)) - integer(1), x);
    REQUIRE(set_equal(s.roots(),
                      {"1", "-1", "-1/2 + sqrt(3)*I/2", "-1/2 - sqrt(3)*I/2",
                       "1/2 + sqrt(3)*I/2", "1/2 - sqrt(3)*I/2"}));
    // Degree-5 reducible with no rational root: (x²+x+1)(x³+2) → 5 roots.
    Poly d5(expand((pow(x, integer(2)) + x + integer(1))
                   * (pow(x, integer(3)) + integer(2))),
            x);
    REQUIRE(d5.roots().size() == 5);
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

// Regression (CANCEL-1): cancel()/factor() used to hang on polynomials with
// symbolic (multivariate) coefficients — the divmod GCD loop never terminated
// because (b+b²) − (b+b²) did not fold to a structural zero. divmod now expands
// the coefficient subtraction (and has a degree-decrease backstop), and factor
// bails to the input when coefficients are non-numeric.
TEST_CASE("Poly: cancel/factor terminate on symbolic coefficients (CANCEL-1)",
          "[4][poly][cancel][factor][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto a = symbol("a");
    auto b = symbol("b");
    // Multivariate cancel now reduces instead of hanging.
    REQUIRE(oracle.equivalent(
        cancel((pow(x, integer(2)) - pow(y, integer(2))) / (x - y), x)->str(),
        "x + y"));
    REQUIRE(oracle.equivalent(
        cancel((pow(a, integer(2)) - pow(b, integer(2))) / (a + b), a)->str(),
        "a - b"));
    // The documented hang case terminates and stays equivalent.
    auto c = cancel((b - a + integer(1)) * (a + b) / integer(2), a);
    REQUIRE(oracle.equivalent(c->str(), "(b - a + 1)*(a + b)/2"));
    // factor over a symbolic coefficient no longer hangs (returns unfactored).
    auto f = factor(pow(x, integer(2)) - pow(y, integer(2)), x);
    REQUIRE(oracle.equivalent(f->str(), "x**2 - y**2"));
}

TEST_CASE("Poly: factor x^2 + 2x + 1 = (x+1)^2", "[4][poly][factor][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto e = pow(x, integer(2)) + integer(2) * x + integer(1);
    auto f = factor(e, x);
    REQUIRE(oracle.equivalent(f->str(), "(x + 1)**2"));
}

// FACTOR-NONMONIC-POW-1: a perfect power of a non-monic linear factor must not
// leak its leading coefficient into the content. 4x²+4x+1 = (2x+1)² (content 1),
// previously mis-factored as 2·(2x+1)² (= 8x²+8x+2, numerically wrong). The
// primitive-part scalar of a multiplicity-m factor enters the content as qᵐ.
TEST_CASE("Poly: factor non-monic perfect powers (FACTOR-NONMONIC-POW-1)",
          "[4][poly][factor][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto sq = [&](const Expr& e) { return pow(e, integer(2)); };
    // (2x+1)² = 4x²+4x+1 — content must be 1, not 2.
    REQUIRE(oracle.equivalent(
        factor(integer(4) * sq(x) + integer(4) * x + integer(1), x)->str(),
        "(2*x + 1)**2"));
    // (3x+1)² and a cube (2x+1)³ = 8x³+12x²+6x+1.
    REQUIRE(oracle.equivalent(
        factor(integer(9) * sq(x) + integer(6) * x + integer(1), x)->str(),
        "(3*x + 1)**2"));
    REQUIRE(oracle.equivalent(
        factor(integer(8) * pow(x, integer(3)) + integer(12) * sq(x)
                   + integer(6) * x + integer(1),
               x)->str(),
        "(2*x + 1)**3"));
    // Genuine content is still pulled out: 2x²+4x+2 = 2·(x+1)².
    REQUIRE(oracle.equivalent(
        factor(integer(2) * sq(x) + integer(4) * x + integer(2), x)->str(),
        "2*(x + 1)**2"));
    // Round-trip: factor(p) expands back to p exactly.
    auto p = integer(16) * pow(x, integer(4)) - integer(8) * sq(x) + integer(1);
    REQUIRE(expand(factor(p, x)) == expand(p));
    // Distinct non-monic linears (multiplicity 1) unaffected: 6x²+11x+3.
    REQUIRE(oracle.equivalent(
        factor(integer(6) * sq(x) + integer(11) * x + integer(3), x)->str(),
        "(2*x + 3)*(3*x + 1)"));
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

// FACTOR-HOM-1: homogeneous bivariate polynomials factor via dehomogenize →
// factor-over-ℚ → re-homogenize, verified by re-expansion. Covers the common
// difference-of-squares / cubes / perfect-square-trinomial multivariate cases.
TEST_CASE("Poly: factor homogeneous bivariate (FACTOR-HOM-1)",
          "[4][poly][factor][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    REQUIRE(oracle.equivalent(
        factor(pow(x, integer(2)) - pow(y, integer(2)), x)->str(),
        "(x - y)*(x + y)"));
    REQUIRE(oracle.equivalent(
        factor(pow(x, integer(2)) + integer(2) * x * y + pow(y, integer(2)), x)
            ->str(),
        "(x + y)**2"));
    REQUIRE(oracle.equivalent(
        factor(pow(x, integer(3)) - pow(y, integer(3)), x)->str(),
        "(x - y)*(x**2 + x*y + y**2)"));
    // A pure-other factor is pulled out: x²y − y³ = y(x−y)(x+y).
    REQUIRE(oracle.equivalent(
        factor(pow(x, integer(2)) * y - pow(y, integer(3)), x)->str(),
        "y*(x - y)*(x + y)"));
    // Irreducible over ℚ stays put (self-verification rejects a bad split).
    REQUIRE(oracle.equivalent(
        factor(pow(x, integer(2)) + pow(y, integer(2)), x)->str(),
        "x**2 + y**2"));
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

// TOGETHER-LCM-1: combine a sum of fractions over the LCM of the denominators,
// not their product — so a shared denominator factor isn't duplicated
// (a/b + c/b → (a+c)/b, not (a·b + b·c)/b²). Fixes as_numer_denom, which feeds
// together / cancel / apart / simplify.
TEST_CASE("together: shared denominator uses the LCM (TOGETHER-LCM-1)",
          "[4][together][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto a = symbol("a");
    auto b = symbol("b");
    auto c = symbol("c");
    // a/b + c/b → (a + c)/b.
    REQUIRE(oracle.equivalent(
        together(a / b + c / b)->str(), "(a + c)/b"));
    // 3 terms over the same denominator.
    auto d = symbol("d");
    REQUIRE(oracle.equivalent(
        together(a / b + c / b + d / b)->str(), "(a + c + d)/b"));
    // Shared (x+1) collapses to 1 instead of an (x+1)² bloat.
    REQUIRE(together(x / (x + integer(1)) + integer(1) / (x + integer(1)))
            == integer(1));
    // Power relationship: LCM is (x-1)², not (x-1)³.
    REQUIRE(oracle.equivalent(
        together(pow(x - integer(1), integer(-1))
                 + pow(x - integer(1), integer(-2)))->str(),
        "x/(x - 1)**2"));
    // Distinct denominators still cross-multiply correctly.
    auto y = symbol("y");
    REQUIRE(oracle.equivalent(
        together(pow(x, integer(-1)) + pow(y, integer(-1)))->str(),
        "(x + y)/(x*y)"));
}

// TOGETHER-NESTED-1: compound (nested) fractions collapse bottom-up. together
// now recurses into Add/Mul/Pow children, so a reciprocal of a sum of fractions
// reduces — 1/(1 + 1/x) → x/(x + 1) — and nests arbitrarily deep. Matches SymPy.
TEST_CASE("together: nested compound fractions (TOGETHER-NESTED-1)",
          "[4][together][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto one = integer(1);
    auto recip = [](const Expr& e) { return pow(e, integer(-1)); };
    // 1/(1 + 1/x) → x/(x + 1) ; 1/(1 − 1/x) → x/(x − 1).
    REQUIRE(oracle.equivalent(
        together(recip(one + recip(x)))->str(), "x/(x + 1)"));
    REQUIRE(oracle.equivalent(
        together(recip(one - recip(x)))->str(), "x/(x - 1)"));
    // Twice-nested: 1/(1 + 1/(1 + 1/x)) → (x + 1)/(2x + 1).
    REQUIRE(oracle.equivalent(
        together(recip(one + recip(one + recip(x))))->str(),
        "(x + 1)/(2*x + 1)"));
    // Mixed symbols: a/(b + c/d) → a*d/(b*d + c).
    auto a = symbol("a"), b = symbol("b"), c = symbol("c"), d = symbol("d");
    REQUIRE(oracle.equivalent(
        together(a * recip(b + c * recip(d)))->str(), "a*d/(b*d + c)"));
    // No regression: a plain sum of reciprocals is unaffected.
    auto y = symbol("y");
    REQUIRE(oracle.equivalent(
        together(recip(x) + recip(y))->str(), "(x + y)/(x*y)"));
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

// POLYOP-1: parser-facing univariate polynomial operations with the variable
// inferred from the single free symbol — degree, quo, rem, and 1-argument
// cancel. Previously these parsed to opaque unevaluated function nodes.
TEST_CASE("poly ops: degree/quo/rem/cancel infer the variable (POLYOP-1)",
          "[4][poly][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto sq = [&](const Expr& e) { return pow(e, integer(2)); };
    // degree.
    REQUIRE(degree(pow(x, integer(3)) + integer(2) * x) == integer(3));
    REQUIRE(degree(integer(5)) == integer(0));
    // quo / rem: x²−1 = (x+1)(x−1) + 0; x² = (x+1)(x−1) + 1.
    REQUIRE(oracle.equivalent(quo(sq(x) - integer(1), x - integer(1))->str(),
                              "x + 1"));
    REQUIRE(rem(sq(x), x - integer(1)) == integer(1));
    REQUIRE(oracle.equivalent(
        quo(pow(x, integer(3)) - integer(1), x - integer(1))->str(),
        "x**2 + x + 1"));
    // 1-argument cancel infers the variable.
    REQUIRE(oracle.equivalent(
        cancel((sq(x) - integer(1)) / (x - integer(1)))->str(), "x + 1"));
}

// POLYOP-2: parser-facing resultant (two-arg) and discriminant (one-arg) with
// the variable inferred from the free symbols. Previously opaque function nodes.
TEST_CASE("poly ops: resultant/discriminant infer the variable (POLYOP-2)",
          "[4][poly][regression]") {
    auto x = symbol("x");
    auto sq = [&](const Expr& e) { return pow(e, integer(2)); };
    // discriminant of a quadratic: b² − 4ac.
    REQUIRE(discriminant(sq(x) + integer(2) * x + integer(1)) == integer(0));
    REQUIRE(discriminant(sq(x) - integer(5) * x + integer(6)) == integer(1));
    REQUIRE(discriminant(sq(x) + integer(1)) == integer(-4));
    REQUIRE(discriminant(sq(x) - integer(2)) == integer(8));
    // resultant: zero iff the polynomials share a root.
    REQUIRE(resultant(sq(x) - integer(1), x - integer(1)) == integer(0));
    REQUIRE(resultant(sq(x) + integer(1), x - integer(2)) == integer(5));
    // Sign convention matches SymPy: res(x−1, x−2) = −1, res(x−2, x−1) = 1.
    REQUIRE(resultant(x - integer(1), x - integer(2)) == integer(-1));
    REQUIRE(resultant(x - integer(2), x - integer(1)) == integer(1));
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

// Regression (APART-1): partial fractions over irreducible quadratic factors.
// apart() previously bailed (returned the input) when the denominator had a
// non-linear irreducible factor; it now decomposes the squarefree case.
TEST_CASE("apart: decomposes over irreducible quadratic factors",
          "[4][apart][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // 1/(x³+1) = 1/(3(x+1)) − (x−2)/(3(x²−x+1)).
    auto a1 = apart(pow(pow(x, integer(3)) + integer(1), integer(-1)), x);
    REQUIRE(oracle.equivalent(
        a1->str(), "1/(3*(x+1)) - (x-2)/(3*(x**2-x+1))"));
    // The decomposition recombines to the original fraction.
    REQUIRE(oracle.equivalent(a1->str(), "1/(x**3+1)"));
    // 1/(x⁴−1) mixes linear and irreducible-quadratic factors.
    auto a2 = apart(pow(pow(x, integer(4)) - integer(1), integer(-1)), x);
    REQUIRE(oracle.equivalent(a2->str(), "1/(x**4-1)"));
    REQUIRE(a2->str().find("x**2 + 1") != std::string::npos);  // has the quad term
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

// WANG-1: bivariate (multivariate) factorization — content extraction plus the
// quadratic-in-var case via a perfect-square discriminant, with denominator
// clearing and numeric-content stripping. Each result is verified to expand back
// to the input and cross-checked against SymPy.
TEST_CASE("factor: bivariate multivariate factorization (WANG-1)",
          "[4][polys][factor][wang][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x"), y = symbol("y");
    auto P = [&](const Expr& b, int e) { return pow(b, integer(e)); };
    auto check = [&](const Expr& e, const char* expected) {
        Expr f = factor(e, x);
        REQUIRE(oracle.equivalent(f->str(), e->str()));         // expands back
        REQUIRE(oracle.equivalent(f->str(), expected));         // matches SymPy form
        REQUIRE(f->str() != e->str());                          // actually factored
    };
    // Polynomial-in-y coefficients (true multivariate), previously unfactored:
    check(mul(P(x, 2), P(y, 2)) - integer(1), "(x*y - 1)*(x*y + 1)");   // x²y²−1
    check(P(x, 2) + mul(x, y) + x + y, "(x + 1)*(x + y)");              // x²+xy+x+y
    check(P(x + integer(1), 2) - P(y, 2), "(x - y + 1)*(x + y + 1)");   // (x+1)²−y²
    check(add(mul(P(x, 2), y), add(mul(integer(2), mul(x, y)), y)),
          "y*(x + 1)**2");                                              // x²y+2xy+y
    check(mul(x, y) + x + y + integer(1), "(x + 1)*(y + 1)");           // content (y+1)

    // A non-factorable bivariate quadratic stays put (no false factorization).
    Expr irr = P(x, 2) + mul(x, y) + integer(1);  // x²+xy+1 irreducible over ℚ[y]
    REQUIRE(factor(irr, x)->str() == irr->str());
}

// WANG-2: higher-degree (cubic and beyond) bivariate factorization — find a
// closed-form monomial root r(y) of the primitive part, deflate by (x − r),
// and recurse on the lower-degree quotient (the quadratic disc path or a
// further deflation). Every result is verified to expand back to the input and
// cross-checked against SymPy; a polynomial with no such root stays put.
TEST_CASE("factor: higher-degree bivariate (WANG-2)",
          "[4][polys][factor][wang][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x"), y = symbol("y");
    auto P = [&](const Expr& b, int e) { return pow(b, integer(e)); };
    auto check = [&](const Expr& e, const char* expected) {
        Expr f = factor(e, x);
        REQUIRE(oracle.equivalent(f->str(), e->str()));   // expands back
        REQUIRE(oracle.equivalent(f->str(), expected));   // matches SymPy form
        REQUIRE(f->str() != e->str());                    // actually factored
    };

    // x³ − y³ → (x − y)(x² + x·y + y²) : root r = y, deflated quotient quadratic.
    check(P(x, 3) - P(y, 3), "(x - y)*(x**2 + x*y + y**2)");

    // x³ + y·x² − x − y → (x − 1)(x + 1)(x + y) : root r = −y, then the
    // quotient quadratic x² − 1 splits via the disc path.
    check(P(x, 3) + mul(y, P(x, 2)) - x - y, "(x - 1)*(x + 1)*(x + y)");

    // x³ − y²·x → x·(x − y)(x + y) : content x factored out, cubic primitive
    // part x² − y² deflated.
    check(P(x, 3) - mul(P(y, 2), x), "x*(x - y)*(x + y)");

    // x³ − y·x² + x − y → (x − y)(x² + 1) : root r = y, quotient quadratic
    // x² + 1 stays irreducible over ℚ (still a valid more-factored result).
    check(P(x, 3) - mul(y, P(x, 2)) + x - y, "(x - y)*(x**2 + 1)");

    // x⁴ − y⁴ → (x − y)(x + y)(x² + y²) : two successive deflations.
    check(P(x, 4) - P(y, 4), "(x - y)*(x + y)*(x**2 + y**2)");

    // A cubic with no rational/monomial root in y stays put (no false split).
    Expr irr = P(x, 3) + x + y;  // x³ + x + y, irreducible over ℚ[y]
    REQUIRE(factor(irr, x)->str() == irr->str());
}

// WANG-3: trivariate symmetric cubic — the root of x³+y³+z³−3xyz in x is the
// linear form x = −(y+z), not a monomial. factor() treats the input as a
// univariate polynomial in x with coefficients in ℚ[y,z], searches candidate
// roots that are small linear forms over {y,z}, deflates by (x − r), and
// recurses on the irreducible quadratic quotient. Cross-checked against SymPy.
TEST_CASE("factor: trivariate symmetric cubic (WANG-3)",
          "[4][polys][factor][wang][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x"), y = symbol("y"), z = symbol("z");
    auto P = [&](const Expr& b, int e) { return pow(b, integer(e)); };

    // x³ + y³ + z³ − 3xyz → (x + y + z)(x² − xy − xz + y² − yz + z²)
    Expr e = add({P(x, 3), P(y, 3), P(z, 3),
                  mul(integer(-3), mul(x, mul(y, z)))});
    Expr f = factor(e, x);
    REQUIRE(oracle.equivalent(f->str(), e->str()));  // expands back to input
    REQUIRE(oracle.equivalent(
        f->str(),
        "(x + y + z)*(x**2 - x*y - x*z + y**2 - y*z + z**2)"));  // SymPy form
    REQUIRE(f->str() != e->str());                               // actually factored

    // A trivariate cubic with no small linear-form root in x stays put.
    Expr irr = add({P(x, 3), x, y, z});  // x³ + x + y + z
    REQUIRE(factor(irr, x)->str() == irr->str());
}

// WANG-4: degree-2 multivariate (>= 2 other free variables). The bivariate-only
// quadratic-discriminant branch does not apply, so the linear-form deflation
// path handles it — find a small linear-form root r(y, z), deflate by (x − r),
// and recurse on the linear quotient. (x + y)² − z² has root x = −y + z (or
// −y − z), both inside find_linear_form_root's {−1,0,1} coefficient grid.
TEST_CASE("factor: degree-2 multivariate (WANG-4)",
          "[4][polys][factor][wang][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x"), y = symbol("y"), z = symbol("z");
    auto P = [&](const Expr& b, int e) { return pow(b, integer(e)); };
    auto check = [&](const Expr& e, const char* expected) {
        Expr f = factor(e, x);
        REQUIRE(oracle.equivalent(f->str(), e->str()));   // expands back
        REQUIRE(oracle.equivalent(f->str(), expected));   // matches SymPy form
        REQUIRE(f->str() != e->str());                    // actually factored
    };

    // (x + y)² − z² → (x + y − z)(x + y + z)
    check(P(add({x, y}), 2) - P(z, 2), "(x + y - z)*(x + y + z)");
    // Expanded form: x² + 2xy + y² − z² → same factorization.
    check(add({P(x, 2), mul(integer(2), mul(x, y)), P(y, 2),
               mul(integer(-1), P(z, 2))}),
          "(x + y - z)*(x + y + z)");

    // A degree-2 trivariate quadratic with no small linear-form root stays put.
    Expr irr = add({P(x, 2), mul(x, y), z, integer(1)});  // x² + xy + z + 1
    REQUIRE(factor(irr, x)->str() == irr->str());
}

// WANG-5: 4-variable factorization. With three other free variables the
// linear-form / quadratic-disc bivariate paths do not apply, so factor() uses
// the general multivariate handler:
//   (1) multivariate CONTENT extraction — when the var-coefficients share a
//       polynomial factor in the other variables, pull it out and factor it
//       recursively (handles the linear/bilinear cases below);
//   (2) generalized difference of squares — A·var² + C with A and −C perfect
//       squares as polynomials in the other vars (handles x²y² − z²w²).
// Each result is verified to expand back to the input and cross-checked against
// SymPy; a polynomial with no such structure stays put.
TEST_CASE("factor: 4-variable (WANG-5)",
          "[4][polys][factor][wang][oracle]") {
    auto& oracle = Oracle::instance();
    auto a = symbol("a"), b = symbol("b"), c = symbol("c"), d = symbol("d");
    auto w = symbol("w"), x = symbol("x"), y = symbol("y"), z = symbol("z");
    auto P = [&](const Expr& bb, int e) { return pow(bb, integer(e)); };
    auto check = [&](const Expr& e, const Expr& var, const char* expected) {
        Expr f = factor(e, var);
        REQUIRE(oracle.equivalent(f->str(), e->str()));   // expands back
        REQUIRE(oracle.equivalent(f->str(), expected));   // matches SymPy form
        REQUIRE(f->str() != e->str());                    // actually factored
    };

    // a·b − a·c + b·d − c·d → (a + d)(b − c) : linear in a, content (b − c).
    check(add({mul(a, b), mul(integer(-1), mul(a, c)),
               mul(b, d), mul(integer(-1), mul(c, d))}),
          a, "(a + d)*(b - c)");

    // w·x + w·y + z·x + z·y → (w + z)(x + y) : linear in w, content (x + y).
    check(add({mul(w, x), mul(w, y), mul(z, x), mul(z, y)}),
          w, "(w + z)*(x + y)");

    // x²y² − z²w² → (x·y − w·z)(x·y + w·z) : difference of squares with
    // composite leading/constant coefficients (perfect squares y² and z²w²).
    check(add({mul(P(x, 2), P(y, 2)),
               mul(integer(-1), mul(P(z, 2), P(w, 2)))}),
          x, "(x*y - w*z)*(x*y + w*z)");

    // A 4-variable polynomial with no shared content and no square structure
    // stays put (no false factorization).
    Expr irr = add({mul(a, b), c, d, integer(1)});  // a·b + c + d + 1
    REQUIRE(factor(irr, a)->str() == irr->str());
}

// x^N − 1 = ∏_{d|N} Φ_d(x). The pure binomial is detected and decomposed into
// cyclotomic polynomials (each irreducible over ℚ) instead of being driven
// through the super-linear Kronecker path — factor(x^60−1) goes from >30s to
// effectively instant. Each result must expand back to x^N−1, match SymPy's
// factor, and have exactly τ(N) (number-of-divisors) irreducible factors.
TEST_CASE("factor: cyclotomic x^N-1 (FACTOR-CYCLO-1)",
          "[polys][factor][cyclotomic][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");

    auto divisor_count = [](long N) {
        long c = 0;
        for (long d = 1; d <= N; ++d)
            if (N % d == 0) ++c;
        return c;
    };

    for (long N : {6, 12, 30, 60}) {
        const std::string xN = "x**" + std::to_string(N) + " - 1";
        Expr e = expand(pow(x, integer(N)) - integer(1));
        Expr f = factor(e, x);

        // Value-equivalence: factored form expands back to x^N − 1.
        REQUIRE(oracle.equivalent(f->str(), xN));
        // Matches SymPy's factorization exactly.
        auto sympy_factor = oracle.send({{"op", "factor"}, {"expr", xN}});
        REQUIRE(sympy_factor.ok);
        REQUIRE(oracle.equivalent(f->str(), sympy_factor.result_str()));

        // Factor count: ∏_{d|N} Φ_d has τ(N) distinct irreducible factors,
        // each of multiplicity 1. Cross-check against SymPy's factor_list.
        auto pf = oracle.send({{"op", "polyfactor"},
                               {"expr", xN},
                               {"var", "x"}});
        REQUIRE(pf.ok);
        REQUIRE(pf.result_str() == std::to_string(divisor_count(N)));

        // Our own factor count from the Mul structure.
        std::size_t our_count = (f->type_id() == sympp::TypeId::Mul)
                                    ? f->args().size()
                                    : 1;
        REQUIRE(our_count == static_cast<std::size_t>(divisor_count(N)));
    }
}
