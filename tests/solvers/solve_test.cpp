// solve / linsolve tests.
//
// References:
//   sympy/solvers/tests/test_solveset.py
//   MATLAB Symbolic Math Toolbox §3 — solve(eqn, var), linsolve

#include <stdexcept>

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/matrices/matrix.hpp>
#include <sympp/polys/rootof.hpp>
#include <sympp/solvers/solve.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- solve (univariate polynomial) ----------------------------------------

TEST_CASE("solve: linear equation", "[10][solve]") {
    auto x = symbol("x");
    // 3x - 7 = 0 -> x = 7/3
    auto roots = solve(integer(3) * x - integer(7), x);
    REQUIRE(roots.size() == 1);
    REQUIRE(roots[0] == rational(7, 3));
}

TEST_CASE("solve: lhs == rhs form", "[10][solve]") {
    auto x = symbol("x");
    // 2x + 1 == 5 -> x = 2
    auto roots = solve(integer(2) * x + integer(1), integer(5), x);
    REQUIRE(roots.size() == 1);
    REQUIRE(roots[0] == integer(2));
}

// SOLVE-DEDUP-1: solve returns the DISTINCT solution set. A repeated factor
// ((x+2)², x²·(x−1)) used to yield duplicate roots; SymPy's solve dedupes.
TEST_CASE("solve: distinct roots for repeated factors (SOLVE-DEDUP-1)",
          "[10][solve][regression]") {
    auto x = symbol("x");
    // (x+2)² → {−2}, not {−2, −2}.
    REQUIRE(solve(pow(x + integer(2), integer(2)), x).size() == 1);
    // x²·(x−1) → {0, 1}.
    REQUIRE(solve(pow(x, integer(2)) * (x - integer(1)), x).size() == 2);
    // (x−1)²·(x+1) → {1, −1}.
    REQUIRE(solve(pow(x - integer(1), integer(2)) * (x + integer(1)), x).size()
            == 2);
    // (x−1)³ → {1}.
    REQUIRE(solve(pow(x - integer(1), integer(3)), x).size() == 1);
    // Distinct roots are unaffected: (x−1)(x−2)(x−3) → 3 roots.
    REQUIRE(solve((x - integer(1)) * (x - integer(2)) * (x - integer(3)), x)
                .size() == 3);
}

TEST_CASE("solve: quadratic with rational roots", "[10][solve]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // x^2 - 5x + 6 = 0 -> x ∈ {2, 3}
    auto roots = solve(pow(x, integer(2)) - integer(5) * x + integer(6), x);
    REQUIRE(roots.size() == 2);
    REQUIRE(oracle.equivalent(roots[0]->str(), "3"));
    REQUIRE(oracle.equivalent(roots[1]->str(), "2"));
}

TEST_CASE("solve: quadratic with irrational roots", "[10][solve][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // x^2 - 2 = 0 -> ±sqrt(2)
    auto roots = solve(pow(x, integer(2)) - integer(2), x);
    REQUIRE(roots.size() == 2);
    REQUIRE(oracle.equivalent(roots[0]->str(), "sqrt(2)"));
    REQUIRE(oracle.equivalent(roots[1]->str(), "-sqrt(2)"));
}

TEST_CASE("solve: quadratic with parametric coefficient",
          "[10][solve][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto a = symbol("a");
    auto b = symbol("b");
    auto c = symbol("c");
    // a*x^2 + b*x + c = 0
    auto roots = solve(a * pow(x, integer(2)) + b * x + c, x);
    REQUIRE(roots.size() == 2);
    INFO("root[0]: " << roots[0]->str());
    INFO("root[1]: " << roots[1]->str());
    REQUIRE(oracle.equivalent(roots[0]->str(),
                              "(-b + sqrt(b**2 - 4*a*c)) / (2*a)"));
    REQUIRE(oracle.equivalent(roots[1]->str(),
                              "(-b - sqrt(b**2 - 4*a*c)) / (2*a)"));
}

TEST_CASE("solve: cubic via rational-root deflation", "[10][solve][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    // x^3 - 8 = (x - 2)(x² + 2x + 4); rational root 2 + two from quadratic.
    auto roots = solve(pow(x, integer(3)) - integer(8), x);
    REQUIRE(roots.size() == 3);
    // Each root must satisfy x³ - 8 = 0.
    for (const auto& r : roots) {
        auto val = pow(r, integer(3)) - integer(8);
        REQUIRE(oracle.equivalent(val->str(), "0"));
    }
}

// ----- transcendental solve (regression, issue #11 / SOLVE-1) ----------------
// solve() used to be polynomial-only and returned [] for equations that
// solveset solves via the _invert chain (log/exp/sinh/…). It now falls back
// to solveset and surfaces a finite solution set.
TEST_CASE("solve: log(x) - 1 = 0 -> x = E", "[10][solve][transcendental][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto roots = solve(log(x) - integer(1), x);
    REQUIRE(roots.size() == 1);
    REQUIRE(oracle.equivalent(roots[0]->str(), "E"));
}

TEST_CASE("solve: exp(x) - 2 = 0 -> x = log(2)", "[10][solve][transcendental][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto roots = solve(exp(x) - integer(2), x);
    REQUIRE(roots.size() == 1);
    REQUIRE(oracle.equivalent(roots[0]->str(), "log(2)"));
}

// SOLVE-RATIONAL-1: solve() of a rational equation clears the denominator,
// solves the numerator, and discards roots that vanish the denominator (poles).
// These reached solve() empty before, since the polynomial path can't build a
// Poly from a 1/x term. Matches SymPy.
TEST_CASE("solve: rational equations clear denominators (SOLVE-RATIONAL-1)",
          "[10][solve][oracle][regression]") {
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
    auto inv = [](const Expr& e) { return pow(e, integer(-1)); };
    // x + 1/x - 2 = 0 → (x-1)²/x → x = 1.
    REQUIRE(set_equal(solve(x + inv(x) - integer(2), x), {"1"}));
    // 1/x - 1/2 = 0 → x = 2.
    REQUIRE(set_equal(solve(inv(x) - rational(1, 2), x), {"2"}));
    // 1/(x-1) + 1/(x+1) = 0 → x = 0.
    REQUIRE(set_equal(
        solve(inv(x - integer(1)) + inv(x + integer(1)), x), {"0"}));
    // (x²-1)/(x-1): x = 1 is a removable pole — only x = -1 survives.
    REQUIRE(set_equal(
        solve((pow(x, integer(2)) - integer(1)) * inv(x - integer(1)), x),
        {"-1"}));
    // 1/(x+1) - 1/(x-1) = 0 has no solution (constant numerator).
    REQUIRE(solve(inv(x + integer(1)) - inv(x - integer(1)), x).empty());
    // 2/(x-1) - 3/(x+2) = 0 → x = 7.
    REQUIRE(set_equal(
        solve(integer(2) * inv(x - integer(1))
                  - integer(3) * inv(x + integer(2)),
              x),
        {"7"}));
}

// SOLVE-EXPLOG-POLY-1: solve() of a polynomial in a single exp or log atom.
// Substitute u = g(x), solve the polynomial, invert each root: exp(x)=c → log(c)
// (skipping c=0), log(x)=c → exp(c). Matches SymPy as a set, including the
// complex root from exp(x)=−1.
TEST_CASE("solve: polynomial in exp/log (SOLVE-EXPLOG-POLY-1)",
          "[10][solve][transcendental][oracle][regression]") {
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
    // exp(x)² − 3exp(x) + 2 → exp(x) = 1, 2 → x = 0, log(2).
    REQUIRE(set_equal(
        solve(pow(exp(x), integer(2)) - integer(3) * exp(x) + integer(2), x),
        {"0", "log(2)"}));
    // exp(x)² − 1 → exp(x) = ±1 → x = 0, iπ (complex root included).
    REQUIRE(set_equal(solve(pow(exp(x), integer(2)) - integer(1), x),
                      {"0", "I*pi"}));
    // log(x)² − 1 → log(x) = ±1 → x = e, 1/e.
    REQUIRE(set_equal(solve(pow(log(x), integer(2)) - integer(1), x),
                      {"E", "exp(-1)"}));
    // log(x)² − 3log(x) + 2 → x = e, e².
    REQUIRE(set_equal(
        solve(pow(log(x), integer(2)) - integer(3) * log(x) + integer(2), x),
        {"E", "exp(2)"}));
    // Scaled log argument: log(2x) = 1 → x = e/2.
    REQUIRE(set_equal(solve(log(integer(2) * x) - integer(1), x), {"E/2"}));
}

// SOLVE-EXPSUM-1: a (Laurent) polynomial in exp(x) whose terms are exp(m·x) for
// integer m — exp(x)+exp(-x)-2, exp(2x)-3·exp(x)+2. Substitute u = exp(x)
// (exp(m·x) → uᵐ), solve in u (clearing negative powers), invert via x = log(u).
// Scaled exponents contribute every complex representative, matching SymPy
// (exp(2x)=1 → {0, iπ}; exp(3x)=1 → the three cube-root logs).
TEST_CASE("solve: Laurent polynomial in exp(x) (SOLVE-EXPSUM-1)",
          "[10][solve][transcendental][oracle][regression]") {
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
    auto ex = [&](const Expr& a) { return exp(a); };
    // 2cosh(x) = 2 → x = 0.
    REQUIRE(set_equal(solve(ex(x) + ex(-x) - integer(2), x), {"0"}));
    // 2sinh(x) = 0 → x = 0, iπ.
    REQUIRE(set_equal(solve(ex(x) - ex(-x), x), {"0", "I*pi"}));
    // Composite exponent: exp(2x) − 3exp(x) + 2 → u=1,2 → x = 0, log(2).
    REQUIRE(set_equal(
        solve(ex(integer(2) * x) - integer(3) * ex(x) + integer(2), x),
        {"0", "log(2)"}));
    // exp(2x) − 5exp(x) + 6 → x = log(2), log(3).
    REQUIRE(set_equal(
        solve(ex(integer(2) * x) - integer(5) * ex(x) + integer(6), x),
        {"log(2)", "log(3)"}));
    // Scaled exponent gives every root: exp(2x)=1 → {0, iπ}.
    REQUIRE(set_equal(solve(ex(integer(2) * x) - integer(1), x), {"0", "I*pi"}));
    // exp(3x)=1 → the three cube roots of unity.
    REQUIRE(set_equal(solve(ex(integer(3) * x) - integer(1), x),
                      {"0", "2*I*pi/3", "-2*I*pi/3"}));
}

// SOLVE-LAMBERT-1: the canonical Lambert-W form a·x·exp(b·x)+c=0 inverts to
// x = W(−c·b/a)/b via W(z)·exp(W(z))=z. Previously returned []. Matches SymPy.
TEST_CASE("solve: Lambert-W equations (SOLVE-LAMBERT-1)",
          "[10][solve][transcendental][oracle][regression]") {
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
    // x·eˣ = c → W(c).
    REQUIRE(set_equal(solve(x * exp(x) - integer(1), x), {"LambertW(1)"}));
    REQUIRE(set_equal(solve(x * exp(x) - integer(2), x), {"LambertW(2)"}));
    // x·eˣ = −1 → W(−1).
    REQUIRE(set_equal(solve(x * exp(x) + integer(1), x), {"LambertW(-1)"}));
    // Scaled exponent: x·e^(2x) = 1 → W(2)/2.
    REQUIRE(set_equal(solve(x * exp(integer(2) * x) - integer(1), x),
                      {"LambertW(2)/2"}));
    // Leading coefficient: 2·x·eˣ = 1 → W(1/2).
    REQUIRE(set_equal(solve(integer(2) * x * exp(x) - integer(1), x),
                      {"LambertW(1/2)"}));
    // Both: x·e^(3x) = 5 → W(15)/3.
    REQUIRE(set_equal(solve(x * exp(integer(3) * x) - integer(5), x),
                      {"LambertW(15)/3"}));
    // Homogeneous: x·eˣ = 0 → W(0) = 0.
    REQUIRE(set_equal(solve(x * exp(x), x), {"0"}));
    // Product-log: x·log(x) − c → exp(W(c)).
    REQUIRE(set_equal(solve(x * log(x) - integer(1), x), {"exp(LambertW(1))"}));
    REQUIRE(set_equal(solve(x * log(x) + integer(1), x), {"exp(LambertW(-1))"}));
    // Additive-log: x + log(x) + c → W(e^(−c)); x+log(x)−1 auto-evaluates to 1.
    REQUIRE(set_equal(solve(x + log(x), x), {"LambertW(1)"}));
    REQUIRE(set_equal(solve(x + log(x) - integer(1), x), {"1"}));
    // Additive-exp: x + eˣ + c → −c − W(e^(−c)); x+eˣ−1 auto-evaluates to 0.
    REQUIRE(set_equal(solve(x + exp(x), x), {"-LambertW(1)"}));
    REQUIRE(set_equal(solve(x + exp(x) - integer(1), x), {"0"}));
    // SOLVE-LAMBERT-2: a non-unit coefficient on the bare-var term, a·var + eᵛᵃʳ + c
    // → var = −W(e^(−c/a)/a) − c/a (and the log analogue → W(a·e^(−c))/a). Covers
    // eˣ = x + 2, i.e. eˣ − x − 2, which previously returned [].
    // eˣ − x − 2 → −2 − W(−e^(−2)).
    REQUIRE(set_equal(solve(exp(x) - x - integer(2), x),
                      {"-2 - LambertW(-exp(-2))"}));
    // eˣ − x → −W(−1).
    REQUIRE(set_equal(solve(exp(x) - x, x), {"-LambertW(-1)"}));
    // 2x + eˣ → −W(1/2).
    REQUIRE(set_equal(solve(integer(2) * x + exp(x), x), {"-LambertW(1/2)"}));
    // eˣ − 2x − 1 → −1/2 − W(−e^(−1/2)/2).
    REQUIRE(set_equal(solve(exp(x) - integer(2) * x - integer(1), x),
                      {"-1/2 - LambertW(-exp(-1/2)/2)"}));
    // Generalized additive-log: 2x + log(x) → W(2)/2.
    REQUIRE(set_equal(solve(integer(2) * x + log(x), x), {"LambertW(2)/2"}));
}

// SOLVE-RADPOLY-1: solve() of a polynomial in a radical x^(1/d) — e.g.
// x − √x − 2 (a quadratic in √x). Substitute t = x^(1/d), solve, back-substitute
// x = t^d, and verify each candidate against the original so extraneous roots
// (√x = −1 ⇒ x = 1) are dropped. Matches SymPy.
TEST_CASE("solve: polynomial in a radical x^(1/d) (SOLVE-RADPOLY-1)",
          "[10][solve][oracle][regression]") {
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
    auto sx = pow(x, rational(1, 2));
    // x − √x − 2 = 0 → √x = 2 (√x = −1 extraneous) → x = 4.
    REQUIRE(set_equal(solve(x - sx - integer(2), x), {"4"}));
    // x − 3√x + 2 → √x = 1, 2 → x = 1, 4.
    REQUIRE(set_equal(solve(x - integer(3) * sx + integer(2), x), {"1", "4"}));
    // x + √x − 6 → √x = 2 (√x = −3 extraneous) → x = 4.
    REQUIRE(set_equal(solve(x + sx - integer(6), x), {"4"}));
    // x − √x → √x(√x − 1) = 0 → x = 0, 1.
    REQUIRE(set_equal(solve(x - sx, x), {"0", "1"}));
    // Cube root: x^(1/3) − 2 → x = 8.
    REQUIRE(set_equal(solve(pow(x, rational(1, 3)) - integer(2), x), {"8"}));
    // A plain polynomial is untouched.
    REQUIRE(set_equal(solve(pow(x, integer(2)) - integer(4), x), {"2", "-2"}));
}

// SOLVE-RADISOLATE-1: a single square root √(g(x)) of a non-trivial radicand
// appearing linearly — isolate and square, then filter the extraneous root the
// squaring introduces (numerically, since denesting √((3−√5)/2) symbolically is
// beyond simplify). solve_radical_poly only handles a polynomial in x^(1/d) of
// the bare variable, so these previously returned [].
TEST_CASE("solve: single √(g(x)) by isolate-and-square (SOLVE-RADISOLATE-1)",
          "[10][solve][oracle][regression]") {
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
    auto sq = [&](const Expr& g) { return pow(g, rational(1, 2)); };
    // √(x+1) − x + 1 = 0 → x = 3 (x = 0 extraneous).
    REQUIRE(set_equal(solve(sq(x + integer(1)) - x + integer(1), x), {"3"}));
    // √(2x+3) − x = 0 → x = 3 (x = −1 extraneous).
    REQUIRE(set_equal(solve(sq(integer(2) * x + integer(3)) - x, x), {"3"}));
    // √(x²+1) − x − 1 = 0 → x = 0.
    REQUIRE(set_equal(
        solve(sq(pow(x, integer(2)) + integer(1)) - x - integer(1), x), {"0"}));
    // √(x+1) + x = 0 → x = (1−√5)/2 (the (1+√5)/2 branch is extraneous).
    REQUIRE(set_equal(solve(sq(x + integer(1)) + x, x), {"1/2 - sqrt(5)/2"}));
    // √(x+1) − x − 1 = 0 → x = −1, 0.
    REQUIRE(set_equal(solve(sq(x + integer(1)) - x - integer(1), x), {"-1", "0"}));
    // No real solution: √(x+1) + 2 = 0 (a positive root would be needed).
    REQUIRE(solve(sq(x + integer(1)) + integer(2), x).empty());
    // The radical = constant case still works.
    REQUIRE(set_equal(solve(sq(x + integer(1)) - integer(2), x), {"3"}));
}

// SOLVE-RADISOLATE-2: a sum/difference of TWO square roots, √g1 ± √g2 = c. Isolate
// one radical and square once, leaving a single radical the same path then clears;
// candidates are filtered against the original to drop the roots squaring adds.
// Previously returned []. Matches SymPy.
TEST_CASE("solve: two square roots by isolate-and-square (SOLVE-RADISOLATE-2)",
          "[10][solve][oracle][regression]") {
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
    auto sq = [&](const Expr& g) { return pow(g, rational(1, 2)); };
    // √x + √(x+1) = 3 → x = 16/9.
    REQUIRE(set_equal(solve(sq(x) + sq(x + integer(1)) - integer(3), x),
                      {"16/9"}));
    // √(2x+1) − √x = 1 → x = 0, 4 (both valid).
    REQUIRE(set_equal(
        solve(sq(integer(2) * x + integer(1)) - sq(x) - integer(1), x),
        {"0", "4"}));
    // √(x−1) + √(x+4) = 5 → x = 5.
    REQUIRE(set_equal(
        solve(sq(x - integer(1)) + sq(x + integer(4)) - integer(5), x), {"5"}));
    // √(x+1) − √(x−1) = 1 → x = 5/4.
    REQUIRE(set_equal(
        solve(sq(x + integer(1)) - sq(x - integer(1)) - integer(1), x),
        {"5/4"}));
}

// SOLVE-EXPBASE-1: constant-base exponential a^x = c (a^x is a Pow with a
// numeric base, so it skips the exp/log transcendental path) → x = log(c)/log(a),
// with an exact integer exponent when c is a power of a (2^x=8 → 3). Restricted
// to a non-perfect-power integer base and bare exponent so the answer matches
// SymPy without enumerating complex representatives.
TEST_CASE("solve: constant-base exponential a^x = c (SOLVE-EXPBASE-1)",
          "[10][solve][transcendental][oracle][regression]") {
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
    auto powx = [&](int a) { return pow(integer(a), x); };
    // Exact integer exponents.
    REQUIRE(set_equal(solve(powx(2) - integer(8), x), {"3"}));
    REQUIRE(set_equal(solve(powx(3) - integer(9), x), {"2"}));
    REQUIRE(set_equal(solve(powx(10) - integer(100), x), {"2"}));
    REQUIRE(set_equal(solve(powx(2) - integer(1), x), {"0"}));   // 2^0 = 1
    REQUIRE(set_equal(solve(powx(6) - integer(36), x), {"2"}));
    // Non-power RHS → log ratio.
    REQUIRE(set_equal(solve(powx(5) - integer(3), x), {"log(3)/log(5)"}));
    REQUIRE(set_equal(solve(powx(2) - integer(7), x), {"log(7)/log(2)"}));
    // Negative RHS → the complex log form, matching SymPy.
    REQUIRE(set_equal(solve(powx(2) + integer(7), x),
                      {"(log(7) + I*pi)/log(2)"}));
    // a^x = 0 has no solution.
    REQUIRE(solve(powx(2), x).empty());
    // Perfect-power base / scaled exponent now solve via solve_const_base_exp_sum
    // (4^x = 2 → 1/2, 2^(2x) = 8 → 3/2), matching SymPy.
    REQUIRE(set_equal(solve(pow(integer(4), x) - integer(2), x), {"1/2"}));
    REQUIRE(set_equal(solve(pow(integer(2), integer(2) * x) - integer(8), x),
                      {"3/2"}));
}

// SOLVE-EXPBASE-SUM-1: sums of constant-base exponentials. Each term reduces to
// coeff·exp(rate·x); commensurate rates substitute u=exp(r₀x) (a polynomial in
// u), two incommensurate rates use the ratio method. Previously returned [].
TEST_CASE("solve: sum of constant-base exponentials (SOLVE-EXPBASE-SUM-1)",
          "[10][solve][transcendental][oracle][regression]") {
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
    auto px = [&](int a) { return pow(integer(a), x); };
    // Incommensurate bases (ratio method): 2^x − 3^x = 0 → 0; 5^x − 2^x = 0 → 0.
    REQUIRE(set_equal(solve(px(2) - px(3), x), {"0"}));
    REQUIRE(set_equal(solve(px(5) - px(2), x), {"0"}));
    // Commensurate (polynomial in u = 2^x): 2^(2x) − 5·2^x + 4 → {0, 2}.
    REQUIRE(set_equal(
        solve(pow(integer(2), integer(2) * x)
                  - integer(5) * px(2) + integer(4),
              x),
        {"0", "2"}));
    // 2^(x+1) − 8 → 2; product of bases 2^x·3^x − 6 = 6^x − 6 → 1.
    REQUIRE(set_equal(solve(pow(integer(2), x + integer(1)) - integer(8), x),
                      {"2"}));
    REQUIRE(set_equal(solve(px(2) * px(3) - integer(6), x), {"1"}));
    // No real solution: 2^x − 3^x + something with negative ratio handled by the
    // ratio sign check — 4^x − 2^(x+1) → 1.
    REQUIRE(set_equal(
        solve(pow(integer(4), x) - pow(integer(2), x + integer(1)), x), {"1"}));
}

// SOLVE-RAD-1: radical equations g^p = c (p a non-integer rational) invert to
// g = c^(1/p). The polynomial path can't see through a fractional power, so
// these used to come back empty. Matches SymPy, including the principal-branch
// convention (negative real RHS → no solution) and integer powers still
// returning every root.
TEST_CASE("solve: radical equations via fractional-power invert (SOLVE-RAD-1)",
          "[10][solve][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto half = rational(1, 2);
    // sqrt(x) - 2 = 0 -> x = 4.
    auto r1 = solve(pow(x, half) - integer(2), x);
    REQUIRE(r1.size() == 1);
    REQUIRE(oracle.equivalent(r1[0]->str(), "4"));
    // x^(1/3) - 2 = 0 -> x = 8.
    auto r2 = solve(pow(x, rational(1, 3)) - integer(2), x);
    REQUIRE(r2.size() == 1);
    REQUIRE(oracle.equivalent(r2[0]->str(), "8"));
    // sqrt(x+1) - 2 = 0 -> x = 3 (recursion through the inner argument).
    auto r3 = solve(pow(x + integer(1), half) - integer(2), x);
    REQUIRE(r3.size() == 1);
    REQUIRE(oracle.equivalent(r3[0]->str(), "3"));
    // Symbolic RHS: sqrt(x) - y = 0 -> x = y**2.
    auto r4 = solve(pow(x, half) - y, x);
    REQUIRE(r4.size() == 1);
    REQUIRE(oracle.equivalent(r4[0]->str(), "y**2"));
    // Principal-branch: sqrt(x) + 2 = 0 has no (real) solution.
    REQUIRE(solve(pow(x, half) + integer(2), x).empty());
    // Integer powers are untouched: x**2 - 4 still returns BOTH roots.
    REQUIRE(solve(pow(x, integer(2)) - integer(4), x).size() == 2);
}

// SOLVE-VAR-1: solve must never return a "solution" that still contains the
// variable being solved for. solve_poly treats a var-dependent coefficient
// (exp(x) in x·exp(x) − 1) as constant and used to hand back the rearrangement
// x = exp(x)**(-1); that is not a solution. The canonical x·exp(x) form now
// returns a (var-free) LambertW root (SOLVE-LAMBERT-1); the other forms here
// stay unsolved — either way, never a var-containing value.
TEST_CASE("solve: never returns a var-dependent rearrangement (SOLVE-VAR-1)",
          "[10][solve][regression]") {
    auto x = symbol("x");
    auto y = symbol("y");
    // Mixed: x·exp(x)−c solves via LambertW; exp(x)+x and x·log(x)−1 don't yet.
    for (const auto& e : {x * exp(x) - integer(1), x * exp(x) - integer(2),
                          exp(x) + x, x * log(x) - integer(1)}) {
        auto roots = solve(e, x);
        for (const auto& r : roots) {
            INFO("offending root: " << r->str());
            REQUIRE_FALSE(has(r, x));
        }
    }
    // Genuine solves are unaffected (sanity guard around the filter).
    REQUIRE(solve(pow(x, integer(2)) - integer(1), x).size() == 2);
    auto par = solve(pow(x, integer(2)) - y, x);  // x = ±sqrt(y): free of x
    REQUIRE(par.size() == 2);
    for (const auto& r : par) REQUIRE_FALSE(has(r, x));
}

TEST_CASE("solve: sinh(x) - 1 = 0 -> x = asinh(1)", "[10][solve][transcendental][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto roots = solve(sinh(x) - integer(1), x);
    REQUIRE(roots.size() == 1);
    REQUIRE(oracle.equivalent(roots[0]->str(), "asinh(1)"));
}

// Guard against the solve <-> solveset recursion: a polynomial must still be
// solved by the (unchanged) polynomial path, and a transcendental fallback
// must terminate.
TEST_CASE("solve: polynomial path intact after transcendental fallback",
          "[10][solve][regression]") {
    auto x = symbol("x");
    auto roots = solve(pow(x, integer(2)) - integer(5) * x + integer(6), x);
    REQUIRE(roots.size() == 2);  // {2, 3} — no regression, no recursion
}

// ----- linsolve --------------------------------------------------------------

TEST_CASE("linsolve: 2x2 system", "[10][linsolve]") {
    // x + 2y = 5
    // 3x + 4y = 11
    // -> x = 1, y = 2
    Matrix A = {{integer(1), integer(2)}, {integer(3), integer(4)}};
    Matrix b = {{integer(5)}, {integer(11)}};
    Matrix x = linsolve(A, b);
    REQUIRE(x.rows() == 2);
    REQUIRE(x.cols() == 1);
    REQUIRE(x.at(0, 0) == integer(1));
    REQUIRE(x.at(1, 0) == integer(2));
}

TEST_CASE("linsolve: singular A throws", "[10][linsolve]") {
    Matrix A = {{integer(1), integer(2)}, {integer(2), integer(4)}};
    Matrix b = {{integer(1)}, {integer(2)}};
    REQUIRE_THROWS(linsolve(A, b));
}

TEST_CASE("linsolve: identity returns b unchanged", "[10][linsolve]") {
    Matrix A = Matrix::identity(3);
    Matrix b = {{integer(7)}, {integer(11)}, {integer(13)}};
    Matrix x = linsolve(A, b);
    REQUIRE(x.equals(b));
}

#include <sympp/core/boolean.hpp>
#include <sympp/sets/sets.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/parsing/parser.hpp>

using sympp::parsing::parse;

// SOLVE-EQ-1: solve() accepts an equation Eq(lhs, rhs), reducing it to
// solve(lhs − rhs). Previously an Eq node fell through to the polynomial path
// as an opaque function and returned []. The parser now builds Eq/Ne/Lt/… as
// Relational nodes, and solve() reduces the Eq case.
TEST_CASE("solve: equation form Eq(lhs, rhs) (SOLVE-EQ-1)",
          "[10][solve][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto set_eq = [&](std::vector<Expr> got,
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
    // Eq(x², 4) → {2, −2}.
    REQUIRE(set_eq(solve(eq(pow(x, integer(2)), integer(4)), x), {"2", "-2"}));
    // Eq(x³, x) → {0, 1, −1}.
    REQUIRE(set_eq(solve(eq(pow(x, integer(3)), x), x), {"0", "1", "-1"}));
    // Eq(2x+1, 5) → {2}.
    REQUIRE(set_eq(solve(eq(integer(2) * x + integer(1), integer(5)), x),
                   {"2"}));
    // The parser builds Eq as a Relational, so solving the parsed string works.
    REQUIRE(set_eq(solve(parse("Eq(x**2, 2*x)"), x), {"0", "2"}));
}

// SOLVE-LOGSUM-1: solve() of a sum of logarithms. Combine Σ cᵢ·log(gᵢ) + K via
// log(∏ gᵢ^cᵢ) = −K ⇒ ∏ gᵢ^cᵢ = exp(−K), solve, and keep only roots in the log
// domain (every gᵢ(root) > 0): log(x)+log(x−1)=0 → x(x−1)=1 → (1+√5)/2 (the
// negative root is dropped). Previously these returned [].
TEST_CASE("solve: sum of logarithms (SOLVE-LOGSUM-1)",
          "[10][solve][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto set_eq = [&](std::vector<Expr> got,
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
    // log(x)+log(x−1)=0 → only the in-domain root (1+√5)/2.
    REQUIRE(set_eq(solve(log(x) + log(x - integer(1)), x),
                   {"1/2 + sqrt(5)/2"}));
    // log(x)+log(x+1)=log(6) → {2} (the root −3 is out of domain).
    REQUIRE(set_eq(
        solve(log(x) + log(x + integer(1)) - log(integer(6)), x), {"2"}));
    // 2·log(x)−log(x+2)=0 → x²=x+2 → {2} (−1 out of domain).
    REQUIRE(set_eq(
        solve(integer(2) * log(x) - log(x + integer(2)), x), {"2"}));
    // log(x+1)+log(x−1)=0 → x²−1=1 → √2 (−√2 out of domain).
    REQUIRE(set_eq(solve(log(x + integer(1)) + log(x - integer(1)), x),
                   {"sqrt(2)"}));
}

// ----- solveset --------------------------------------------------------------

TEST_CASE("solveset: x^2 - 4 = 0 → {2, -2}", "[10][solveset]") {
    auto x = symbol("x");
    auto s = solveset(pow(x, integer(2)) - integer(4), x);
    REQUIRE(s->kind() == SetKind::FiniteSet);
    auto fs = std::static_pointer_cast<const FiniteSet>(s);
    REQUIRE(fs->size() == 2);
}

// SOLVE-TRIG-1: solve() of a single trig equation A*f(B*x)+C=0 returns
// representative roots over one principal period (matching SymPy's solve),
// where solveset would yield an infinite ImageSet. Compared as a set, since
// the root order is presentation-dependent.
TEST_CASE("solve: representative roots of trig equations (SOLVE-TRIG-1)",
          "[10][solve][trig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto set_equal = [&](const std::vector<Expr>& got,
                         const std::vector<std::string>& want) {
        if (got.size() != want.size()) return false;
        // Every wanted root must be matched (oracle-equivalent) by a distinct
        // returned root.
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
    REQUIRE(set_equal(solve(sin(x), x), {"0", "pi"}));
    REQUIRE(set_equal(solve(cos(x), x), {"pi/2", "3*pi/2"}));
    REQUIRE(set_equal(solve(tan(x), x), {"0"}));
    REQUIRE(set_equal(solve(sin(x) - integer(1), x), {"pi/2"}));
    REQUIRE(set_equal(solve(cos(x) + integer(1), x), {"pi"}));
    REQUIRE(set_equal(solve(integer(2) * sin(x) - integer(1), x),
                      {"pi/6", "5*pi/6"}));
    REQUIRE(set_equal(solve(cos(x) - integer(1), x), {"0", "2*pi"}));
    // scaled argument: invert then divide by B
    REQUIRE(set_equal(solve(sin(integer(2) * x), x), {"0", "pi/2"}));
    REQUIRE(set_equal(solve(cos(integer(3) * x), x), {"pi/6", "pi/2"}));
    REQUIRE(set_equal(solve(tan(x) - integer(1), x), {"pi/4"}));
    // homogeneous powers f(x)^n = 0 reduce to f(x) = 0
    REQUIRE(set_equal(solve(pow(sin(x), integer(2)), x), {"0", "pi"}));
    REQUIRE(set_equal(solve(pow(cos(x), integer(2)), x), {"pi/2", "3*pi/2"}));
    REQUIRE(set_equal(solve(pow(sin(x), integer(3)), x), {"0", "pi"}));
    REQUIRE(set_equal(solve(pow(tan(x), integer(2)), x), {"0"}));
    REQUIRE(set_equal(
        solve(integer(2) * pow(sin(x), integer(2)), x), {"0", "pi"}));
    // zero-product: a product of trig factors vanishes iff some factor does
    REQUIRE(set_equal(solve(sin(x) * cos(x), x),
                      {"0", "pi/2", "pi", "3*pi/2"}));
    REQUIRE(set_equal(solve(sin(x) * (cos(x) - integer(1)), x),
                      {"0", "pi", "2*pi"}));
    REQUIRE(set_equal(solve(sin(x) * tan(x), x), {"0", "pi"}));
    REQUIRE(set_equal(
        solve((sin(x) - integer(1)) * (cos(x) + integer(1)), x),
        {"pi/2", "pi"}));
    REQUIRE(set_equal(solve(sin(x) * cos(x) * tan(x), x),
                      {"0", "pi/2", "pi", "3*pi/2"}));
}

// SOLVE-ZEROPROD-1: zero-product over factors that mix polynomial and
// transcendental parts, including an Add with a common factor. solve_poly reads
// x·cos(x) as a linear polynomial with coefficient cos(x) (partial root {0}), and
// the common-factor Add x²·eˣ − eˣ as non-polynomial ({}). Zero-product solves
// each factor and unions, skipping the never-zero eˣ (which also removes the
// spurious zoo that solving eˣ = 0 injected).
TEST_CASE("solve: zero-product with transcendental factors (SOLVE-ZEROPROD-1)",
          "[10][solve][oracle][regression]") {
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
    // Add with a common eˣ factor: eˣ·(x²−1). eˣ never zero → roots of x²−1.
    REQUIRE(set_equal(
        solve(pow(x, integer(2)) * exp(x) - exp(x), x), {"1", "-1"}));
    // Mul with a never-zero factor: no spurious zoo from eˣ = 0.
    REQUIRE(set_equal(
        solve(exp(x) * (pow(x, integer(2)) - integer(4)), x), {"2", "-2"}));
    // Common factor that IS solvable (sin): both it and the cofactor contribute.
    REQUIRE(set_equal(solve(pow(x, integer(2)) * sin(x) - sin(x), x),
                      {"0", "pi", "1", "-1"}));
    // Common x·eˣ factor: x³·eˣ − x·eˣ = x·eˣ·(x²−1).
    REQUIRE(set_equal(
        solve(pow(x, integer(3)) * exp(x) - x * exp(x), x), {"0", "1", "-1"}));
    // Polynomial × trig product, complete root set (solve_poly gave only {0}).
    REQUIRE(set_equal(solve(x * cos(x), x), {"0", "pi/2", "3*pi/2"}));
    REQUIRE(set_equal(solve(sin(x) * (x - integer(1)), x), {"0", "1", "pi"}));
    // Several factors, mixed: eˣ·(x²−1)·(x−3).
    REQUIRE(set_equal(
        solve(exp(x) * (pow(x, integer(2)) - integer(1)) * (x - integer(3)), x),
        {"1", "-1", "3"}));
}

// SOLVE-TRIG-PHASE-1: solve() of a trig equation whose argument carries an
// additive phase — sin(x+1), cos(2x+π/3), tan(x+1) — inverts the inner function
// and solves the affine argument B·x+P=θ for x=(θ−P)/B. Previously these
// returned [] because the solver only accepted a bare B·x argument.
TEST_CASE("solve: trig equations with an additive phase (SOLVE-TRIG-PHASE-1)",
          "[10][solve][trig][oracle][regression]") {
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
    const Expr pi = S::Pi();
    // sin(x+1) = 1/2  ->  x+1 = pi/6 or 5pi/6
    REQUIRE(set_equal(solve(sin(x + integer(1)) - rational(1, 2), x),
                      {"pi/6 - 1", "5*pi/6 - 1"}));
    // cos(2x + pi/3) = 0  ->  2x+pi/3 = pi/2 or 3pi/2
    REQUIRE(set_equal(solve(cos(integer(2) * x + pi / integer(3)), x),
                      {"pi/12", "7*pi/12"}));
    // sin(2x + 1) = 0  ->  2x+1 = 0 or pi
    REQUIRE(set_equal(solve(sin(integer(2) * x + integer(1)), x),
                      {"-1/2", "(pi - 1)/2"}));
    // tan(x + 1) = 1  ->  x+1 = pi/4
    REQUIRE(set_equal(solve(tan(x + integer(1)) - integer(1), x),
                      {"pi/4 - 1"}));
    // phase with a polynomial-in-trig form: sin(x+1)^2 = 1/4
    REQUIRE(set_equal(
        solve(pow(sin(x + integer(1)), integer(2)) - rational(1, 4), x),
        {"pi/6 - 1", "5*pi/6 - 1", "-pi/6 - 1", "7*pi/6 - 1"}));
    // bare argument (P = 0) still works through the same path
    REQUIRE(set_equal(solve(integer(2) * sin(x) - integer(1), x),
                      {"pi/6", "5*pi/6"}));
}

// SOLVE-TRIG-POLY-1: solve() of a polynomial in a single trig function
// (sin(x)²−1, 2cos(x)²−cos(x)−1, …) substitutes u = f(B·x), solves the
// polynomial in u, and inverts each in-range root to representative angles —
// matching SymPy's solve as a set. Out-of-range roots (|c|>1 for sin/cos)
// contribute no real solution.
TEST_CASE("solve: polynomial in a single trig function (SOLVE-TRIG-POLY-1)",
          "[10][solve][trig][oracle][regression]") {
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
    // sin² − 1 = 0 → sin = ±1.
    REQUIRE(set_equal(solve(pow(sin(x), integer(2)) - integer(1), x),
                      {"-pi/2", "pi/2", "3*pi/2"}));
    // 2·sin² − 1 = 0 → sin = ±1/√2 (four representatives over [0, 2π)).
    REQUIRE(set_equal(solve(integer(2) * pow(sin(x), integer(2)) - integer(1), x),
                      {"-pi/4", "pi/4", "3*pi/4", "5*pi/4"}));
    // sin² − sin = 0 → sin = 0 or 1.
    REQUIRE(set_equal(solve(pow(sin(x), integer(2)) - sin(x), x),
                      {"0", "pi/2", "pi"}));
    // cos² − 1/4 = 0 → cos = ±1/2.
    REQUIRE(set_equal(
        solve(pow(cos(x), integer(2)) - rational(1, 4), x),
        {"pi/3", "2*pi/3", "4*pi/3", "5*pi/3"}));
    // tan² − 1 = 0 → tan = ±1 (no range limit on tan).
    REQUIRE(set_equal(solve(pow(tan(x), integer(2)) - integer(1), x),
                      {"-pi/4", "pi/4"}));
    // Quadratic in cos with two distinct roots, one of them ±1 (endpoints).
    REQUIRE(set_equal(
        solve(integer(2) * pow(cos(x), integer(2)) - cos(x) - integer(1), x),
        {"0", "2*pi/3", "4*pi/3", "2*pi"}));
    // Scaled argument: sin(2x)² − 1 inverts then divides by B = 2.
    REQUIRE(set_equal(solve(pow(sin(integer(2) * x), integer(2)) - integer(1), x),
                      {"-pi/4", "pi/4", "3*pi/4"}));
    // Out-of-range: sin² = 4 has no real solution → empty.
    REQUIRE(solve(pow(sin(x), integer(2)) - integer(4), x).empty());
}

// SOLVE-TRIG-MULTIANGLE-1: solve() of an equation mixing different integer
// multiples of one angle — sin(x) − cos(2x), cos(2x) + cos(x), sin(2x) − sin(x).
// expand_trig rewrites every multiple angle in terms of sin/cos of the base
// angle, after which the equation is a polynomial in sin θ and cos θ; reducing
// sin² → 1 − cos² yields P(c) + s·Q(c) = 0, solved through the cosine. Each
// returned root must satisfy the original equation (verified via SymPy); the
// representative-angle form is presentation-dependent so the count and validity
// are checked rather than an exact set. Previously these returned [].
TEST_CASE("solve: mixed multiple-angle trig (SOLVE-TRIG-MULTIANGLE-1)",
          "[10][solve][trig][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto satisfies = [&](const Expr& eq, std::size_t n) {
        auto roots = solve(eq, x);
        if (roots.size() != n) return false;
        for (const auto& r : roots) {
            if (has(r, x)) return false;
            if (!oracle.equivalent(subs(eq, x, r)->str(), "0")) return false;
        }
        return true;
    };
    const Expr two = integer(2);
    // cos(2x) reduces to a polynomial in sin x (cos(2x) = 1 − 2 sin²x).
    REQUIRE(satisfies(sin(x) - cos(two * x), 3));
    REQUIRE(satisfies(cos(two * x) - sin(x), 3));
    REQUIRE(satisfies(cos(two * x) + integer(3) * sin(x) - two, 3));
    REQUIRE(satisfies(cos(two * x) - integer(1) + sin(x), 4));
    // cos(2x) reduces to a polynomial in cos x (cos(2x) = 2 cos²x − 1).
    REQUIRE(satisfies(cos(two * x) + cos(x), 3));
    REQUIRE(satisfies(cos(two * x) - cos(x), 4));
    // sin(2x) = 2 sin x cos x — the mixed s·Q(c) factor cases.
    REQUIRE(satisfies(sin(two * x) - sin(x), 4));
    REQUIRE(satisfies(sin(two * x) - cos(x), 4));
    // A triple angle, via the Chebyshev expansion of cos(3x)/sin(3x).
    REQUIRE(satisfies(cos(integer(3) * x) - cos(x), 5));
    // No real solution: |sin x + cos x| ≤ √2 < 5, so the resultant has only
    // complex cosine roots and the result is empty (not a spurious angle).
    REQUIRE(solve(sin(x) + cos(x) - integer(5), x).empty());
}

// SOLVE-INVFN-1: solve() of an inverse trig/hyperbolic equation inverts via the
// forward function — g⁻¹(B·x)=c → x = g(c)/B — and rejects a root outside the
// inverse function's range (asin(x)=2 → []). Matches SymPy exactly.
TEST_CASE("solve: inverse trig/hyperbolic equations (SOLVE-INVFN-1)",
          "[10][solve][transcendental][oracle][regression]") {
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
    REQUIRE(set_equal(solve(asin(x) - integer(1), x), {"sin(1)"}));
    REQUIRE(set_equal(solve(acos(x) - integer(1), x), {"cos(1)"}));
    REQUIRE(set_equal(solve(atan(x) - integer(1), x), {"tan(1)"}));
    REQUIRE(set_equal(solve(asinh(x) - integer(2), x), {"sinh(2)"}));
    REQUIRE(set_equal(solve(acosh(x) - integer(1), x), {"cosh(1)"}));
    REQUIRE(set_equal(solve(atanh(x) - rational(1, 2), x), {"tanh(1/2)"}));
    // Auto-evaluating RHS: asin(x) = π/6 → x = sin(π/6) = 1/2.
    REQUIRE(set_equal(solve(asin(x) - S::Pi() / integer(6), x), {"1/2"}));
    // Scaled argument: atan(2x) = 1 → x = tan(1)/2.
    REQUIRE(set_equal(solve(atan(integer(2) * x) - integer(1), x),
                      {"tan(1)/2"}));
    // Quadratic in asin: asin(x)² = 1 → x = ±sin(1).
    REQUIRE(set_equal(solve(pow(asin(x), integer(2)) - integer(1), x),
                      {"sin(1)", "-sin(1)"}));
    // Out of range → no solution.
    REQUIRE(solve(asin(x) - integer(2), x).empty());       // 2 > π/2
    REQUIRE(solve(acos(x) - integer(4), x).empty());       // 4 > π
    REQUIRE(solve(acosh(x) + integer(1), x).empty());      // c = -1 < 0
    // Symbolic RHS → the formal principal-branch inverse (matching SymPy). The
    // bounded-range check previously rejected a non-numeric c, returning [].
    auto a = symbol("a");
    REQUIRE(set_equal(solve(atan(x) - a, x), {"tan(a)"}));
    REQUIRE(set_equal(solve(asin(x) - a, x), {"sin(a)"}));
    REQUIRE(set_equal(solve(acos(x) - a, x), {"cos(a)"}));
    REQUIRE(set_equal(solve(atanh(x) - a, x), {"tanh(a)"}));
    REQUIRE(set_equal(solve(asinh(x) - a, x), {"sinh(a)"}));
    // Scaled argument with a symbolic RHS: atan(2x) = a → x = tan(a)/2.
    REQUIRE(set_equal(solve(atan(integer(2) * x) - a, x), {"tan(a)/2"}));
}

// SOLVE-TRIG-LINEAR-1: solve() of a linear combination a·sin(B·x)+b·cos(B·x)+c
// (the R-method). Homogeneous (c=0) reduces to tan(B·x)=-b/a (one root);
// otherwise a·sin+b·cos = R·sin(B·x+φ) with R=√(a²+b²), φ=atan2(b,a), giving
// two representatives. Matches SymPy's solve as a set.
TEST_CASE("solve: linear combination a*sin+b*cos (SOLVE-TRIG-LINEAR-1)",
          "[10][solve][trig][oracle][regression]") {
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
    // Homogeneous: tan(x) = -b/a, a single representative.
    REQUIRE(set_equal(solve(sin(x) + cos(x), x), {"-pi/4"}));
    REQUIRE(set_equal(solve(sin(x) - cos(x), x), {"pi/4"}));
    REQUIRE(set_equal(solve(cos(x) - sin(x), x), {"pi/4"}));
    REQUIRE(set_equal(
        solve(pow(integer(3), rational(1, 2)) * sin(x) + cos(x), x), {"-pi/6"}));
    REQUIRE(set_equal(solve(integer(3) * sin(x) + integer(4) * cos(x), x),
                      {"-atan(4/3)"}));
    // Non-homogeneous: R-method gives two representatives.
    REQUIRE(set_equal(solve(sin(x) + cos(x) - integer(1), x), {"0", "pi/2"}));
    // Scaled argument: divide through by B = 2.
    REQUIRE(set_equal(solve(sin(integer(2) * x) + cos(integer(2) * x), x),
                      {"-pi/8"}));
    // Out of range: |c| > sqrt(a²+b²) has no real solution.
    REQUIRE(solve(sin(x) + cos(x) - integer(5), x).empty());
}

// SOLVESET-POW0-1: g^n = 0 (n a positive integer) has the same solution set as
// g = 0. sin(x)² = 0 used to come back EmptySet because the polynomial path
// can't see through sin; it now recurses on the base.
TEST_CASE("solveset: f(x)^n = 0 reduces to f(x) = 0 (SOLVESET-POW0-1)",
          "[10][solveset][regression]") {
    auto x = symbol("x");
    // sin(x)² = 0 → {n·π} (a periodic ImageSet), NOT EmptySet.
    auto s = solveset(pow(sin(x), integer(2)), x);
    REQUIRE(s->kind() != SetKind::Empty);
    REQUIRE(s->kind() == SetKind::ImageSet);
    // tan(x)² = 0 → {n·π}, matching SymPy exactly.
    auto t = solveset(pow(tan(x), integer(2)), x);
    REQUIRE(t->kind() == SetKind::ImageSet);
    // Higher power, same base.
    REQUIRE(solveset(pow(sin(x), integer(4)), x)->kind() == SetKind::ImageSet);
    // A polynomial base is unaffected: (x-1)² = 0 → {1}.
    auto p = solveset(pow(x - integer(1), integer(2)), x);
    REQUIRE(p->kind() == SetKind::FiniteSet);
    REQUIRE(std::static_pointer_cast<const FiniteSet>(p)->size() == 1);
}

// SOLVESET-TRIG-SCALE-1: trig equations with a scaled argument a·x and a
// nonzero RHS (cos(2x)=1) now invert — the periodic image divides through by a.
// Also dedupes the cos union when acos(c) ∈ {0, π} (c = ±1), so cos(x)=1 gives a
// single {2nπ} instead of two identical ImageSets.
TEST_CASE("solveset: scaled-argument trig (SOLVESET-TRIG-SCALE-1)",
          "[10][solveset][regression]") {
    auto x = symbol("x");
    // cos(2x) = 1 → {nπ} (was EmptySet).
    auto s = solveset(cos(integer(2) * x) - integer(1), x);
    REQUIRE(s->kind() == SetKind::ImageSet);
    REQUIRE(s->str().find("__n*pi") != std::string::npos);
    // cos(x) = 1 → a single {2nπ}, not a Union of two identical sets.
    auto c1 = solveset(cos(x) - integer(1), x);
    REQUIRE(c1->kind() == SetKind::ImageSet);
    // cos(x) = -1 → {π + 2nπ} (dedup via acos(-1) = π).
    REQUIRE(solveset(cos(x) + integer(1), x)->kind() == SetKind::ImageSet);
    // tan(2x) = 1 stays a single ImageSet.
    REQUIRE(solveset(tan(integer(2) * x) - integer(1), x)->kind()
            == SetKind::ImageSet);
    // A generic RHS keeps the two-branch union.
    REQUIRE(solveset(cos(x) - rational(1, 2), x)->kind() == SetKind::Union);
}

TEST_CASE("solveset: no solutions → EmptySet",
          "[10][solveset]") {
    auto x = symbol("x");
    // Constant 5 = 0 has no solution.
    auto s = solveset(integer(5), x);
    REQUIRE(s->kind() == SetKind::Empty);
}

TEST_CASE("solveset: filters via domain",
          "[10][solveset]") {
    auto x = symbol("x");
    auto eq = pow(x, integer(2)) - integer(4);
    auto domain = interval(integer(0), integer(10));
    auto s = solveset(eq, x, domain);
    // Only x=2 falls in [0, 10]; x=-2 is filtered out.
    REQUIRE(s->kind() == SetKind::FiniteSet);
    auto fs = std::static_pointer_cast<const FiniteSet>(s);
    REQUIRE(fs->size() == 1);
    REQUIRE(fs->elements()[0] == integer(2));
}

// ----- solveset: _invert chain (transcendental) -----------------------------

#include <sympp/functions/hyperbolic.hpp>
#include <sympp/functions/miscellaneous.hpp>

TEST_CASE("solveset: log(x) = 2 → {exp(2)}",
          "[10][solveset][invert][oracle]") {
    auto x = symbol("x");
    auto eq = log(x) - integer(2);
    auto s = solveset(eq, x);
    REQUIRE(s->kind() == SetKind::FiniteSet);
    auto fs = std::static_pointer_cast<const FiniteSet>(s);
    REQUIRE(fs->size() == 1);
    auto& oracle = sympp::testing::Oracle::instance();
    REQUIRE(oracle.equivalent(fs->elements()[0]->str(), "exp(2)"));
}

TEST_CASE("solveset: exp(x) = 5 → {log(5)}",
          "[10][solveset][invert][oracle]") {
    auto x = symbol("x");
    auto eq = exp(x) - integer(5);
    auto s = solveset(eq, x);
    REQUIRE(s->kind() == SetKind::FiniteSet);
    auto fs = std::static_pointer_cast<const FiniteSet>(s);
    REQUIRE(fs->size() == 1);
    auto& oracle = sympp::testing::Oracle::instance();
    REQUIRE(oracle.equivalent(fs->elements()[0]->str(), "log(5)"));
}

TEST_CASE("solveset: sinh(x) = 3 → {asinh(3)}",
          "[10][solveset][invert][oracle]") {
    auto x = symbol("x");
    auto eq = sinh(x) - integer(3);
    auto s = solveset(eq, x);
    REQUIRE(s->kind() == SetKind::FiniteSet);
    auto fs = std::static_pointer_cast<const FiniteSet>(s);
    REQUIRE(fs->size() == 1);
    auto& oracle = sympp::testing::Oracle::instance();
    REQUIRE(oracle.equivalent(fs->elements()[0]->str(), "asinh(3)"));
}

TEST_CASE("solveset: cosh(x) = 5 → {acosh(5), -acosh(5)}",
          "[10][solveset][invert]") {
    auto x = symbol("x");
    auto eq = cosh(x) - integer(5);
    auto s = solveset(eq, x);
    REQUIRE(s->kind() == SetKind::FiniteSet);
    auto fs = std::static_pointer_cast<const FiniteSet>(s);
    REQUIRE(fs->size() == 2);
}

TEST_CASE("solveset: |x| = 3 → {3, -3}",
          "[10][solveset][invert]") {
    auto x = symbol("x");
    auto eq = abs(x) - integer(3);
    auto s = solveset(eq, x);
    REQUIRE(s->kind() == SetKind::FiniteSet);
    auto fs = std::static_pointer_cast<const FiniteSet>(s);
    REQUIRE(fs->size() == 2);
}

// SOLVE-ABS-1: solve() of an absolute-value equation |g(x)| = c. solveset
// returns {3} ∪ {−1} for |x−1| = 2 (a Union of finite sets); solve() now
// flattens a Union of finite sets into the root list. |g| = c with c < 0 has no
// real solution (the previously-unsound case that the flattening exposed).
TEST_CASE("solve: absolute-value equations (SOLVE-ABS-1)",
          "[10][solve][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto set_eq = [&](std::vector<Expr> got,
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
    // |x−1| = 2 → {3, −1}; |2x−1| = 5 → {3, −2}; |x²−1| = 3 → {2, −2}.
    REQUIRE(set_eq(solve(abs(x - integer(1)) - integer(2), x), {"3", "-1"}));
    REQUIRE(set_eq(solve(abs(integer(2) * x - integer(1)) - integer(5), x),
                   {"3", "-2"}));
    REQUIRE(set_eq(solve(abs(pow(x, integer(2)) - integer(1)) - integer(3), x),
                   {"2", "-2"}));
    // |x| = 0 → {0}.
    REQUIRE(set_eq(solve(abs(x), x), {"0"}));
    // |g| = c with c < 0 has no real solution (soundness).
    REQUIRE(solve(abs(x + integer(1)) + integer(2), x).empty());
    REQUIRE(solve(abs(x) + integer(5), x).empty());
}

TEST_CASE("solveset: sin(x) = 1/2 → ImageSet over ℤ",
          "[10][solveset][invert]") {
    auto x = symbol("x");
    auto eq = sin(x) - rational(1, 2);
    auto s = solveset(eq, x);
    REQUIRE(s->kind() == SetKind::ImageSet);
}

TEST_CASE("solveset: cos(x) = 1/2 → Union of two ImageSets",
          "[10][solveset][invert]") {
    auto x = symbol("x");
    auto eq = cos(x) - rational(1, 2);
    auto s = solveset(eq, x);
    REQUIRE(s->kind() == SetKind::Union);
}

TEST_CASE("solveset: tan(x) = 1 → ImageSet over ℤ",
          "[10][solveset][invert]") {
    auto x = symbol("x");
    auto eq = tan(x) - integer(1);
    auto s = solveset(eq, x);
    REQUIRE(s->kind() == SetKind::ImageSet);
}

// ----- nsolve ----------------------------------------------------------------

TEST_CASE("nsolve: x^3 - 2 ≈ 1.2599",
          "[10][nsolve][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto root = nsolve(pow(x, integer(3)) - integer(2), x, integer(1), 15);
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", root->str() + " - 1.2599210498948732"},
                             {"prec", 30}, {"tol", 10}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

TEST_CASE("nsolve: cos(x) - x ≈ 0.7390851332151607",
          "[10][nsolve][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto root = nsolve(cos(x) - x, x, integer(1), 15);
    auto resp = oracle.send({{"op", "evalf_is_zero"},
                             {"expr", root->str() + " - 0.7390851332151607"},
                             {"prec", 30}, {"tol", 10}});
    REQUIRE(resp.ok);
    REQUIRE(resp.raw.at("result").get<bool>());
}

// ----- solve_univariate_inequality -------------------------------------------

TEST_CASE("inequality: x^2 - 4 < 0 → (-2, 2)",
          "[10][inequality]") {
    auto x = symbol("x");
    auto op = symbol("<");
    auto s = solve_univariate_inequality(
        pow(x, integer(2)) - integer(4), op, S::Zero(), x);
    // The solution should be a single open interval roughly (-2, 2).
    REQUIRE(s->kind() != SetKind::Empty);
    REQUIRE(s->contains(integer(0)) == std::optional<bool>{true});
    REQUIRE(s->contains(integer(3)) == std::optional<bool>{false});
}

TEST_CASE("inequality: x^2 - 4 > 0 outside (-2, 2)",
          "[10][inequality]") {
    auto x = symbol("x");
    auto op = symbol(">");
    auto s = solve_univariate_inequality(
        pow(x, integer(2)) - integer(4), op, S::Zero(), x);
    REQUIRE(s->contains(integer(3)) == std::optional<bool>{true});
    REQUIRE(s->contains(integer(0)) == std::optional<bool>{false});
}

// INEQ-EXACT-1: inequality solutions use EXACT root endpoints (−2, not −2.0…)
// and the real ±∞ (not a 1e30 proxy). solve() also now expands a factored
// polynomial, so (x−1)(x−2) < 0 resolves to (1, 2) instead of EmptySet.
TEST_CASE("inequality: exact endpoints, real infinity (INEQ-EXACT-1)",
          "[10][inequality][regression]") {
    auto x = symbol("x");
    auto lt = symbol("<");
    auto gt = symbol(">");
    auto le = symbol("<=");
    // x²−4 < 0 → (−2, 2) with EXACT integer endpoints.
    {
        auto s = solve_univariate_inequality(pow(x, integer(2)) - integer(4),
                                             lt, S::Zero(), x);
        REQUIRE(s->kind() == SetKind::Interval);
        const auto& iv = static_cast<const Interval&>(*s);
        REQUIRE(iv.lo() == integer(-2));
        REQUIRE(iv.hi() == integer(2));
        REQUIRE(iv.left_open());
        REQUIRE(iv.right_open());
    }
    // x²−4 > 0 → (−∞,−2) ∪ (2,∞): a Union of rays, real infinity.
    REQUIRE(solve_univariate_inequality(pow(x, integer(2)) - integer(4), gt,
                                        S::Zero(), x)
                ->kind() == SetKind::Union);
    // x²+1 > 0 → Reals (no real roots, positive everywhere).
    REQUIRE(solve_univariate_inequality(pow(x, integer(2)) + integer(1), gt,
                                        S::Zero(), x)
                ->kind() == SetKind::Reals);
    // (x−1)(x−2) < 0 → (1, 2) — the solve-expand fix (was EmptySet).
    {
        auto s = solve_univariate_inequality(
            (x - integer(1)) * (x - integer(2)), lt, S::Zero(), x);
        REQUIRE(s->kind() == SetKind::Interval);
        const auto& iv = static_cast<const Interval&>(*s);
        REQUIRE(iv.lo() == integer(1));
        REQUIRE(iv.hi() == integer(2));
    }
    // x²−4 ≤ 0 → [−2, 2] (closed endpoints).
    {
        auto s = solve_univariate_inequality(pow(x, integer(2)) - integer(4),
                                             le, S::Zero(), x);
        const auto& iv = static_cast<const Interval&>(*s);
        REQUIRE_FALSE(iv.left_open());
        REQUIRE_FALSE(iv.right_open());
    }
}

// ----- rsolve ----------------------------------------------------------------

TEST_CASE("rsolve: y(n+1) - 2*y(n) = 0 → y = C * 2^n",
          "[10][rsolve][oracle]") {
    auto& oracle = Oracle::instance();
    auto n = symbol("n");
    // Coeffs: c_0 = -2, c_1 = 1 (i.e., -2 y(n) + y(n+1) = 0)
    auto sol = rsolve({integer(-2), integer(1)}, n);
    // Should be C * 2^n for some symbolic C.
    REQUIRE(oracle.equivalent(sol->str(), "__C0 * 2**n"));
}

TEST_CASE("rsolve: Fibonacci recurrence y(n+2) = y(n+1) + y(n)",
          "[10][rsolve]") {
    auto n = symbol("n");
    // y(n+2) - y(n+1) - y(n) = 0  → coeffs [-1, -1, 1]
    auto sol = rsolve({integer(-1), integer(-1), integer(1)}, n);
    // Solution: C0 * φ^n + C1 * (1-φ)^n with φ = (1+√5)/2.
    REQUIRE(sol->type_id() == TypeId::Add);
}

// ----- nonlinsolve -----------------------------------------------------------

TEST_CASE("nonlinsolve: x^2 + y^2 = 1 ∧ x = y",
          "[10][nonlinsolve][oracle]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto eq1 = pow(x, integer(2)) + pow(y, integer(2)) - integer(1);
    auto eq2 = x - y;
    auto sols = nonlinsolve({eq1, eq2}, {x, y});
    // Two solutions: (sqrt(2)/2, sqrt(2)/2) and (-sqrt(2)/2, -sqrt(2)/2).
    REQUIRE(sols.size() == 2);
}

// ----- linear_diophantine ----------------------------------------------------

TEST_CASE("diophantine: 3x + 5y = 1 has solution",
          "[10][diophantine]") {
    auto sol = linear_diophantine(integer(3), integer(5), integer(1));
    REQUIRE(sol.has_value());
    auto [x, y] = *sol;
    // Verify 3*x + 5*y = 1.
    auto check = integer(3) * x + integer(5) * y;
    REQUIRE(check == integer(1));
}

TEST_CASE("diophantine: 6x + 9y = 5 has no solution (gcd = 3 ∤ 5)",
          "[10][diophantine]") {
    auto sol = linear_diophantine(integer(6), integer(9), integer(5));
    REQUIRE(!sol.has_value());
}

TEST_CASE("diophantine: 4x + 6y = 10 (gcd = 2 ∣ 10)",
          "[10][diophantine]") {
    auto sol = linear_diophantine(integer(4), integer(6), integer(10));
    REQUIRE(sol.has_value());
    auto [x, y] = *sol;
    auto check = integer(4) * x + integer(6) * y;
    REQUIRE(check == integer(10));
}

// ----- Trig solveset → ImageSet ----------------------------------------------

TEST_CASE("solveset: sin(x) = 0 → ImageSet of n*pi",
          "[10d][solveset][trig]") {
    auto x = symbol("x");
    auto s = solveset(sin(x), x);
    REQUIRE(s->kind() == SetKind::ImageSet);
}

TEST_CASE("solveset: cos(x) = 0 → ImageSet of (2n+1)*pi/2",
          "[10d][solveset][trig]") {
    auto x = symbol("x");
    auto s = solveset(cos(x), x);
    REQUIRE(s->kind() == SetKind::ImageSet);
}

TEST_CASE("solveset: sin(2*x) = 0 → ImageSet of n*pi/2",
          "[10d][solveset][trig]") {
    auto x = symbol("x");
    auto s = solveset(sin(integer(2) * x), x);
    REQUIRE(s->kind() == SetKind::ImageSet);
}

TEST_CASE("solveset: tan(x) = 0 → ImageSet of n*pi",
          "[10d][solveset][trig]") {
    auto x = symbol("x");
    auto s = solveset(tan(x), x);
    REQUIRE(s->kind() == SetKind::ImageSet);
}

// ----- reduce_inequalities ---------------------------------------------------

TEST_CASE("reduce_inequalities: single x < 5",
          "[10d][reduce_inequalities]") {
    auto x = symbol("x");
    auto rel = lt(x, integer(5));
    auto s = reduce_inequalities(rel, x);
    REQUIRE(s->contains(integer(0)) == std::optional<bool>{true});
    REQUIRE(s->contains(integer(10)) == std::optional<bool>{false});
}

TEST_CASE("reduce_inequalities: AND combines with intersection",
          "[10d][reduce_inequalities]") {
    auto x = symbol("x");
    // x > 0 AND x < 5
    std::vector<Expr> rels = {gt(x, integer(0)), lt(x, integer(5))};
    auto s = reduce_inequalities(rels, x, true);
    REQUIRE(s->contains(integer(3)) == std::optional<bool>{true});
    REQUIRE(s->contains(integer(-1)) == std::optional<bool>{false});
    REQUIRE(s->contains(integer(10)) == std::optional<bool>{false});
}

TEST_CASE("reduce_inequalities: OR combines with union",
          "[10d][reduce_inequalities]") {
    auto x = symbol("x");
    // x < -5 OR x > 5 — outside [-5, 5]
    std::vector<Expr> rels = {lt(x, integer(-5)), gt(x, integer(5))};
    auto s = reduce_inequalities(rels, x, false);
    REQUIRE(s->contains(integer(-10)) == std::optional<bool>{true});
    REQUIRE(s->contains(integer(10)) == std::optional<bool>{true});
    REQUIRE(s->contains(integer(0)) == std::optional<bool>{false});
}

// ----- Pythagorean Diophantine -----------------------------------------------

TEST_CASE("pythagorean_triples: max_z=20 contains (3,4,5), (5,12,13)",
          "[10d][pythagorean]") {
    auto triples = pythagorean_triples(20);
    REQUIRE(!triples.empty());
    bool has_345 = false, has_5_12_13 = false;
    for (auto& t : triples) {
        if (t[0] == integer(3) && t[1] == integer(4) && t[2] == integer(5)) {
            has_345 = true;
        }
        if (t[0] == integer(5) && t[1] == integer(12) && t[2] == integer(13)) {
            has_5_12_13 = true;
        }
    }
    REQUIRE(has_345);
    REQUIRE(has_5_12_13);
    // Verify x² + y² = z² for each.
    for (auto& t : triples) {
        REQUIRE(pow(t[0], integer(2)) + pow(t[1], integer(2))
                == pow(t[2], integer(2)));
    }
}

// ----- Gröbner basis ---------------------------------------------------------

TEST_CASE("groebner: simple linear system",
          "[10d][groebner]") {
    auto x = symbol("x");
    auto y = symbol("y");
    // x - y, x + y - 2 → unique solution (1, 1).
    auto basis = groebner({x - y, x + y - integer(2)}, {x, y});
    REQUIRE(!basis.empty());
}

TEST_CASE("groebner: x²+y²=1 ∧ x=y returns triangular form",
          "[10d][groebner]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto basis = groebner({pow(x, integer(2)) + pow(y, integer(2)) - integer(1),
                            x - y}, {x, y});
    REQUIRE(!basis.empty());
}

TEST_CASE("nonlinsolve_groebner: 2-variable polynomial system",
          "[10d][nonlinsolve_groebner]") {
    auto x = symbol("x");
    auto y = symbol("y");
    // x*y = 1, x + y = 3 → roots of t² - 3t + 1, two solutions.
    auto sols = nonlinsolve_groebner(
        {x * y - integer(1), x + y - integer(3)}, {x, y});
    REQUIRE(!sols.empty());
}

// SOLVE-ROOTOF-1: an irreducible polynomial of degree >= 5 is not solvable by
// radicals, so the closed-form solver left it unrepresented and solve() returned
// an empty list — implying "no solutions" for x^5 - x - 1, which has a real root
// near 1.1673. solve() now supplements with RootOf (rendered CRootOf) for the
// real roots of such factors. Each returned root is checked numerically by
// substituting its high-precision value back into the polynomial.
TEST_CASE("solve: irreducible quintic returns real CRootOf (SOLVE-ROOTOF-1)",
          "[10][solve][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");

    // x^5 - x - 1: one real root, four complex (the real one is represented).
    {
        Expr p = pow(x, integer(5)) - x - integer(1);
        auto roots = solve(p, x);
        REQUIRE(roots.size() == 1);
        REQUIRE(roots[0]->type_id() == TypeId::RootOf);
        auto v = static_cast<const RootOf&>(*roots[0]).try_evalf(30);
        REQUIRE(v.has_value());
        // p(root) ≈ 0
        Expr resid = subs(p, x, *v);
        auto resp = oracle.send({{"op", "evalf_is_zero"},
                                 {"expr", resid->str()},
                                 {"prec", 40},
                                 {"tol", 20}});
        REQUIRE(resp.ok);
        REQUIRE(resp.raw.at("result").get<bool>());
    }

    // 2*x^5 - 10*x + 5: three real roots, all returned as CRootOf.
    {
        Expr p = integer(2) * pow(x, integer(5)) - integer(10) * x + integer(5);
        auto roots = solve(p, x);
        REQUIRE(roots.size() == 3);
        for (const auto& r : roots) {
            REQUIRE(r->type_id() == TypeId::RootOf);
            auto v = static_cast<const RootOf&>(*r).try_evalf(30);
            REQUIRE(v.has_value());
            Expr resid = subs(p, x, *v);
            auto resp = oracle.send({{"op", "evalf_is_zero"},
                                     {"expr", resid->str()},
                                     {"prec", 40},
                                     {"tol", 20}});
            REQUIRE(resp.ok);
            REQUIRE(resp.raw.at("result").get<bool>());
        }
    }

    // Mixed: (x - 2)·(x^5 - x - 1) keeps the rational root and adds the CRootOf.
    {
        Expr p = (x - integer(2)) * (pow(x, integer(5)) - x - integer(1));
        auto roots = solve(p, x);
        REQUIRE(roots.size() == 2);
        bool has_two = false, has_rootof = false;
        for (const auto& r : roots) {
            if (r == integer(2)) has_two = true;
            if (r->type_id() == TypeId::RootOf) has_rootof = true;
        }
        REQUIRE(has_two);
        REQUIRE(has_rootof);
    }
}
