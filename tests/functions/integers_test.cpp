// floor / ceiling tests.
//
// Reference: sympy/functions/elementary/tests/test_integers.py
//   * test_floor, test_ceiling

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/functions/integers.hpp>
#include <sympp/parsing/parser.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- floor -----------------------------------------------------------------

TEST_CASE("floor: integer is identity", "[3g][floor]") {
    REQUIRE(floor(integer(0)) == integer(0));
    REQUIRE(floor(integer(7)) == integer(7));
    REQUIRE(floor(integer(-3)) == integer(-3));
}

TEST_CASE("floor: rational rounds toward -inf", "[3g][floor]") {
    REQUIRE(floor(rational(7, 2)) == integer(3));     // 3.5 -> 3
    REQUIRE(floor(rational(-7, 2)) == integer(-4));   // -3.5 -> -4
    REQUIRE(floor(rational(1, 3)) == integer(0));     // 0.333.. -> 0
    REQUIRE(floor(rational(-1, 3)) == integer(-1));   // -0.333 -> -1
}

TEST_CASE("floor: float matches mpfr_floor", "[3g][floor]") {
    auto e = floor(float_value(2.7));
    REQUIRE(e->type_id() == TypeId::Float);
}

TEST_CASE("floor: integer-tagged symbol is identity",
          "[3g][floor][assumptions]") {
    auto n = symbol("n", AssumptionMask{}.set_integer(true));
    REQUIRE(floor(n) == n);
}

TEST_CASE("floor: stays unevaluated on a generic symbol", "[3g][floor]") {
    auto x = symbol("x");
    REQUIRE(floor(x)->str() == "floor(x)");
}

TEST_CASE("floor/ceiling: integer-shift invariance (ASSUME-13)",
          "[3g][floor][ceiling][assumptions][regression]") {
    auto n = symbol("n", AssumptionMask{}.set_integer(true));
    auto x = symbol("x");
    // floor(n + x) = n + floor(x); the integer part pulls out.
    REQUIRE(floor(n + rational(1, 2)) == n);            // floor(1/2)=0
    REQUIRE(ceiling(n + rational(1, 2)) == n + integer(1));
    REQUIRE(floor(integer(2) * n + x) == integer(2) * n + floor(x));
    REQUIRE(floor(n + x) == n + floor(x));
    // No integer part to pull → unchanged.
    REQUIRE(floor(x + rational(1, 2))->str().find("floor(") != std::string::npos);
}

// FLOOR-IDEMPOTENT-1: floor/ceiling of an already integer-valued floor/ceiling
// is the identity — floor(floor x)=floor x, floor(ceil x)=ceil x, and the
// cross/mirror pairs. Holds for a generic argument (no reality assumption needed),
// matching SymPy.
TEST_CASE("floor/ceiling: idempotence (FLOOR-IDEMPOTENT-1)",
          "[3g][floor][ceiling][regression]") {
    auto x = symbol("x");  // generic — reality unknown
    REQUIRE(floor(floor(x)) == floor(x));
    REQUIRE(ceiling(ceiling(x)) == ceiling(x));
    REQUIRE(floor(ceiling(x)) == ceiling(x));
    REQUIRE(ceiling(floor(x)) == floor(x));
    // Composes with the integer-shift pull-out: floor(floor(x) + 2) = floor(x) + 2.
    REQUIRE(floor(floor(x) + integer(2)) == floor(x) + integer(2));
    // A non-trivial multiple of a floor is left intact (as SymPy does).
    REQUIRE(floor(integer(2) * floor(x))->str().find("floor(2") != std::string::npos);
}

// Regression (FLOOR-CONST): floor/ceiling of a real constant expression
// evaluate numerically — floor(pi)=3, floor(2*pi)=6, floor(sqrt(2))=1.
TEST_CASE("floor/ceiling: real constants evaluate numerically",
          "[3g][floor][ceiling][regression]") {
    REQUIRE(floor(S::Pi()) == integer(3));
    REQUIRE(ceiling(S::Pi()) == integer(4));
    REQUIRE(floor(S::E()) == integer(2));
    REQUIRE(floor(integer(2) * S::Pi()) == integer(6));
    REQUIRE(floor(mul(integer(-1), S::Pi())) == integer(-4));
    REQUIRE(ceiling(mul(integer(-1), S::Pi())) == integer(-3));
    REQUIRE(floor(pow(integer(2), rational(1, 2))) == integer(1));
    REQUIRE(ceiling(pow(integer(2), rational(1, 2))) == integer(2));
    // Imaginary / infinite arguments must NOT fold.
    REQUIRE(floor(S::I())->str() == "floor(I)");
    REQUIRE(floor(S::Infinity())->str() == "floor(oo)");
}

// ----- ceiling ---------------------------------------------------------------

TEST_CASE("ceiling: integer is identity", "[3g][ceiling]") {
    REQUIRE(ceiling(integer(7)) == integer(7));
    REQUIRE(ceiling(integer(-3)) == integer(-3));
}

TEST_CASE("ceiling: rational rounds toward +inf", "[3g][ceiling]") {
    REQUIRE(ceiling(rational(7, 2)) == integer(4));    // 3.5 -> 4
    REQUIRE(ceiling(rational(-7, 2)) == integer(-3));  // -3.5 -> -3
    REQUIRE(ceiling(rational(1, 3)) == integer(1));
    REQUIRE(ceiling(rational(-1, 3)) == integer(0));
}

TEST_CASE("ceiling: float matches mpfr_ceil", "[3g][ceiling]") {
    auto e = ceiling(float_value(2.3));
    REQUIRE(e->type_id() == TypeId::Float);
}

TEST_CASE("ceiling: stays unevaluated on a generic symbol", "[3g][ceiling]") {
    auto x = symbol("x");
    REQUIRE(ceiling(x)->str() == "ceiling(x)");
}

// ----- Substitution ----------------------------------------------------------

TEST_CASE("floor/ceiling: subs collapses on numeric substitution",
          "[3g][subs]") {
    auto x = symbol("x");
    REQUIRE(subs(floor(x), x, rational(11, 3)) == integer(3));
    REQUIRE(subs(ceiling(x), x, rational(11, 3)) == integer(4));
}

// ----- Oracle structural -----------------------------------------------------

TEST_CASE("floor / ceiling: structural and value match SymPy",
          "[3g][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    REQUIRE(oracle.equivalent(floor(x)->str(), "floor(x)"));
    REQUIRE(oracle.equivalent(ceiling(x)->str(), "ceiling(x)"));

    REQUIRE(oracle.equivalent(floor(rational(7, 2))->str(), "3"));
    REQUIRE(oracle.equivalent(ceiling(rational(7, 2))->str(), "4"));
    REQUIRE(oracle.equivalent(floor(rational(-7, 2))->str(), "-4"));
    REQUIRE(oracle.equivalent(ceiling(rational(-7, 2))->str(), "-3"));
}

// ----- Frac (fractional part, FRAC-1) ----------------------------------------
// frac(x) = x − floor(x), always in [0, 1): frac(7/2)=1/2, frac(-7/2)=1/2,
// frac(int)=0, frac(pi)=pi−3. Symbolic args stay unevaluated.
TEST_CASE("frac: rational and integer args", "[3g][frac]") {
    REQUIRE(frac(rational(7, 2)) == rational(1, 2));
    REQUIRE(frac(rational(-7, 2)) == rational(1, 2));  // in [0,1), not -1/2
    REQUIRE(frac(integer(5)) == integer(0));
    REQUIRE(frac(integer(-3)) == integer(0));
    REQUIRE(frac(rational(1, 3)) == rational(1, 3));
}

TEST_CASE("frac: real constant and symbolic", "[3g][frac][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    REQUIRE(oracle.equivalent(frac(S::Pi())->str(), "-3 + pi"));  // pi - 3
    REQUIRE(frac(x)->type_id() == TypeId::Function);
    REQUIRE(parsing::parse("frac(x)") == frac(x));
    REQUIRE(frac(x)->str() == "frac(x)");
}

// ----- Mod (floored modulo, MOD-1) -------------------------------------------
// Mod(p,q) = p − q·floor(p/q); the result takes the sign of q (floored, not
// truncated): Mod(-7,3)=2, Mod(7,-3)=-2. Integer and rational args evaluate;
// Mod(0,q)=0, Mod(x,x)=0; q==0 stays unevaluated (no throw).
TEST_CASE("mod: integer floored modulo", "[3g][mod]") {
    REQUIRE(mod(integer(7), integer(3)) == integer(1));
    REQUIRE(mod(integer(-7), integer(3)) == integer(2));     // sign of divisor
    REQUIRE(mod(integer(7), integer(-3)) == integer(-2));
    REQUIRE(mod(integer(-7), integer(-3)) == integer(-1));
    REQUIRE(mod(integer(0), integer(3)) == integer(0));
    REQUIRE(mod(integer(10), integer(5)) == integer(0));
}

TEST_CASE("mod: rational args", "[3g][mod]") {
    REQUIRE(mod(rational(1, 2), rational(1, 3)) == rational(1, 6));
}

TEST_CASE("mod: structural identities and zero divisor", "[3g][mod]") {
    auto x = symbol("x");
    REQUIRE(mod(integer(0), x) == integer(0));   // Mod(0, q) = 0
    REQUIRE(mod(x, x) == integer(0));            // Mod(x, x) = 0
    // q == 0 is undefined — kept unevaluated rather than throwing.
    REQUIRE(mod(x, integer(0))->type_id() == TypeId::Function);
    // Generic symbolic args stay unevaluated.
    REQUIRE(mod(x, integer(3))->type_id() == TypeId::Function);
}

TEST_CASE("mod: Mod(integer, 1) = 0 (ASSUME-14)",
          "[3g][mod][assumptions][regression]") {
    auto n = symbol("n", AssumptionMask{}.set_integer(true));
    auto x = symbol("x");  // generic
    REQUIRE(mod(n, integer(1)) == integer(0));
    // A non-integer (generic) argument keeps Mod(x, 1) (= frac x).
    REQUIRE(mod(x, integer(1))->type_id() == TypeId::Function);
}

// MOD-DIVIDEND-REDUCE-1: the numeric constant of an Add dividend is reduced
// modulo a numeric modulus — Mod(x+5,3) → Mod(x+2,3), Mod(x+2,2) → Mod(x,2)
// (the term drops when c mod q = 0). Non-numeric / symbolic-integer terms are
// kept. Matches SymPy.
TEST_CASE("mod: numeric dividend constant reduces mod q (MOD-DIVIDEND-REDUCE-1)",
          "[3g][mod][oracle][regression]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto n = symbol("n", AssumptionMask{}.set_integer(true));
    REQUIRE(oracle.equivalent(mod(x + integer(5), integer(3))->str(),
                              "Mod(x + 2, 3)"));
    REQUIRE(oracle.equivalent(mod(x + integer(2), integer(2))->str(),
                              "Mod(x, 2)"));
    REQUIRE(oracle.equivalent(mod(x - integer(1), integer(3))->str(),
                              "Mod(x + 2, 3)"));
    REQUIRE(oracle.equivalent(
        mod(integer(2) * x + integer(4), integer(3))->str(), "Mod(2*x + 1, 3)"));
    // Constant already in range → unchanged.
    REQUIRE(mod(x + integer(1), integer(3))->str() == "Mod(x + 1, 3)");
    // Symbolic-integer term is not numeric → not reduced.
    REQUIRE(mod(x + n, integer(1))->type_id() == TypeId::Function);
    REQUIRE(has(mod(x + n, integer(1)), n));
}

TEST_CASE("mod: parse round-trip and SymPy agreement", "[3g][mod][parser][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    REQUIRE(parsing::parse("Mod(x, 3)") == mod(x, integer(3)));
    REQUIRE(mod(x, integer(3))->str() == "Mod(x, 3)");
    REQUIRE(oracle.equivalent(mod(integer(-7), integer(3))->str(), "Mod(-7, 3)"));
    REQUIRE(oracle.equivalent(mod(integer(7), integer(-3))->str(), "Mod(7, -3)"));
}
