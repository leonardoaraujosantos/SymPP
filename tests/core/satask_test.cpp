// Boolean / SAT-style `ask` over predicate propositions.
//
// Reference behaviour from sympy/assumptions/satask.py — `ask` over boolean
// combinations of predicates (Q.positive(x) | Q.negative(x), etc.).

#include <optional>

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/satask.hpp>
#include <sympp/core/symbol.hpp>

using namespace sympp;

namespace {
constexpr auto Pos = AssumptionKey::Positive;
constexpr auto Neg = AssumptionKey::Negative;
constexpr auto Zer = AssumptionKey::Zero;
constexpr auto Even = AssumptionKey::Even;
constexpr auto Odd = AssumptionKey::Odd;
constexpr auto Real = AssumptionKey::Real;
constexpr auto Nonneg = AssumptionKey::Nonnegative;
constexpr auto Nonpos = AssumptionKey::Nonpositive;
}  // namespace

TEST_CASE("satask: Kleene evaluation of decidable combinations", "[satask]") {
    auto p = symbol("p", AssumptionMask{}.set_positive(true));
    auto n = symbol("n", AssumptionMask{}.set_negative(true));

    // Both decidable directly.
    REQUIRE(ask(Q(p, Pos) || Q(n, Pos)) == true);     // p>0 true
    REQUIRE(ask(Q(p, Pos) && Q(n, Neg)) == true);     // p>0 ∧ n<0
    REQUIRE(ask(Q(p, Neg) || Q(n, Pos)) == false);    // both false
    REQUIRE(ask(!Q(p, Neg)) == true);                  // ¬(p<0)
    REQUIRE(ask(Q(p, Pos) && Q(n, Pos)) == false);    // n>0 is false
}

TEST_CASE("satask: exhaustive disjunction over a single expression", "[satask]") {
    // A real, nonzero symbol is positive ∨ negative even though neither
    // disjunct is individually decidable — SymPy returns True here.
    auto x = symbol("x", AssumptionMask{}.set_real(true).set_nonzero(true));
    REQUIRE(!ask(x, Pos).has_value());  // individually unknown
    REQUIRE(!ask(x, Neg).has_value());
    REQUIRE(ask(Q(x, Pos) || Q(x, Neg)) == true);

    // A real symbol is positive ∨ negative ∨ zero (trichotomy).
    auto y = symbol("y", AssumptionMask{}.set_real(true));
    REQUIRE(ask(Q(y, Pos) || Q(y, Neg) || Q(y, Zer)) == true);
    // ... and ≥0 ∨ ≤0.
    REQUIRE(ask(Q(y, Nonneg) || Q(y, Nonpos)) == true);
    // But without nonzero, positive ∨ negative is NOT entailed.
    REQUIRE(!ask(Q(y, Pos) || Q(y, Neg)).has_value());

    // An integer is even ∨ odd.
    auto k = symbol("k", AssumptionMask{}.set_integer(true));
    REQUIRE(ask(Q(k, Even) || Q(k, Odd)) == true);
}

TEST_CASE("satask: contradictions are refuted to false", "[satask]") {
    auto x = symbol("x", AssumptionMask{}.set_real(true).set_nonzero(true));
    // Can't be both positive and negative.
    REQUIRE(ask(Q(x, Pos) && Q(x, Neg)) == false);
    // A nonzero real is not zero.
    REQUIRE(ask(Q(x, Zer)) == false);

    auto k = symbol("k", AssumptionMask{}.set_integer(true));
    REQUIRE(ask(Q(k, Even) && Q(k, Odd)) == false);
}

TEST_CASE("satask: concrete numbers", "[satask]") {
    auto two = integer(2);
    REQUIRE(ask(Q(two, Pos) || Q(two, Neg)) == true);
    REQUIRE(ask(Q(two, Even) || Q(two, Odd)) == true);
    REQUIRE(ask(Q(two, Odd)) == false);
}

TEST_CASE("satask: parity with SymPy's definite answers", "[satask]") {
    // Mirrors a battery of combinations where SymPy's `ask` is definite; SymPP
    // must agree (and never contradict a definite SymPy answer).
    auto p = symbol("p", AssumptionMask{}.set_positive(true));
    auto i = symbol("i", AssumptionMask{}.set_integer(true));

    REQUIRE(ask(Q(p, Real) && Q(p, AssumptionKey::Finite)) == true);   // True
    REQUIRE(ask(!Q(p, Neg)) == true);                                  // True
    REQUIRE(ask(!Q(i, AssumptionKey::Irrational)) == true);            // True
}

TEST_CASE("satask: generic symbol stays unknown", "[satask]") {
    auto x = symbol("x");  // generic complex
    REQUIRE(!ask(Q(x, Pos) || Q(x, Neg)).has_value());
    REQUIRE(!ask(Q(x, Real)).has_value());
}

TEST_CASE("assumptions_consistent: genuine contradictions", "[satask]") {
    REQUIRE(assumptions_consistent(AssumptionMask{}));  // empty is consistent
    REQUIRE(assumptions_consistent(AssumptionMask{}.set_real(true)));

    AssumptionMask realNotSigned;  // real ∧ ¬pos ∧ ¬neg ∧ ¬zero — impossible
    realNotSigned.set_real(true);
    realNotSigned.set(Pos, false);
    realNotSigned.set(Neg, false);
    realNotSigned.set(Zer, false);
    REQUIRE_FALSE(assumptions_consistent(realNotSigned));

    AssumptionMask posAndNeg;
    posAndNeg.set(Pos, true);
    posAndNeg.set(Neg, true);
    REQUIRE_FALSE(assumptions_consistent(posAndNeg));
}

TEST_CASE("commutative: default true; non-commutative is structural", "[commutative]") {
    // Ordinary symbols and numbers are commutative by default (as in SymPy).
    REQUIRE(is_commutative(symbol("x")) == true);
    REQUIRE(is_commutative(integer(7)) == true);

    // A symbol declared commutative=false is non-commutative ...
    auto a = symbol("a", AssumptionMask{}.set_commutative(false));
    auto b = symbol("b", AssumptionMask{}.set_commutative(false));
    auto c = symbol("c");  // commutative
    REQUIRE(is_commutative(a) == false);

    // ... and that propagates through every compound it appears in.
    REQUIRE(is_commutative(mul(a, b)) == false);
    REQUIRE(is_commutative(mul(c, a)) == false);
    REQUIRE(is_commutative(c + a) == false);
    REQUIRE(is_commutative(pow(a, integer(2))) == false);

    // An expression built only from commutative pieces stays commutative.
    REQUIRE(is_commutative(mul(c, integer(3))) == true);
    REQUIRE(is_commutative(pow(c, integer(2)) + integer(1)) == true);
}
