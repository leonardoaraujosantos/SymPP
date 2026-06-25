// Second quantization — normal ordering of bosonic / fermionic operator words.
// Verifies the canonical (anti)commutation rewrites
//   boson:    a a† = a† a + 1
//   fermion:  a a† = 1 − a† a,   a† a† = 0
//   cross:    a_i a_j† = a_j† a_i        (i ≠ j)
// and cross-checks the boson expansions a a a† a† / a a† a† a numerically by
// evaluating the original word and its normal-ordered form on a Fock state |n⟩.

#include <catch2/catch_test_macros.hpp>

#include <vector>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/physics/physics.hpp>
#include <sympp/simplify/simplify.hpp>

using namespace sympp;
namespace ph = sympp::physics;

namespace {
bool coeff_eq(const Expr& a, const Expr& b) { return simplify(a - b) == S::Zero(); }

ph::LadderOp ann(std::size_t m = 0) { return ph::LadderOp{m, false}; }
ph::LadderOp cre(std::size_t m = 0) { return ph::LadderOp{m, true}; }

// Find the coefficient attached to `target` in a normal-ordered result; zero if
// the word is absent.
Expr coeff_of(const std::vector<ph::OperatorTerm>& terms, const ph::OperatorWord& target) {
    for (const auto& t : terms) {
        if (t.word == target) return t.coefficient;
    }
    return S::Zero();
}

// Apply a single ladder operator to a single-mode bosonic Fock state.
ph::FockState apply_one(const ph::LadderOp& op, const ph::FockState& s) {
    return op.dagger ? ph::apply_creation(s) : ph::apply_annihilation(s);
}

// Evaluate a single-mode operator word on |n⟩ (right-to-left application).
ph::FockState eval_word(const ph::OperatorWord& w, std::size_t n) {
    ph::FockState s{n, S::One()};
    for (auto it = w.rbegin(); it != w.rend(); ++it) s = apply_one(*it, s);
    return s;
}

// Evaluate a normal-ordered linear combination on |n⟩, summing terms that land
// back on occupation n (all our test words are number-conserving).
Expr eval_terms(const std::vector<ph::OperatorTerm>& terms, std::size_t n) {
    Expr total = S::Zero();
    for (const auto& t : terms) {
        ph::FockState s = eval_word(t.word, n);
        if (s.is_zero()) continue;
        total = total + t.coefficient * (s.occupation == n ? s.coefficient : S::Zero());
    }
    return simplify(total);
}
}  // namespace

TEST_CASE("normal_order: boson a a† = a† a + 1", "[physics][quantum][normal_order]") {
    auto terms = ph::normal_order({ann(), cre()}, ph::Statistics::Boson);
    REQUIRE(coeff_eq(coeff_of(terms, {cre(), ann()}), S::One()));  // a† a term
    REQUIRE(coeff_eq(coeff_of(terms, {}), S::One()));              // +1 (identity)
    REQUIRE(terms.size() == 2);
}

TEST_CASE("normal_order: fermion a a† = 1 − a† a", "[physics][quantum][normal_order]") {
    auto terms = ph::normal_order({ann(), cre()}, ph::Statistics::Fermion);
    REQUIRE(coeff_eq(coeff_of(terms, {cre(), ann()}), integer(-1)));  // −a† a
    REQUIRE(coeff_eq(coeff_of(terms, {}), S::One()));                 // +1
    REQUIRE(terms.size() == 2);
}

TEST_CASE("normal_order: fermion a† a† = 0", "[physics][quantum][normal_order]") {
    auto terms = ph::normal_order({cre(), cre()}, ph::Statistics::Fermion);
    REQUIRE(terms.empty());
}

TEST_CASE("normal_order: fermion a a = 0", "[physics][quantum][normal_order]") {
    auto terms = ph::normal_order({ann(), ann()}, ph::Statistics::Fermion);
    REQUIRE(terms.empty());
}

TEST_CASE("normal_order: cross-mode a_i a_j† = a_j† a_i (i != j)",
          "[physics][quantum][normal_order]") {
    // Bosons commute across modes: a_0 a_1† = a_1† a_0 (coeff +1).
    auto bosons = ph::normal_order({ann(0), cre(1)}, ph::Statistics::Boson);
    REQUIRE(bosons.size() == 1);
    REQUIRE(coeff_eq(coeff_of(bosons, {cre(1), ann(0)}), S::One()));
    // Fermions anticommute across modes: a_0 a_1† = −a_1† a_0 (coeff −1).
    auto fermions = ph::normal_order({ann(0), cre(1)}, ph::Statistics::Fermion);
    REQUIRE(fermions.size() == 1);
    REQUIRE(coeff_eq(coeff_of(fermions, {cre(1), ann(0)}), integer(-1)));
}

TEST_CASE("normal_order: already-ordered word is unchanged",
          "[physics][quantum][normal_order]") {
    ph::OperatorWord w{cre(), cre(), ann(), ann()};  // a† a† a a
    auto terms = ph::normal_order(w, ph::Statistics::Boson);
    REQUIRE(terms.size() == 1);
    REQUIRE(coeff_eq(coeff_of(terms, w), S::One()));
}

TEST_CASE("normal_order: boson a a a† a† expands and matches |n>",
          "[physics][quantum][normal_order]") {
    // a a a† a† = a† a† a a + 4 a† a + 2  (standard expansion).
    ph::OperatorWord w{ann(), ann(), cre(), cre()};
    auto terms = ph::normal_order(w, ph::Statistics::Boson);

    REQUIRE(coeff_eq(coeff_of(terms, {cre(), cre(), ann(), ann()}), S::One()));
    REQUIRE(coeff_eq(coeff_of(terms, {cre(), ann()}), integer(4)));
    REQUIRE(coeff_eq(coeff_of(terms, {}), integer(2)));

    // Numerical cross-check on |n⟩ for several n: original vs normal-ordered.
    for (std::size_t n : {0u, 1u, 2u, 3u, 4u}) {
        Expr lhs = eval_terms({ph::OperatorTerm{S::One(), w}}, n);
        Expr rhs = eval_terms(terms, n);
        REQUIRE(coeff_eq(lhs, rhs));
    }
}

TEST_CASE("normal_order: boson a a† a† a matches |n>",
          "[physics][quantum][normal_order]") {
    ph::OperatorWord w{ann(), cre(), cre(), ann()};
    auto terms = ph::normal_order(w, ph::Statistics::Boson);
    for (std::size_t n : {0u, 1u, 2u, 3u, 4u}) {
        Expr lhs = eval_terms({ph::OperatorTerm{S::One(), w}}, n);
        Expr rhs = eval_terms(terms, n);
        REQUIRE(coeff_eq(lhs, rhs));
    }
}
