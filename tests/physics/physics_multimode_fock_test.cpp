// Second quantization — multi-mode bosonic ladder operators acting on Fock
// states c·|n₀ n₁ …⟩. Verifies aᵢ|…nᵢ…⟩=√nᵢ|…nᵢ−1…⟩, aᵢ†|…nᵢ…⟩=√(nᵢ+1)|…nᵢ+1…⟩,
// Nᵢ|n⟩=nᵢ|n⟩, the canonical relations [aᵢ,aⱼ†]=δᵢⱼ and [aᵢ,aⱼ]=0, that aᵢ acts
// only on mode i, and vacuum annihilation aᵢ|0⟩=0.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/physics/physics.hpp>
#include <sympp/simplify/simplify.hpp>

using namespace sympp;
namespace ph = sympp::physics;

namespace {
bool coeff_eq(const Expr& a, const Expr& b) { return simplify(a - b) == S::Zero(); }
}  // namespace

TEST_CASE("multimode fock: annihilation lowers only mode i by sqrt(n_i)",
          "[physics][quantum][fock][multimode]") {
    ph::MultiModeFockState s{{3, 5}, S::One()};  // |3 5⟩
    auto out = ph::apply_annihilation(s, 0);
    REQUIRE(out.occupation == std::vector<int>{2, 5});  // mode 1 untouched
    REQUIRE(coeff_eq(out.coefficient, sqrt(integer(3))));

    auto out1 = ph::apply_annihilation(s, 1);
    REQUIRE(out1.occupation == std::vector<int>{3, 4});  // mode 0 untouched
    REQUIRE(coeff_eq(out1.coefficient, sqrt(integer(5))));
}

TEST_CASE("multimode fock: creation raises only mode i by sqrt(n_i+1)",
          "[physics][quantum][fock][multimode]") {
    ph::MultiModeFockState s{{2, 0, 7}, S::One()};  // |2 0 7⟩
    auto out = ph::apply_creation(s, 1);
    REQUIRE(out.occupation == std::vector<int>{2, 1, 7});
    REQUIRE(coeff_eq(out.coefficient, sqrt(integer(1))));
}

TEST_CASE("multimode fock: number operator eigenvalue is n_i",
          "[physics][quantum][fock][multimode]") {
    ph::MultiModeFockState s{{4, 0, 9}, S::One()};
    REQUIRE(coeff_eq(ph::apply_number(s, 0).coefficient, integer(4)));
    REQUIRE(ph::apply_number(s, 1).is_zero());  // n_1 = 0 -> 0
    REQUIRE(coeff_eq(ph::apply_number(s, 2).coefficient, integer(9)));
}

TEST_CASE("multimode fock: a_i|0> = 0 (vacuum annihilation)",
          "[physics][quantum][fock][multimode]") {
    ph::MultiModeFockState vac{{0, 0}, S::One()};
    REQUIRE(ph::apply_annihilation(vac, 0).is_zero());
    REQUIRE(ph::apply_annihilation(vac, 1).is_zero());
}

TEST_CASE("multimode fock: [a_i, a_i^dagger] = 1 (same mode)",
          "[physics][quantum][fock][multimode]") {
    for (const std::vector<int>& occ :
         {std::vector<int>{0, 0}, std::vector<int>{2, 5}, std::vector<int>{9, 1}}) {
        ph::MultiModeFockState s{occ, S::One()};
        auto c0 = ph::apply_boson_commutator(s, 0, 0);
        REQUIRE(c0.occupation == occ);
        REQUIRE(coeff_eq(c0.coefficient, S::One()));
        auto c1 = ph::apply_boson_commutator(s, 1, 1);
        REQUIRE(coeff_eq(c1.coefficient, S::One()));
    }
}

TEST_CASE("multimode fock: [a_1, a_2^dagger] = 0 (different modes commute)",
          "[physics][quantum][fock][multimode]") {
    ph::MultiModeFockState s{{3, 4}, S::One()};  // |3 4⟩
    auto cross = ph::apply_boson_commutator(s, 0, 1);
    REQUIRE(cross.is_zero());
    auto cross2 = ph::apply_boson_commutator(s, 1, 0);
    REQUIRE(cross2.is_zero());
}

TEST_CASE("multimode fock: [a_i, a_j] = 0 (annihilators commute)",
          "[physics][quantum][fock][multimode]") {
    ph::MultiModeFockState s{{3, 4, 5}, S::One()};
    REQUIRE(ph::apply_annihilation_commutator(s, 0, 1).is_zero());  // i != j
    REQUIRE(ph::apply_annihilation_commutator(s, 0, 0).is_zero());  // i == j trivially
    REQUIRE(ph::apply_annihilation_commutator(s, 1, 2).is_zero());
}

TEST_CASE("multimode fock: symbolic coefficient carries through",
          "[physics][quantum][fock][multimode]") {
    auto c = symbol("c");
    ph::MultiModeFockState s{{2, 1}, c};  // c·|2 1⟩
    auto out = ph::apply_creation(s, 0);
    REQUIRE(out.occupation == std::vector<int>{3, 1});
    REQUIRE(coeff_eq(out.coefficient, c * sqrt(integer(3))));
    // [a_0, a_0^dagger] preserves the coefficient.
    auto comm = ph::apply_boson_commutator(s, 0, 0);
    REQUIRE(coeff_eq(comm.coefficient, c));
}
