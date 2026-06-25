// Second quantization — single-mode fermionic ladder operators acting on Fock
// states c·|n⟩ with Pauli exclusion (n ∈ {0,1}). Verifies a|1⟩=|0⟩, a|0⟩=0,
// a†|0⟩=|1⟩, a†|1⟩=0, N|n⟩=n|n⟩, a†a|n⟩=n|n⟩ and the canonical anticommutation
// relation {a,a†}|n⟩=|n⟩.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/physics/physics.hpp>
#include <sympp/simplify/simplify.hpp>

using namespace sympp;
namespace ph = sympp::physics;

namespace {
bool coeff_eq(const Expr& a, const Expr& b) { return simplify(a - b) == S::Zero(); }
}  // namespace

TEST_CASE("fermion: a|1> = |0>", "[physics][quantum][fermion]") {
    ph::FermionState s{1, S::One()};
    auto out = ph::apply_fermion_annihilation(s);
    REQUIRE(out.occupation == 0);
    REQUIRE(coeff_eq(out.coefficient, S::One()));
}

TEST_CASE("fermion: a|0> = 0", "[physics][quantum][fermion]") {
    ph::FermionState vac{0, S::One()};
    auto out = ph::apply_fermion_annihilation(vac);
    REQUIRE(out.is_zero());
}

TEST_CASE("fermion: a-dagger|0> = |1>", "[physics][quantum][fermion]") {
    ph::FermionState vac{0, S::One()};
    auto out = ph::apply_fermion_creation(vac);
    REQUIRE(out.occupation == 1);
    REQUIRE(coeff_eq(out.coefficient, S::One()));
}

TEST_CASE("fermion: a-dagger|1> = 0 (Pauli exclusion)", "[physics][quantum][fermion]") {
    ph::FermionState s{1, S::One()};
    auto out = ph::apply_fermion_creation(s);
    REQUIRE(out.is_zero());
}

TEST_CASE("fermion: number operator eigenvalue is n", "[physics][quantum][fermion]") {
    for (std::size_t n : {0u, 1u}) {
        ph::FermionState s{n, S::One()};
        auto out = ph::apply_fermion_number(s);
        REQUIRE(out.occupation == n);
        REQUIRE(coeff_eq(out.coefficient, integer(static_cast<long>(n))));
    }
}

TEST_CASE("fermion: a-dagger a |1> = |1>", "[physics][quantum][fermion]") {
    ph::FermionState s{1, S::One()};
    auto out = ph::apply_fermion_creation(ph::apply_fermion_annihilation(s));
    REQUIRE(out.occupation == 1);
    REQUIRE(coeff_eq(out.coefficient, S::One()));
}

TEST_CASE("fermion: a-dagger a |0> = 0", "[physics][quantum][fermion]") {
    ph::FermionState s{0, S::One()};
    auto out = ph::apply_fermion_creation(ph::apply_fermion_annihilation(s));
    REQUIRE(out.is_zero());
}

TEST_CASE("fermion: a a-dagger |0> = |0>", "[physics][quantum][fermion]") {
    ph::FermionState s{0, S::One()};
    auto out = ph::apply_fermion_annihilation(ph::apply_fermion_creation(s));
    REQUIRE(out.occupation == 0);
    REQUIRE(coeff_eq(out.coefficient, S::One()));
}

TEST_CASE("fermion: {a, a-dagger} |n> = |n> (canonical anticommutation)",
          "[physics][quantum][fermion]") {
    for (std::size_t n : {0u, 1u}) {
        ph::FermionState s{n, S::One()};
        auto out = ph::apply_fermion_anticommutator(s);
        REQUIRE(out.occupation == n);
        REQUIRE(coeff_eq(out.coefficient, S::One()));  // eigenvalue 1
    }
}

TEST_CASE("fermion: coefficient carries through ladder operators",
          "[physics][quantum][fermion]") {
    auto c = symbol("c");
    ph::FermionState s0{0, c};  // c·|0⟩
    auto up = ph::apply_fermion_creation(s0);
    REQUIRE(up.occupation == 1);
    REQUIRE(coeff_eq(up.coefficient, c));

    ph::FermionState s1{1, c};  // c·|1⟩
    auto down = ph::apply_fermion_annihilation(s1);
    REQUIRE(down.occupation == 0);
    REQUIRE(coeff_eq(down.coefficient, c));

    // {a,a†} preserves the symbolic coefficient: {a,a†}(c·|n⟩) = c·|n⟩.
    auto anti0 = ph::apply_fermion_anticommutator(s0);
    REQUIRE(anti0.occupation == 0);
    REQUIRE(coeff_eq(anti0.coefficient, c));
    auto anti1 = ph::apply_fermion_anticommutator(s1);
    REQUIRE(anti1.occupation == 1);
    REQUIRE(coeff_eq(anti1.coefficient, c));
}
