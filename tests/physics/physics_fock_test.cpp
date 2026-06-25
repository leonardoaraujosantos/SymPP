// Second quantization — single-mode bosonic ladder operators acting on Fock
// states c·|n⟩. Verifies the canonical relations a|n⟩=√n|n−1⟩, a†|n⟩=√(n+1)|n+1⟩,
// N|n⟩=n|n⟩, a†a|n⟩=n|n⟩, aa†|n⟩=(n+1)|n⟩, [a,a†]|n⟩=|n⟩ and a|0⟩=0.

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

TEST_CASE("fock: annihilation lowers and scales by sqrt(n)", "[physics][quantum][fock]") {
    ph::FockState s{3, S::One()};  // |3⟩
    auto out = ph::apply_annihilation(s);
    REQUIRE(out.occupation == 2);
    REQUIRE(coeff_eq(out.coefficient, sqrt(integer(3))));
}

TEST_CASE("fock: creation raises and scales by sqrt(n+1)", "[physics][quantum][fock]") {
    ph::FockState s{3, S::One()};  // |3⟩
    auto out = ph::apply_creation(s);
    REQUIRE(out.occupation == 4);
    REQUIRE(coeff_eq(out.coefficient, sqrt(integer(4))));
}

TEST_CASE("fock: number operator eigenvalue is n", "[physics][quantum][fock]") {
    for (std::size_t n : {0u, 1u, 5u}) {
        ph::FockState s{n, S::One()};
        auto out = ph::apply_number(s);
        REQUIRE(out.occupation == n);
        REQUIRE(coeff_eq(out.coefficient, integer(static_cast<long>(n))));
    }
}

TEST_CASE("fock: a|0> = 0 (vacuum annihilation)", "[physics][quantum][fock]") {
    ph::FockState vac{0, S::One()};
    auto out = ph::apply_annihilation(vac);
    REQUIRE(out.is_zero());
}

TEST_CASE("fock: a-dagger a |n> = n |n>", "[physics][quantum][fock]") {
    for (std::size_t n : {0u, 1u, 4u, 7u}) {
        ph::FockState s{n, S::One()};
        auto out = ph::apply_creation(ph::apply_annihilation(s));
        if (n == 0) {
            REQUIRE(out.is_zero());  // a†a|0⟩ = 0 = 0·|0⟩
        } else {
            REQUIRE(out.occupation == n);
            REQUIRE(coeff_eq(out.coefficient, integer(static_cast<long>(n))));
        }
    }
}

TEST_CASE("fock: a a-dagger |n> = (n+1) |n>", "[physics][quantum][fock]") {
    for (std::size_t n : {0u, 1u, 4u, 7u}) {
        ph::FockState s{n, S::One()};
        auto out = ph::apply_annihilation(ph::apply_creation(s));
        REQUIRE(out.occupation == n);
        REQUIRE(coeff_eq(out.coefficient, integer(static_cast<long>(n + 1))));
    }
}

TEST_CASE("fock: [a, a-dagger] |n> = |n> (canonical commutation)", "[physics][quantum][fock]") {
    for (std::size_t n : {0u, 1u, 2u, 9u}) {
        ph::FockState s{n, S::One()};
        auto out = ph::apply_boson_commutator(s);
        REQUIRE(out.occupation == n);
        REQUIRE(coeff_eq(out.coefficient, S::One()));  // eigenvalue 1
    }
}

TEST_CASE("fock: coefficient carries through ladder operators", "[physics][quantum][fock]") {
    auto c = symbol("c");
    ph::FockState s{2, c};  // c·|2⟩
    auto out = ph::apply_creation(s);
    REQUIRE(out.occupation == 3);
    REQUIRE(coeff_eq(out.coefficient, c * sqrt(integer(3))));
    // [a,a†] preserves the symbolic coefficient: [a,a†](c·|2⟩) = c·|2⟩.
    auto comm = ph::apply_boson_commutator(s);
    REQUIRE(comm.occupation == 2);
    REQUIRE(coeff_eq(comm.coefficient, c));
}
