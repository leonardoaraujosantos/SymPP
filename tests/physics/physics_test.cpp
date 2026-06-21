// Physics core — quantum operator identities, ABCD optics, conjugate momentum.

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/matrices/matrix.hpp>
#include <sympp/physics/physics.hpp>
#include <sympp/simplify/simplify.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
namespace ph = sympp::physics;

namespace {
bool mat_equiv(const Matrix& a, const Matrix& b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) return false;
    for (std::size_t i = 0; i < a.rows(); ++i) {
        for (std::size_t j = 0; j < a.cols(); ++j) {
            if (!(simplify(a.at(i, j) - b.at(i, j)) == S::Zero())) return false;
        }
    }
    return true;
}
}  // namespace

TEST_CASE("quantum: Pauli matrix identities", "[physics][quantum]") {
    auto I2 = Matrix::identity(2);
    REQUIRE(mat_equiv(ph::pauli_x() * ph::pauli_x(), I2));      // σx² = I
    REQUIRE(mat_equiv(ph::pauli_y() * ph::pauli_y(), I2));
    REQUIRE(mat_equiv(ph::pauli_z() * ph::pauli_z(), I2));
    // [σx, σy] = 2i σz
    Matrix expected = ph::pauli_z().scalar_mul(mul(integer(2), S::I()));
    REQUIRE(mat_equiv(ph::commutator(ph::pauli_x(), ph::pauli_y()), expected));
    // {σx, σx} = 2I
    REQUIRE(mat_equiv(ph::anticommutator(ph::pauli_x(), ph::pauli_x()),
                      I2.scalar_mul(integer(2))));
}

TEST_CASE("quantum: ladder and number operators", "[physics][quantum]") {
    auto a = ph::annihilation_operator(4);
    auto ad = ph::creation_operator(4);
    auto N = ph::number_operator(4);
    REQUIRE(mat_equiv(ad, a.transpose()));        // a† = aᵀ (real basis)
    REQUIRE(mat_equiv(ad * a, N));                // a†a = N = diag(0,1,2,3)
    REQUIRE(N.at(2, 2) == integer(2));
    REQUIRE(a.at(0, 1) == integer(1));            // a|1> = √1 |0>
}

TEST_CASE("optics: ABCD ray-transfer matrices", "[physics][optics]") {
    // Free space composes additively.
    REQUIRE(mat_equiv(ph::abcd_free_space(integer(2)) * ph::abcd_free_space(integer(3)),
                      ph::abcd_free_space(integer(5))));
    // A thin lens is unit-determinant.
    auto f = symbol("f");
    REQUIRE(simplify(ph::abcd_thin_lens(f).det()) == integer(1));
}

TEST_CASE("mechanics: conjugate momentum", "[physics][mechanics]") {
    auto m = symbol("m"), v = symbol("v");
    // L = ½ m v²  →  p = ∂L/∂v = m v.
    Expr L = mul(rational(1, 2), mul(m, pow(v, integer(2))));
    REQUIRE(sympp::testing::Oracle::instance().equivalent(
        ph::conjugate_momentum(L, v)->str(), (m * v)->str()));
}
