// General Jordan form (chains > 2) — verified by reconstruction A = P·J·P⁻¹
// and Jordan-structure checks.

#include <vector>

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/matrices/matrix.hpp>
#include <sympp/simplify/simplify.hpp>

using namespace sympp;

namespace {
// A valid Jordan decomposition: P invertible, A = P·J·P⁻¹, and J is upper
// bidiagonal with super-diagonal entries in {0, 1} and zeros below the diagonal.
void check_jordan(const Matrix& A) {
    auto [P, J] = A.jordan_form();
    std::size_t n = A.rows();
    REQUIRE_FALSE(simplify(P.det()) == S::Zero());
    Matrix recon = P * J * P.inverse();
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            REQUIRE(simplify(recon.at(i, j) - A.at(i, j)) == S::Zero());
            if (j > i + 1 || j < i) {  // strictly below diagonal or above super-diagonal
                REQUIRE(simplify(J.at(i, j)) == S::Zero());
            }
            if (j == i + 1) {
                Expr s = simplify(J.at(i, j));
                REQUIRE((s == S::Zero() || s == integer(1)));
            }
        }
    }
}
}  // namespace

TEST_CASE("Jordan: 3×3 single chain of length 3", "[matrix][jordan]") {
    // Similarity transform of a single λ=2 Jordan block by a unimodular S.
    Matrix Jb = {{integer(2), integer(1), integer(0)},
                 {integer(0), integer(2), integer(1)},
                 {integer(0), integer(0), integer(2)}};
    Matrix S = {{integer(1), integer(1), integer(0)},
                {integer(0), integer(1), integer(1)},
                {integer(0), integer(0), integer(1)}};
    Matrix A = S * Jb * S.inverse();
    check_jordan(A);
}

TEST_CASE("Jordan: 4×4 with two length-2 blocks (same eigenvalue)", "[matrix][jordan]") {
    Matrix A = {{integer(3), integer(1), integer(0), integer(0)},
                {integer(0), integer(3), integer(0), integer(0)},
                {integer(0), integer(0), integer(3), integer(1)},
                {integer(0), integer(0), integer(0), integer(3)}};
    check_jordan(A);
}

TEST_CASE("Jordan: 3×3 with a length-2 and a length-1 block", "[matrix][jordan]") {
    // λ=1, algebraic mult 3, geometric mult 2.
    Matrix A = {{integer(1), integer(1), integer(0)},
                {integer(0), integer(1), integer(0)},
                {integer(0), integer(0), integer(1)}};
    check_jordan(A);
}

TEST_CASE("Jordan: 4×4 length-3 block + length-1 block", "[matrix][jordan]") {
    Matrix A = {{integer(5), integer(1), integer(0), integer(0)},
                {integer(0), integer(5), integer(1), integer(0)},
                {integer(0), integer(0), integer(5), integer(0)},
                {integer(0), integer(0), integer(0), integer(5)}};
    check_jordan(A);
}
