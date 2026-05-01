// Matrix tests.
//
// References:
//   sympy/matrices/tests/test_commonmatrix.py — det/inv/transpose/trace
//   sympy/matrices/tests/test_matrices.py
//   MATLAB Symbolic Math Toolbox §3 — Hilbert matrix, jacobian

#include <stdexcept>

#include <catch2/catch_test_macros.hpp>

#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/matrices/matrix.hpp>

#include "oracle/oracle.hpp"

using namespace sympp;
using sympp::testing::Oracle;

// ----- Construction ----------------------------------------------------------

TEST_CASE("Matrix: zeros and identity", "[9a][matrix]") {
    auto z = Matrix::zeros(2, 3);
    REQUIRE(z.rows() == 2);
    REQUIRE(z.cols() == 3);
    REQUIRE(z.at(0, 0) == S::Zero());
    REQUIRE(z.at(1, 2) == S::Zero());

    auto i = Matrix::identity(3);
    REQUIRE(i.at(0, 0) == S::One());
    REQUIRE(i.at(1, 1) == S::One());
    REQUIRE(i.at(2, 2) == S::One());
    REQUIRE(i.at(0, 1) == S::Zero());
}

TEST_CASE("Matrix: initializer-list construction", "[9a][matrix]") {
    Matrix m = {{integer(1), integer(2)}, {integer(3), integer(4)}};
    REQUIRE(m.rows() == 2);
    REQUIRE(m.cols() == 2);
    REQUIRE(m.at(0, 0) == integer(1));
    REQUIRE(m.at(1, 1) == integer(4));
}

TEST_CASE("Matrix: ragged initializer throws", "[9a][matrix]") {
    REQUIRE_THROWS_AS((Matrix{{integer(1), integer(2)}, {integer(3)}}),
                      std::invalid_argument);
}

// ----- Arithmetic ------------------------------------------------------------

TEST_CASE("Matrix: add", "[9a][matrix]") {
    Matrix a = {{integer(1), integer(2)}, {integer(3), integer(4)}};
    Matrix b = {{integer(5), integer(6)}, {integer(7), integer(8)}};
    Matrix c = a + b;
    REQUIRE(c.at(0, 0) == integer(6));
    REQUIRE(c.at(1, 1) == integer(12));
}

TEST_CASE("Matrix: multiply", "[9a][matrix]") {
    // [[1,2],[3,4]] * [[5,6],[7,8]] = [[19,22],[43,50]]
    Matrix a = {{integer(1), integer(2)}, {integer(3), integer(4)}};
    Matrix b = {{integer(5), integer(6)}, {integer(7), integer(8)}};
    Matrix c = a * b;
    REQUIRE(c.at(0, 0) == integer(19));
    REQUIRE(c.at(0, 1) == integer(22));
    REQUIRE(c.at(1, 0) == integer(43));
    REQUIRE(c.at(1, 1) == integer(50));
}

TEST_CASE("Matrix: multiply identity", "[9a][matrix]") {
    Matrix a = {{integer(1), integer(2)}, {integer(3), integer(4)}};
    Matrix i = Matrix::identity(2);
    REQUIRE((a * i).equals(a));
    REQUIRE((i * a).equals(a));
}

TEST_CASE("Matrix: scalar mul", "[9a][matrix]") {
    Matrix a = {{integer(1), integer(2)}};
    Matrix b = a.scalar_mul(integer(3));
    REQUIRE(b.at(0, 0) == integer(3));
    REQUIRE(b.at(0, 1) == integer(6));
}

// ----- Transpose / trace -----------------------------------------------------

TEST_CASE("Matrix: transpose", "[9a][matrix]") {
    Matrix a = {{integer(1), integer(2), integer(3)},
                {integer(4), integer(5), integer(6)}};
    Matrix t = a.transpose();
    REQUIRE(t.rows() == 3);
    REQUIRE(t.cols() == 2);
    REQUIRE(t.at(0, 0) == integer(1));
    REQUIRE(t.at(2, 1) == integer(6));
}

TEST_CASE("Matrix: trace", "[9a][matrix]") {
    Matrix a = {{integer(1), integer(2)}, {integer(3), integer(4)}};
    REQUIRE(a.trace() == integer(5));
    REQUIRE(Matrix::identity(4).trace() == integer(4));
}

// ----- Determinant -----------------------------------------------------------

TEST_CASE("Matrix: 1x1 det", "[9a][det]") {
    Matrix a = {{integer(7)}};
    REQUIRE(a.det() == integer(7));
}

TEST_CASE("Matrix: 2x2 det", "[9a][det]") {
    // [[1,2],[3,4]] -> 1*4 - 2*3 = -2
    Matrix a = {{integer(1), integer(2)}, {integer(3), integer(4)}};
    REQUIRE(a.det() == integer(-2));
}

TEST_CASE("Matrix: 3x3 det matches SymPy", "[9a][det][oracle]") {
    auto& oracle = Oracle::instance();
    // [[1,2,3],[4,5,6],[7,8,10]] -> det = -3
    Matrix a = {{integer(1), integer(2), integer(3)},
                {integer(4), integer(5), integer(6)},
                {integer(7), integer(8), integer(10)}};
    auto d = a.det();
    REQUIRE(oracle.equivalent(d->str(), "-3"));
}

TEST_CASE("Matrix: identity has det 1", "[9a][det]") {
    REQUIRE(Matrix::identity(3).det() == integer(1));
    REQUIRE(Matrix::identity(5).det() == integer(1));
}

// ----- Inverse ---------------------------------------------------------------

TEST_CASE("Matrix: 2x2 inverse", "[9a][inverse]") {
    Matrix a = {{integer(1), integer(2)}, {integer(3), integer(4)}};
    Matrix inv = a.inverse();
    Matrix prod = a * inv;
    Matrix i = Matrix::identity(2);
    REQUIRE(prod.equals(i));
}

TEST_CASE("Matrix: singular matrix throws", "[9a][inverse]") {
    Matrix a = {{integer(1), integer(2)}, {integer(2), integer(4)}};
    REQUIRE_THROWS_AS(a.inverse(), std::invalid_argument);
}

// ----- Symbolic --------------------------------------------------------------

TEST_CASE("Matrix: symbolic det", "[9a][matrix][oracle]") {
    auto& oracle = Oracle::instance();
    auto a = symbol("a");
    auto b = symbol("b");
    auto c = symbol("c");
    auto d = symbol("d");
    Matrix m = {{a, b}, {c, d}};
    auto det_expr = m.det();
    REQUIRE(oracle.equivalent(det_expr->str(), "a*d - b*c"));
}

// ----- Jacobian --------------------------------------------------------------

TEST_CASE("Matrix: jacobian on simple polynomials", "[9a][jacobian]") {
    auto x = symbol("x");
    auto y = symbol("y");
    // F = [x^2*y, sin(... no — let's keep elementary)]
    // F = [x*y, x + y]
    auto j = jacobian({x * y, x + y}, {x, y});
    REQUIRE(j.rows() == 2);
    REQUIRE(j.cols() == 2);
    REQUIRE(j.at(0, 0) == y);
    REQUIRE(j.at(0, 1) == x);
    REQUIRE(j.at(1, 0) == integer(1));
    REQUIRE(j.at(1, 1) == integer(1));
}

TEST_CASE("Matrix: jacobian matches SymPy", "[9a][jacobian][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto y = symbol("y");
    auto j = jacobian({pow(x, integer(2)) + y, integer(2) * x * y}, {x, y});
    REQUIRE(oracle.equivalent(j.at(0, 0)->str(), "2*x"));
    REQUIRE(oracle.equivalent(j.at(0, 1)->str(), "1"));
    REQUIRE(oracle.equivalent(j.at(1, 0)->str(), "2*y"));
    REQUIRE(oracle.equivalent(j.at(1, 1)->str(), "2*x"));
}
