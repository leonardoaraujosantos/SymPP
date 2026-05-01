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
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/matrices/matrix.hpp>
#include <sympp/simplify/simplify.hpp>

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

// ----- Standard constructors -------------------------------------------------

TEST_CASE("matrix: hilbert(3) has 1/(i+j+1) entries", "[9][matrix][hilbert]") {
    auto h = hilbert(3);
    REQUIRE(h.at(0, 0) == integer(1));
    REQUIRE(h.at(0, 1) == rational(1, 2));
    REQUIRE(h.at(1, 1) == rational(1, 3));
    REQUIRE(h.at(2, 2) == rational(1, 5));
}

TEST_CASE("matrix: vandermonde with [a, b, c] cols=3",
          "[9][matrix][vandermonde]") {
    auto a = symbol("a");
    auto b = symbol("b");
    auto c = symbol("c");
    auto v = vandermonde({a, b, c}, 3);
    REQUIRE(v.rows() == 3);
    REQUIRE(v.cols() == 3);
    REQUIRE(v.at(0, 0) == integer(1));
    REQUIRE(v.at(1, 1) == b);
    REQUIRE(v.at(2, 2) == pow(c, integer(2)));
}

TEST_CASE("matrix: companion of x^2 + 2x + 3 has 0,-3 / 1,-2",
          "[9][matrix][companion][oracle]") {
    auto& oracle = Oracle::instance();
    auto comp = companion({integer(3), integer(2)});  // x^2 + 2x + 3
    REQUIRE(comp.rows() == 2);
    REQUIRE(comp.at(0, 0) == S::Zero());
    REQUIRE(comp.at(1, 0) == S::One());
    REQUIRE(comp.at(0, 1) == integer(-3));
    REQUIRE(comp.at(1, 1) == integer(-2));
    // Characteristic polynomial: det(xI - C) = x^2 + 2x + 3.
    auto x = symbol("x");
    Matrix lambda_I = Matrix::identity(2).scalar_mul(x);
    Matrix delta = lambda_I - comp;
    auto cp = delta.det();
    REQUIRE(oracle.equivalent(cp->str(), "x**2 + 2*x + 3"));
}

TEST_CASE("matrix: rotation_matrix_2d", "[9][matrix][rotation][oracle]") {
    auto& oracle = Oracle::instance();
    auto theta = symbol("theta");
    auto r = rotation_matrix_2d(theta);
    // det should be cos²θ + sin²θ = 1.
    auto d = r.det();
    REQUIRE(oracle.equivalent(d->str(), "1"));
}

TEST_CASE("matrix: rotation_matrix_z preserves z-axis",
          "[9][matrix][rotation]") {
    auto theta = symbol("theta");
    auto r = rotation_matrix_z(theta);
    REQUIRE(r.at(2, 2) == integer(1));
    REQUIRE(r.at(0, 2) == S::Zero());
    REQUIRE(r.at(2, 0) == S::Zero());
}

// ----- Hadamard / Kronecker -------------------------------------------------

TEST_CASE("matrix: hadamard product element-wise",
          "[9][matrix][hadamard]") {
    Matrix a = {{integer(1), integer(2)}, {integer(3), integer(4)}};
    Matrix b = {{integer(5), integer(6)}, {integer(7), integer(8)}};
    auto h = hadamard_product(a, b);
    REQUIRE(h.at(0, 0) == integer(5));
    REQUIRE(h.at(0, 1) == integer(12));
    REQUIRE(h.at(1, 0) == integer(21));
    REQUIRE(h.at(1, 1) == integer(32));
}

TEST_CASE("matrix: kronecker product 2x2 ⊗ 2x2 → 4x4",
          "[9][matrix][kronecker]") {
    Matrix a = {{integer(1), integer(2)}, {integer(0), integer(1)}};
    Matrix b = {{integer(0), integer(1)}, {integer(1), integer(0)}};
    auto k = kronecker_product(a, b);
    REQUIRE(k.rows() == 4);
    REQUIRE(k.cols() == 4);
    // Top-left block: 1 * b
    REQUIRE(k.at(0, 0) == S::Zero());
    REQUIRE(k.at(0, 1) == integer(1));
    REQUIRE(k.at(1, 0) == integer(1));
    // Top-right block: 2 * b
    REQUIRE(k.at(0, 2) == S::Zero());
    REQUIRE(k.at(0, 3) == integer(2));
    REQUIRE(k.at(1, 2) == integer(2));
    // Bottom-right block: 1 * b (since a[1,1] = 1)
    REQUIRE(k.at(2, 2) == S::Zero());
    REQUIRE(k.at(2, 3) == integer(1));
    REQUIRE(k.at(3, 2) == integer(1));
}

#include <sympp/calculus/diff.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/trigonometric.hpp>

// ----- adjugate / conjugate_transpose / norm --------------------------------

TEST_CASE("matrix: adjugate of 2x2", "[9][matrix][adjugate][oracle]") {
    Matrix m = {{integer(1), integer(2)}, {integer(3), integer(4)}};
    auto adj = m.adjugate();
    REQUIRE(adj.at(0, 0) == integer(4));
    REQUIRE(adj.at(0, 1) == integer(-2));
    REQUIRE(adj.at(1, 0) == integer(-3));
    REQUIRE(adj.at(1, 1) == integer(1));
}

TEST_CASE("matrix: A·adj(A) = det(A)·I", "[9][matrix][adjugate]") {
    Matrix m = {{integer(2), integer(1), integer(3)},
                {integer(0), integer(4), integer(1)},
                {integer(5), integer(2), integer(0)}};
    auto adj = m.adjugate();
    auto prod = m * adj;
    auto d = m.det();
    auto expected = Matrix::identity(3).scalar_mul(d);
    REQUIRE(prod.equals(expected));
}

TEST_CASE("matrix: conjugate_transpose of real == transpose",
          "[9][matrix][conj_transpose]") {
    Matrix m = {{integer(1), integer(2)}, {integer(3), integer(4)}};
    auto ct = m.conjugate_transpose();
    auto t = m.transpose();
    REQUIRE(ct.equals(t));
}

TEST_CASE("matrix: norm_frobenius of [[1,2],[3,4]] = sqrt(30)",
          "[9][matrix][norm][oracle]") {
    auto& oracle = Oracle::instance();
    Matrix m = {{integer(1), integer(2)}, {integer(3), integer(4)}};
    auto n = m.norm_frobenius();
    REQUIRE(oracle.equivalent(n->str(), "sqrt(30)"));
}

// ----- gradient / hessian / wronskian --------------------------------------

TEST_CASE("matrix: gradient of x²+xy is column [2x+y, x]",
          "[9][matrix][gradient][oracle]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto f = pow(x, integer(2)) + x * y;
    auto g = gradient(f, {x, y});
    auto& oracle = Oracle::instance();
    REQUIRE(g.rows() == 2);
    REQUIRE(g.cols() == 1);
    REQUIRE(oracle.equivalent(g.at(0, 0)->str(), "2*x + y"));
    REQUIRE(g.at(1, 0) == x);
}

TEST_CASE("matrix: hessian of x²+y² is 2I",
          "[9][matrix][hessian][oracle]") {
    auto x = symbol("x");
    auto y = symbol("y");
    auto f = pow(x, integer(2)) + pow(y, integer(2));
    auto H = hessian(f, {x, y});
    REQUIRE(H.at(0, 0) == integer(2));
    REQUIRE(H.at(1, 1) == integer(2));
    REQUIRE(H.at(0, 1) == S::Zero());
    REQUIRE(H.at(1, 0) == S::Zero());
}

TEST_CASE("matrix: wronskian of {1, x, x²} = 2",
          "[9][matrix][wronskian][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto w = wronskian({integer(1), x, pow(x, integer(2))}, x);
    REQUIRE(oracle.equivalent(w->str(), "2"));
}

TEST_CASE("matrix: wronskian of {sin(x), cos(x)} = -1",
          "[9][matrix][wronskian][oracle]") {
    auto& oracle = Oracle::instance();
    auto x = symbol("x");
    auto w = wronskian({sin(x), cos(x)}, x);
    REQUIRE(oracle.equivalent(w->str(), "-1"));
}

// ----- rref / rank / null / column / row spaces -----------------------------

TEST_CASE("matrix: rref of identity is identity", "[9][matrix][rref]") {
    auto I = Matrix::identity(3);
    auto [r, pivots] = I.rref();
    REQUIRE(pivots.size() == 3);
    REQUIRE(r.equals(I));
}

TEST_CASE("matrix: rref of rank-2 3x3", "[9][matrix][rref]") {
    Matrix m = {{integer(1), integer(2), integer(3)},
                {integer(2), integer(4), integer(6)},
                {integer(0), integer(0), integer(1)}};
    // Row 2 is 2 * row 1. Rank should be 2.
    auto [r, pivots] = m.rref();
    REQUIRE(pivots.size() == 2);
    REQUIRE(pivots[0] == 0);
    REQUIRE(pivots[1] == 2);
}

TEST_CASE("matrix: rank reflects row dependence", "[9][matrix][rank]") {
    Matrix m = {{integer(1), integer(2)}, {integer(2), integer(4)}};
    REQUIRE(m.rank() == 1);
    Matrix n = {{integer(1), integer(0)}, {integer(0), integer(1)}};
    REQUIRE(n.rank() == 2);
}

TEST_CASE("matrix: nullspace of rank-deficient 2x2", "[9][matrix][nullspace]") {
    Matrix m = {{integer(1), integer(2)}, {integer(2), integer(4)}};
    auto ns = m.nullspace();
    REQUIRE(ns.size() == 1);
    // Verify A · v = 0 for the basis vector.
    auto Av = m * ns[0];
    REQUIRE(Av.at(0, 0) == S::Zero());
    REQUIRE(Av.at(1, 0) == S::Zero());
}

TEST_CASE("matrix: nullspace of full-rank is empty",
          "[9][matrix][nullspace]") {
    auto I = Matrix::identity(3);
    auto ns = I.nullspace();
    REQUIRE(ns.empty());
}

TEST_CASE("matrix: columnspace of rank-2 3x3 has 2 vectors",
          "[9][matrix][columnspace]") {
    Matrix m = {{integer(1), integer(2), integer(3)},
                {integer(2), integer(4), integer(7)},
                {integer(0), integer(0), integer(1)}};
    auto cs = m.columnspace();
    REQUIRE(cs.size() == 2);
    REQUIRE(cs[0].at(0, 0) == integer(1));  // first pivot column
}

TEST_CASE("matrix: rowspace of rank-1 has 1 vector",
          "[9][matrix][rowspace]") {
    Matrix m = {{integer(1), integer(2)}, {integer(2), integer(4)}};
    auto rs = m.rowspace();
    REQUIRE(rs.size() == 1);
}

// ----- charpoly / eigenvals / eigenvects / diagonalize ----------------------

TEST_CASE("matrix: charpoly of [[2,0],[0,3]] = (λ-2)(λ-3)",
          "[9][matrix][charpoly][oracle]") {
    auto& oracle = Oracle::instance();
    Matrix m = {{integer(2), integer(0)}, {integer(0), integer(3)}};
    auto lam = symbol("lam");
    auto cp = m.charpoly(lam);
    REQUIRE(oracle.equivalent(cp->str(), "(lam - 2)*(lam - 3)"));
}

TEST_CASE("matrix: eigenvals of diagonal [[2, 5]] = {2, 5}",
          "[9][matrix][eigenvals]") {
    Matrix m = {{integer(2), integer(0)}, {integer(0), integer(5)}};
    auto vals = m.eigenvals();
    REQUIRE(vals.size() == 2);
}

TEST_CASE("matrix: eigenvals of [[1, 1], [0, 1]] is {1, 1}",
          "[9][matrix][eigenvals]") {
    // Defective — single eigenvalue 1 with algebraic multiplicity 2.
    Matrix m = {{integer(1), integer(1)}, {integer(0), integer(1)}};
    auto vals = m.eigenvals();
    REQUIRE(vals.size() == 2);
    REQUIRE(vals[0] == integer(1));
    REQUIRE(vals[1] == integer(1));
}

TEST_CASE("matrix: eigenvects of [[2,0],[0,3]] gives e1, e2",
          "[9][matrix][eigenvects]") {
    Matrix m = {{integer(2), integer(0)}, {integer(0), integer(3)}};
    auto evs = m.eigenvects();
    REQUIRE(evs.size() == 2);
    // Each eigenspace is 1-dimensional.
    REQUIRE(evs[0].second.size() == 1);
    REQUIRE(evs[1].second.size() == 1);
}

TEST_CASE("matrix: defective [[1,1],[0,1]] has 1D eigenspace",
          "[9][matrix][eigenvects]") {
    Matrix m = {{integer(1), integer(1)}, {integer(0), integer(1)}};
    auto evs = m.eigenvects();
    // One distinct eigenvalue with geometric multiplicity 1.
    REQUIRE(evs.size() == 1);
    REQUIRE(evs[0].second.size() == 1);
}

TEST_CASE("matrix: is_diagonalizable: diagonal yes, defective no",
          "[9][matrix][diagonalizable]") {
    Matrix d = {{integer(2), integer(0)}, {integer(0), integer(3)}};
    REQUIRE(d.is_diagonalizable());
    Matrix def = {{integer(1), integer(1)}, {integer(0), integer(1)}};
    REQUIRE(!def.is_diagonalizable());
}

TEST_CASE("matrix: diagonalize 2x2 reproduces A = P·D·P⁻¹",
          "[9][matrix][diagonalize]") {
    Matrix m = {{integer(4), integer(1)}, {integer(2), integer(3)}};
    auto [P, D] = m.diagonalize();
    auto reconstructed = P * D * P.inverse();
    // Verify each entry matches via simplify().
    for (std::size_t i = 0; i < 2; ++i) {
        for (std::size_t j = 0; j < 2; ++j) {
            Expr diff = simplify(reconstructed.at(i, j) - m.at(i, j));
            REQUIRE(diff == S::Zero());
        }
    }
}

TEST_CASE("matrix: diagonalize throws on defective matrix",
          "[9][matrix][diagonalize]") {
    Matrix def = {{integer(1), integer(1)}, {integer(0), integer(1)}};
    REQUIRE_THROWS_AS(def.diagonalize(), std::invalid_argument);
}

// ----- LU / QR / Cholesky --------------------------------------------------

TEST_CASE("matrix: LU of 2x2", "[9][matrix][lu]") {
    Matrix m = {{integer(4), integer(3)}, {integer(6), integer(3)}};
    auto [L, U] = m.lu();
    // L should be unit-diagonal lower triangular.
    REQUIRE(L.at(0, 0) == integer(1));
    REQUIRE(L.at(1, 1) == integer(1));
    REQUIRE(L.at(0, 1) == S::Zero());
    // Reconstruction: L · U = A.
    auto rec = L * U;
    REQUIRE(rec.equals(m));
}

TEST_CASE("matrix: LU of 3x3", "[9][matrix][lu]") {
    Matrix m = {{integer(2), integer(1), integer(1)},
                {integer(4), integer(3), integer(3)},
                {integer(8), integer(7), integer(9)}};
    auto [L, U] = m.lu();
    auto rec = L * U;
    REQUIRE(rec.equals(m));
}

TEST_CASE("matrix: QR of 2x2 reconstructs", "[9][matrix][qr]") {
    Matrix m = {{integer(1), integer(0)}, {integer(0), integer(1)}};
    auto [Q, R] = m.qr();
    auto rec = Q * R;
    // Identity in, identity out (up to simplification of sqrt(1)→1).
    for (std::size_t i = 0; i < 2; ++i) {
        for (std::size_t j = 0; j < 2; ++j) {
            Expr d = simplify(rec.at(i, j) - m.at(i, j));
            REQUIRE(d == S::Zero());
        }
    }
}

TEST_CASE("matrix: QR of [[1,1],[0,1]]", "[9][matrix][qr]") {
    Matrix m = {{integer(1), integer(1)}, {integer(0), integer(1)}};
    auto [Q, R] = m.qr();
    auto rec = Q * R;
    for (std::size_t i = 0; i < 2; ++i) {
        for (std::size_t j = 0; j < 2; ++j) {
            Expr d = simplify(rec.at(i, j) - m.at(i, j));
            REQUIRE(d == S::Zero());
        }
    }
}

TEST_CASE("matrix: Cholesky of identity", "[9][matrix][cholesky]") {
    auto I = Matrix::identity(3);
    auto L = I.cholesky();
    // Identity is its own Cholesky factor.
    for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
            Expr d = simplify(L.at(i, j) - I.at(i, j));
            REQUIRE(d == S::Zero());
        }
    }
}

TEST_CASE("matrix: Cholesky of SPD reconstructs A = L·Lᵀ",
          "[9][matrix][cholesky]") {
    Matrix m = {{integer(4), integer(2)}, {integer(2), integer(2)}};
    auto L = m.cholesky();
    auto rec = L * L.transpose();
    for (std::size_t i = 0; i < 2; ++i) {
        for (std::size_t j = 0; j < 2; ++j) {
            Expr d = simplify(rec.at(i, j) - m.at(i, j));
            REQUIRE(d == S::Zero());
        }
    }
}

TEST_CASE("matrix: Cholesky throws on asymmetric matrix",
          "[9][matrix][cholesky]") {
    Matrix m = {{integer(1), integer(2)}, {integer(3), integer(4)}};
    REQUIRE_THROWS_AS(m.cholesky(), std::invalid_argument);
}

#include <sympp/matrices/matrix_symbol.hpp>

// ----- MatrixSymbol + matrix expression tree -------------------------------

TEST_CASE("matrix_symbol: shape and printing", "[9][matrix_symbol]") {
    auto A = matrix_symbol("A", 3, 4);
    REQUIRE(A->rows() == 3);
    REQUIRE(A->cols() == 4);
    REQUIRE(A->str() == "A");
}

TEST_CASE("matrix_symbol: matmul shape product",
          "[9][matrix_symbol]") {
    auto A = matrix_symbol("A", 3, 4);
    auto B = matrix_symbol("B", 4, 2);
    auto C = matmul(A, B);
    REQUIRE(C->rows() == 3);
    REQUIRE(C->cols() == 2);
    REQUIRE(C->str() == "A*B");
}

TEST_CASE("matrix_symbol: matmul throws on inner dim mismatch",
          "[9][matrix_symbol]") {
    auto A = matrix_symbol("A", 3, 4);
    auto B = matrix_symbol("B", 5, 2);
    REQUIRE_THROWS_AS(matmul(A, B), std::invalid_argument);
}

TEST_CASE("matrix_symbol: matadd shape",
          "[9][matrix_symbol]") {
    auto A = matrix_symbol("A", 2, 3);
    auto B = matrix_symbol("B", 2, 3);
    auto C = matadd(A, B);
    REQUIRE(C->rows() == 2);
    REQUIRE(C->cols() == 3);
    REQUIRE(C->str() == "A + B");
}

TEST_CASE("matrix_symbol: transpose flips shape",
          "[9][matrix_symbol]") {
    auto A = matrix_symbol("A", 3, 4);
    auto At = mtranspose(A);
    REQUIRE(At->rows() == 4);
    REQUIRE(At->cols() == 3);
    REQUIRE(At->str() == "A.T");
}

TEST_CASE("matrix_symbol: inverse only on square",
          "[9][matrix_symbol]") {
    auto A = matrix_symbol("A", 3, 4);
    REQUIRE_THROWS_AS(minverse(A), std::invalid_argument);
    auto B = matrix_symbol("B", 3, 3);
    auto Binv = minverse(B);
    REQUIRE(Binv->str() == "B**(-1)");
}

TEST_CASE("matrix_symbol: trace and determinant are scalar",
          "[9][matrix_symbol]") {
    auto A = matrix_symbol("A", 3, 3);
    auto t = mtrace(A);
    REQUIRE(t->rows() == 1);
    REQUIRE(t->cols() == 1);
    REQUIRE(t->str() == "Trace(A)");
    auto d = mdeterminant(A);
    REQUIRE(d->str() == "Determinant(A)");
}

TEST_CASE("matrix_symbol: matpow square only",
          "[9][matrix_symbol]") {
    auto A = matrix_symbol("A", 3, 3);
    auto A3 = matpow(A, 3);
    REQUIRE(A3->str() == "A**3");
    auto B = matrix_symbol("B", 2, 3);
    REQUIRE_THROWS_AS(matpow(B, 2), std::invalid_argument);
}

TEST_CASE("matrix_symbol: nested expressions print with parens",
          "[9][matrix_symbol]") {
    auto A = matrix_symbol("A", 2, 2);
    auto B = matrix_symbol("B", 2, 2);
    auto C = matmul(matadd(A, B), A);
    // (A + B) * A
    REQUIRE(C->str() == "(A + B)*A");
}

TEST_CASE("matrix_symbol: BlockMatrix shape and printing",
          "[9][matrix_symbol][block]") {
    auto A = matrix_symbol("A", 2, 2);
    auto B = matrix_symbol("B", 2, 3);
    auto C = matrix_symbol("C", 4, 2);
    auto D = matrix_symbol("D", 4, 3);
    auto block = block_matrix({{A, B}, {C, D}});
    REQUIRE(block->rows() == 6);
    REQUIRE(block->cols() == 5);
    REQUIRE(block->str() == "BlockMatrix([[A, B], [C, D]])");
}

TEST_CASE("matrix_symbol: BlockMatrix shape mismatch throws",
          "[9][matrix_symbol][block]") {
    auto A = matrix_symbol("A", 2, 2);
    auto B = matrix_symbol("B", 3, 2);  // wrong rows for row 0
    REQUIRE_THROWS_AS(
        block_matrix({{A, B}}),
        std::invalid_argument);
}
