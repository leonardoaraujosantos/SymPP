#pragma once

// Matrix<Expr> — dense symbolic matrix.
//
// The Matrix is *not* a Basic node — it sits one level above Expr because
// matrix algorithms benefit from random access to mutable storage. For
// operations that produce another Matrix (transpose, inverse, sum) we
// return a fresh Matrix.
//
// Subset shipped here:
//   * construction from rows or by lambda
//   * identity, zeros
//   * shape, element access
//   * +, -, *, scalar *
//   * transpose, trace
//   * det (recursive cofactor — fine for small matrices, replaceable later)
//   * inv via cofactor / det
//   * jacobian
//
// Reference: sympy/matrices/{matrices,dense,determinant,inverse}.py

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class SYMPP_EXPORT Matrix {
public:
    Matrix(std::size_t rows, std::size_t cols);
    Matrix(std::size_t rows, std::size_t cols, std::vector<Expr> data);
    Matrix(std::initializer_list<std::initializer_list<Expr>> rows);

    [[nodiscard]] static Matrix zeros(std::size_t rows, std::size_t cols);
    [[nodiscard]] static Matrix identity(std::size_t n);

    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t cols() const noexcept { return cols_; }
    [[nodiscard]] bool is_square() const noexcept { return rows_ == cols_; }

    [[nodiscard]] const Expr& at(std::size_t r, std::size_t c) const;
    void set(std::size_t r, std::size_t c, Expr v);

    [[nodiscard]] Matrix operator+(const Matrix& other) const;
    [[nodiscard]] Matrix operator-(const Matrix& other) const;
    [[nodiscard]] Matrix operator*(const Matrix& other) const;
    [[nodiscard]] Matrix scalar_mul(const Expr& s) const;

    [[nodiscard]] Matrix transpose() const;
    [[nodiscard]] Expr trace() const;
    [[nodiscard]] Expr det() const;
    [[nodiscard]] Matrix inverse() const;
    [[nodiscard]] bool equals(const Matrix& other) const;

    // Apply a unary transform to every element.
    template <typename F>
    [[nodiscard]] Matrix map(F&& fn) const {
        std::vector<Expr> out;
        out.reserve(data_.size());
        for (const auto& v : data_) out.push_back(fn(v));
        return Matrix{rows_, cols_, std::move(out)};
    }

    [[nodiscard]] std::string str() const;

private:
    std::size_t rows_ = 0;
    std::size_t cols_ = 0;
    std::vector<Expr> data_;  // row-major
};

// jacobian: F = [f1, ..., fm] w.r.t. vars = [x1, ..., xn] -> m×n matrix
// of partial derivatives.
[[nodiscard]] SYMPP_EXPORT Matrix jacobian(const std::vector<Expr>& fs,
                                            const std::vector<Expr>& vars);

// ---- Standard matrix constructors -----------------------------------------

// Hilbert matrix H[i, j] = 1 / (i + j + 1) for i, j in [0, n).
[[nodiscard]] SYMPP_EXPORT Matrix hilbert(std::size_t n);

// Vandermonde matrix V[i, j] = xs[i]^j.
[[nodiscard]] SYMPP_EXPORT Matrix vandermonde(const std::vector<Expr>& xs,
                                                std::size_t cols);

// Companion matrix of a monic polynomial with coefficients
// [c_0, c_1, ..., c_{n-1}] (lowest-degree first; the polynomial is
// x^n + c_{n-1}*x^(n-1) + ... + c_0). Returns the n×n matrix whose
// characteristic polynomial equals the input.
[[nodiscard]] SYMPP_EXPORT Matrix companion(const std::vector<Expr>& coeffs);

// 2D rotation by angle θ.
[[nodiscard]] SYMPP_EXPORT Matrix rotation_matrix_2d(const Expr& theta);

// 3D rotation about an axis ('x', 'y', or 'z') by angle θ.
[[nodiscard]] SYMPP_EXPORT Matrix rotation_matrix_x(const Expr& theta);
[[nodiscard]] SYMPP_EXPORT Matrix rotation_matrix_y(const Expr& theta);
[[nodiscard]] SYMPP_EXPORT Matrix rotation_matrix_z(const Expr& theta);

// ---- Element-wise / structural products -----------------------------------

// Hadamard (element-wise) product: A ⊙ B with A, B same shape.
[[nodiscard]] SYMPP_EXPORT Matrix hadamard_product(const Matrix& a,
                                                     const Matrix& b);

// Kronecker product: A ⊗ B has shape (rows(A)*rows(B), cols(A)*cols(B)).
[[nodiscard]] SYMPP_EXPORT Matrix kronecker_product(const Matrix& a,
                                                     const Matrix& b);

}  // namespace sympp
