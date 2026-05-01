#include <sympp/matrices/matrix.hpp>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/add.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/singletons.hpp>

namespace sympp {

Matrix::Matrix(std::size_t rows, std::size_t cols)
    : rows_(rows), cols_(cols), data_(rows * cols, S::Zero()) {}

Matrix::Matrix(std::size_t rows, std::size_t cols, std::vector<Expr> data)
    : rows_(rows), cols_(cols), data_(std::move(data)) {
    if (data_.size() != rows_ * cols_) {
        throw std::invalid_argument("Matrix: data size mismatch");
    }
}

Matrix::Matrix(std::initializer_list<std::initializer_list<Expr>> rows) {
    rows_ = rows.size();
    cols_ = rows_ ? rows.begin()->size() : 0;
    data_.reserve(rows_ * cols_);
    for (const auto& r : rows) {
        if (r.size() != cols_) {
            throw std::invalid_argument("Matrix: ragged initializer");
        }
        for (const auto& v : r) data_.push_back(v);
    }
}

Matrix Matrix::zeros(std::size_t rows, std::size_t cols) {
    return Matrix(rows, cols);
}

Matrix Matrix::identity(std::size_t n) {
    Matrix out(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        out.set(i, i, S::One());
    }
    return out;
}

const Expr& Matrix::at(std::size_t r, std::size_t c) const {
    if (r >= rows_ || c >= cols_) {
        throw std::out_of_range("Matrix: index out of range");
    }
    return data_[r * cols_ + c];
}

void Matrix::set(std::size_t r, std::size_t c, Expr v) {
    if (r >= rows_ || c >= cols_) {
        throw std::out_of_range("Matrix: index out of range");
    }
    data_[r * cols_ + c] = std::move(v);
}

Matrix Matrix::operator+(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::invalid_argument("Matrix: shape mismatch in +");
    }
    Matrix out(rows_, cols_);
    for (std::size_t i = 0; i < data_.size(); ++i) {
        out.data_[i] = add(data_[i], other.data_[i]);
    }
    return out;
}

Matrix Matrix::operator-(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::invalid_argument("Matrix: shape mismatch in -");
    }
    Matrix out(rows_, cols_);
    for (std::size_t i = 0; i < data_.size(); ++i) {
        out.data_[i] = data_[i] - other.data_[i];
    }
    return out;
}

Matrix Matrix::operator*(const Matrix& other) const {
    if (cols_ != other.rows_) {
        throw std::invalid_argument("Matrix: shape mismatch in *");
    }
    Matrix out(rows_, other.cols_);
    for (std::size_t i = 0; i < rows_; ++i) {
        for (std::size_t j = 0; j < other.cols_; ++j) {
            std::vector<Expr> terms;
            terms.reserve(cols_);
            for (std::size_t k = 0; k < cols_; ++k) {
                terms.push_back(mul(at(i, k), other.at(k, j)));
            }
            out.set(i, j, add(std::move(terms)));
        }
    }
    return out;
}

Matrix Matrix::scalar_mul(const Expr& s) const {
    return map([&](const Expr& v) { return mul(s, v); });
}

Matrix Matrix::transpose() const {
    Matrix out(cols_, rows_);
    for (std::size_t i = 0; i < rows_; ++i) {
        for (std::size_t j = 0; j < cols_; ++j) {
            out.set(j, i, at(i, j));
        }
    }
    return out;
}

Expr Matrix::trace() const {
    if (!is_square()) {
        throw std::invalid_argument("Matrix: trace requires square matrix");
    }
    std::vector<Expr> terms;
    terms.reserve(rows_);
    for (std::size_t i = 0; i < rows_; ++i) terms.push_back(at(i, i));
    return add(std::move(terms));
}

namespace {

[[nodiscard]] Matrix minor_matrix(const Matrix& m, std::size_t row, std::size_t col) {
    std::size_t n = m.rows();
    Matrix out(n - 1, n - 1);
    for (std::size_t i = 0, oi = 0; i < n; ++i) {
        if (i == row) continue;
        for (std::size_t j = 0, oj = 0; j < n; ++j) {
            if (j == col) continue;
            out.set(oi, oj, m.at(i, j));
            ++oj;
        }
        ++oi;
    }
    return out;
}

[[nodiscard]] Expr det_recursive(const Matrix& m) {
    std::size_t n = m.rows();
    if (n == 0) return S::One();   // det of empty matrix is 1 (vacuous product)
    if (n == 1) return m.at(0, 0);
    if (n == 2) {
        return m.at(0, 0) * m.at(1, 1) - m.at(0, 1) * m.at(1, 0);
    }
    // Cofactor expansion along the first row.
    std::vector<Expr> terms;
    terms.reserve(n);
    for (std::size_t j = 0; j < n; ++j) {
        Expr sign = (j % 2 == 0) ? S::One() : S::NegativeOne();
        Expr cofactor = det_recursive(minor_matrix(m, 0, j));
        terms.push_back(mul({sign, m.at(0, j), cofactor}));
    }
    return add(std::move(terms));
}

}  // namespace

Expr Matrix::det() const {
    if (!is_square()) {
        throw std::invalid_argument("Matrix: det requires square matrix");
    }
    if (rows_ == 0) return S::One();
    return det_recursive(*this);
}

Matrix Matrix::inverse() const {
    if (!is_square()) {
        throw std::invalid_argument("Matrix: inverse requires square matrix");
    }
    Expr d = det();
    if (d == S::Zero()) {
        throw std::invalid_argument("Matrix: matrix is singular (det = 0)");
    }
    std::size_t n = rows_;
    Matrix adj(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            Expr sign = ((i + j) % 2 == 0) ? S::One() : S::NegativeOne();
            // Cofactor of (j, i) gives the (i, j) entry of adjugate (transposed).
            Expr c = det_recursive(minor_matrix(*this, j, i));
            adj.set(i, j, mul(sign, c));
        }
    }
    return adj.scalar_mul(pow(d, S::NegativeOne()));
}

bool Matrix::equals(const Matrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) return false;
    for (std::size_t i = 0; i < data_.size(); ++i) {
        if (!(data_[i] == other.data_[i])) return false;
    }
    return true;
}

std::string Matrix::str() const {
    std::string out = "Matrix([\n";
    for (std::size_t i = 0; i < rows_; ++i) {
        out += "  [";
        for (std::size_t j = 0; j < cols_; ++j) {
            if (j > 0) out += ", ";
            out += at(i, j)->str();
        }
        out += "]";
        if (i + 1 < rows_) out += ",";
        out += "\n";
    }
    out += "])";
    return out;
}

// ----- jacobian --------------------------------------------------------------

Matrix jacobian(const std::vector<Expr>& fs, const std::vector<Expr>& vars) {
    Matrix out(fs.size(), vars.size());
    for (std::size_t i = 0; i < fs.size(); ++i) {
        for (std::size_t j = 0; j < vars.size(); ++j) {
            out.set(i, j, diff(fs[i], vars[j]));
        }
    }
    return out;
}

}  // namespace sympp
