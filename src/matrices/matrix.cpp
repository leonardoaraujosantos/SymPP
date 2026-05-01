#include <sympp/matrices/matrix.hpp>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/add.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/polys/poly.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/simplify/simplify.hpp>

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

Matrix Matrix::adjugate() const {
    if (!is_square()) {
        throw std::invalid_argument("Matrix: adjugate requires square matrix");
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
    return adj;
}

Matrix Matrix::inverse() const {
    if (!is_square()) {
        throw std::invalid_argument("Matrix: inverse requires square matrix");
    }
    Expr d = det();
    if (d == S::Zero()) {
        throw std::invalid_argument("Matrix: matrix is singular (det = 0)");
    }
    return adjugate().scalar_mul(pow(d, S::NegativeOne()));
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

Matrix Matrix::conjugate_transpose() const {
    Matrix out(cols_, rows_);
    for (std::size_t i = 0; i < rows_; ++i) {
        for (std::size_t j = 0; j < cols_; ++j) {
            out.set(j, i, conjugate(at(i, j)));
        }
    }
    return out;
}

// ----- charpoly / eigenvals / eigenvects / diagonalize ----------------------

Expr Matrix::charpoly(const Expr& lambda_var) const {
    if (!is_square()) {
        throw std::invalid_argument("Matrix: charpoly requires square matrix");
    }
    Matrix lambda_I = Matrix::identity(rows_).scalar_mul(lambda_var);
    return (lambda_I - *this).det();
}

std::vector<Expr> Matrix::eigenvals() const {
    if (!is_square()) {
        throw std::invalid_argument("Matrix: eigenvals requires square matrix");
    }
    auto lam = symbol("__eigval_lambda");
    Expr cp = charpoly(lam);
    Poly p(expand(cp), lam);
    return p.roots();
}

std::vector<std::pair<Expr, std::vector<Matrix>>>
Matrix::eigenvects() const {
    if (!is_square()) {
        throw std::invalid_argument("Matrix: eigenvects requires square matrix");
    }
    std::vector<std::pair<Expr, std::vector<Matrix>>> result;
    auto vals = eigenvals();

    // Group identical eigenvalues — algebraic multiplicity matters for
    // detecting diagonalizability later.
    std::vector<std::pair<Expr, std::size_t>> grouped;
    for (const auto& v : vals) {
        bool merged = false;
        for (auto& [val, count] : grouped) {
            if (val == v) {
                ++count;
                merged = true;
                break;
            }
        }
        if (!merged) grouped.emplace_back(v, 1);
    }

    for (auto& pair : grouped) {
        const Expr& val = pair.first;
        Matrix M(rows_, cols_);
        for (std::size_t i = 0; i < rows_; ++i) {
            for (std::size_t j = 0; j < cols_; ++j) {
                Expr entry = (i == j) ? (at(i, j) - val) : at(i, j);
                M.set(i, j, simplify(entry));
            }
        }
        auto basis = M.nullspace();
        result.emplace_back(val, std::move(basis));
    }
    return result;
}

bool Matrix::is_diagonalizable() const {
    if (!is_square()) return false;
    auto evs = eigenvects();
    std::size_t total = 0;
    for (const auto& pair : evs) total += pair.second.size();
    return total == rows_;
}

std::pair<Matrix, Matrix> Matrix::diagonalize() const {
    auto evs = eigenvects();
    std::vector<Matrix> cols;
    std::vector<Expr> diag_entries;
    for (const auto& pair : evs) {
        const Expr& val = pair.first;
        for (const auto& v : pair.second) {
            cols.push_back(v);
            diag_entries.push_back(val);
        }
    }
    if (cols.size() != rows_) {
        throw std::invalid_argument(
            "Matrix: not diagonalizable in current expression domain");
    }
    Matrix P(rows_, rows_);
    for (std::size_t j = 0; j < cols.size(); ++j) {
        for (std::size_t i = 0; i < rows_; ++i) {
            P.set(i, j, cols[j].at(i, 0));
        }
    }
    Matrix D = Matrix::zeros(rows_, rows_);
    for (std::size_t i = 0; i < rows_; ++i) {
        D.set(i, i, diag_entries[i]);
    }
    return {P, D};
}

// ----- rref / rank / nullspace / columnspace / rowspace ---------------------

std::pair<Matrix, std::vector<std::size_t>> Matrix::rref() const {
    Matrix m = *this;  // working copy
    std::vector<std::size_t> pivots;
    std::size_t r = 0;  // current row
    for (std::size_t c = 0; c < cols_ && r < rows_; ++c) {
        // Find a pivot row >= r whose entry at column c is nonzero.
        std::size_t pivot_row = rows_;
        for (std::size_t i = r; i < rows_; ++i) {
            // simplify here folds e.g. (a-a) to 0 so the test is robust.
            Expr v = simplify(m.at(i, c));
            if (!(v == S::Zero())) {
                pivot_row = i;
                break;
            }
        }
        if (pivot_row == rows_) continue;  // no pivot in this column
        // Swap pivot_row and r.
        if (pivot_row != r) {
            for (std::size_t j = 0; j < cols_; ++j) {
                Expr tmp = m.at(r, j);
                m.set(r, j, m.at(pivot_row, j));
                m.set(pivot_row, j, tmp);
            }
        }
        // Scale row r so the pivot becomes 1.
        Expr pivot = m.at(r, c);
        for (std::size_t j = 0; j < cols_; ++j) {
            m.set(r, j, simplify(m.at(r, j) / pivot));
        }
        // Eliminate this column in every other row.
        for (std::size_t i = 0; i < rows_; ++i) {
            if (i == r) continue;
            Expr factor = m.at(i, c);
            if (factor == S::Zero()) continue;
            for (std::size_t j = 0; j < cols_; ++j) {
                Expr v = m.at(i, j) - factor * m.at(r, j);
                m.set(i, j, simplify(v));
            }
        }
        pivots.push_back(c);
        ++r;
    }
    return {m, pivots};
}

std::size_t Matrix::rank() const {
    return rref().second.size();
}

std::vector<Matrix> Matrix::nullspace() const {
    auto [r, pivots] = rref();
    // Free columns are those not in pivots.
    std::vector<std::size_t> free_cols;
    {
        std::size_t pi = 0;
        for (std::size_t c = 0; c < cols_; ++c) {
            if (pi < pivots.size() && pivots[pi] == c) ++pi;
            else free_cols.push_back(c);
        }
    }
    std::vector<Matrix> basis;
    basis.reserve(free_cols.size());
    for (std::size_t free : free_cols) {
        // Construct a vector with x_free = 1, other free vars = 0,
        // pivot vars = -r[i, free].
        Matrix v(cols_, 1);
        for (std::size_t i = 0; i < cols_; ++i) v.set(i, 0, S::Zero());
        v.set(free, 0, S::One());
        for (std::size_t pi = 0; pi < pivots.size(); ++pi) {
            v.set(pivots[pi], 0, mul(S::NegativeOne(), r.at(pi, free)));
        }
        basis.push_back(std::move(v));
    }
    return basis;
}

std::vector<Matrix> Matrix::columnspace() const {
    auto [_, pivots] = rref();
    std::vector<Matrix> basis;
    basis.reserve(pivots.size());
    for (std::size_t c : pivots) {
        Matrix v(rows_, 1);
        for (std::size_t i = 0; i < rows_; ++i) v.set(i, 0, at(i, c));
        basis.push_back(std::move(v));
    }
    return basis;
}

std::vector<Matrix> Matrix::rowspace() const {
    auto [r, pivots] = rref();
    std::vector<Matrix> basis;
    basis.reserve(pivots.size());
    for (std::size_t i = 0; i < pivots.size(); ++i) {
        Matrix row_vec(1, cols_);
        for (std::size_t j = 0; j < cols_; ++j) row_vec.set(0, j, r.at(i, j));
        basis.push_back(std::move(row_vec));
    }
    return basis;
}

Expr Matrix::norm_frobenius() const {
    // sqrt(Σ |aᵢⱼ|²) — uses |z|² = z·conj(z) so the result handles
    // complex entries correctly. For real entries this reduces to
    // sqrt(Σ aᵢⱼ²).
    Expr sum = S::Zero();
    for (const auto& v : data_) {
        Expr term = mul(v, conjugate(v));
        sum = add(sum, term);
    }
    return sqrt(sum);
}

// ----- jacobian / gradient / hessian / wronskian -----------------------------

Matrix jacobian(const std::vector<Expr>& fs, const std::vector<Expr>& vars) {
    Matrix out(fs.size(), vars.size());
    for (std::size_t i = 0; i < fs.size(); ++i) {
        for (std::size_t j = 0; j < vars.size(); ++j) {
            out.set(i, j, diff(fs[i], vars[j]));
        }
    }
    return out;
}

Matrix gradient(const Expr& f, const std::vector<Expr>& vars) {
    Matrix g(vars.size(), 1);
    for (std::size_t i = 0; i < vars.size(); ++i) {
        g.set(i, 0, diff(f, vars[i]));
    }
    return g;
}

Matrix hessian(const Expr& f, const std::vector<Expr>& vars) {
    const std::size_t n = vars.size();
    Matrix h(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        Expr di = diff(f, vars[i]);
        for (std::size_t j = 0; j < n; ++j) {
            h.set(i, j, diff(di, vars[j]));
        }
    }
    return h;
}

Expr wronskian(const std::vector<Expr>& fs, const Expr& var) {
    const std::size_t n = fs.size();
    Matrix w(n, n);
    std::vector<Expr> current = fs;  // row 0 = fs themselves
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            w.set(i, j, current[j]);
        }
        if (i + 1 < n) {
            for (auto& v : current) v = diff(v, var);
        }
    }
    return w.det();
}

// ----- Constructors ---------------------------------------------------------

Matrix hilbert(std::size_t n) {
    Matrix h(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            h.set(i, j, rational(1, static_cast<long>(i + j + 1)));
        }
    }
    return h;
}

Matrix vandermonde(const std::vector<Expr>& xs, std::size_t cols) {
    Matrix v(xs.size(), cols);
    for (std::size_t i = 0; i < xs.size(); ++i) {
        for (std::size_t j = 0; j < cols; ++j) {
            v.set(i, j, pow(xs[i], integer(static_cast<long>(j))));
        }
    }
    return v;
}

Matrix companion(const std::vector<Expr>& coeffs) {
    // For monic p(x) = x^n + c_{n-1} x^(n-1) + ... + c_0, the companion has
    // 1's on the sub-diagonal and -c_i in the last column.
    const std::size_t n = coeffs.size();
    if (n == 0) {
        throw std::invalid_argument("companion: empty coefficient list");
    }
    Matrix c = Matrix::zeros(n, n);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        c.set(i + 1, i, S::One());
    }
    for (std::size_t i = 0; i < n; ++i) {
        c.set(i, n - 1, mul(S::NegativeOne(), coeffs[i]));
    }
    return c;
}

Matrix rotation_matrix_2d(const Expr& theta) {
    Matrix r(2, 2);
    Expr c = cos(theta);
    Expr s = sin(theta);
    r.set(0, 0, c);
    r.set(0, 1, mul(S::NegativeOne(), s));
    r.set(1, 0, s);
    r.set(1, 1, c);
    return r;
}

Matrix rotation_matrix_x(const Expr& theta) {
    Matrix r = Matrix::identity(3);
    Expr c = cos(theta);
    Expr s = sin(theta);
    r.set(1, 1, c);
    r.set(1, 2, mul(S::NegativeOne(), s));
    r.set(2, 1, s);
    r.set(2, 2, c);
    return r;
}

Matrix rotation_matrix_y(const Expr& theta) {
    Matrix r = Matrix::identity(3);
    Expr c = cos(theta);
    Expr s = sin(theta);
    r.set(0, 0, c);
    r.set(0, 2, s);
    r.set(2, 0, mul(S::NegativeOne(), s));
    r.set(2, 2, c);
    return r;
}

Matrix rotation_matrix_z(const Expr& theta) {
    Matrix r = Matrix::identity(3);
    Expr c = cos(theta);
    Expr s = sin(theta);
    r.set(0, 0, c);
    r.set(0, 1, mul(S::NegativeOne(), s));
    r.set(1, 0, s);
    r.set(1, 1, c);
    return r;
}

// ----- Hadamard / Kronecker -------------------------------------------------

Matrix hadamard_product(const Matrix& a, const Matrix& b) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) {
        throw std::invalid_argument("hadamard_product: shape mismatch");
    }
    Matrix out(a.rows(), a.cols());
    for (std::size_t i = 0; i < a.rows(); ++i) {
        for (std::size_t j = 0; j < a.cols(); ++j) {
            out.set(i, j, mul(a.at(i, j), b.at(i, j)));
        }
    }
    return out;
}

Matrix kronecker_product(const Matrix& a, const Matrix& b) {
    const std::size_t ar = a.rows(), ac = a.cols();
    const std::size_t br = b.rows(), bc = b.cols();
    Matrix out(ar * br, ac * bc);
    for (std::size_t i = 0; i < ar; ++i) {
        for (std::size_t j = 0; j < ac; ++j) {
            const Expr& aij = a.at(i, j);
            for (std::size_t k = 0; k < br; ++k) {
                for (std::size_t l = 0; l < bc; ++l) {
                    out.set(i * br + k, j * bc + l, mul(aij, b.at(k, l)));
                }
            }
        }
    }
    return out;
}

}  // namespace sympp
