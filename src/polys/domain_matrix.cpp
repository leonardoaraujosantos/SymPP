#include <sympp/polys/domain_matrix.hpp>

#include <stdexcept>
#include <utility>

#include <sympp/core/basic.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {
[[nodiscard]] Expr mpq_to_expr(const mpq_class& q) {
    if (q.get_den() == 1) return make<Integer>(q.get_num());
    return rational(q.get_num(), q.get_den());
}
[[nodiscard]] mpq_class entry_to_mpq(const Expr& e) {
    if (e && e->type_id() == TypeId::Integer) {
        return mpq_class(static_cast<const Integer&>(*e).value());
    }
    if (e && e->type_id() == TypeId::Rational) {
        return static_cast<const Rational&>(*e).value();
    }
    throw std::invalid_argument("DomainMatrix: entries must be Integer or Rational");
}
}  // namespace

DomainMatrix::DomainMatrix(std::size_t rows, std::size_t cols, Domain domain)
    : rows_(rows), cols_(cols), domain_(domain), data_(rows * cols, mpq_class(0)) {}

DomainMatrix DomainMatrix::from_matrix(const Matrix& m) {
    bool all_int = true;
    for (std::size_t i = 0; i < m.rows(); ++i) {
        for (std::size_t j = 0; j < m.cols(); ++j) {
            if (!m.at(i, j) || m.at(i, j)->type_id() != TypeId::Integer) all_int = false;
        }
    }
    DomainMatrix d(m.rows(), m.cols(), all_int ? Domain::ZZ : Domain::QQ);
    for (std::size_t i = 0; i < m.rows(); ++i) {
        for (std::size_t j = 0; j < m.cols(); ++j) {
            mpq_class q = entry_to_mpq(m.at(i, j));
            q.canonicalize();
            d.data_[i * d.cols_ + j] = std::move(q);
        }
    }
    return d;
}

Matrix DomainMatrix::to_matrix() const {
    Matrix m = Matrix::zeros(rows_, cols_);
    for (std::size_t i = 0; i < rows_; ++i) {
        for (std::size_t j = 0; j < cols_; ++j) m.set(i, j, mpq_to_expr(at(i, j)));
    }
    return m;
}

const mpq_class& DomainMatrix::at(std::size_t i, std::size_t j) const {
    return data_[i * cols_ + j];
}
void DomainMatrix::set(std::size_t i, std::size_t j, mpq_class v) {
    v.canonicalize();
    data_[i * cols_ + j] = std::move(v);
}

DomainMatrix DomainMatrix::operator+(const DomainMatrix& o) const {
    if (rows_ != o.rows_ || cols_ != o.cols_) {
        throw std::invalid_argument("DomainMatrix: shape mismatch in +");
    }
    DomainMatrix r(rows_, cols_, domain_ == Domain::ZZ && o.domain_ == Domain::ZZ
                                     ? Domain::ZZ : Domain::QQ);
    for (std::size_t k = 0; k < data_.size(); ++k) r.data_[k] = data_[k] + o.data_[k];
    return r;
}

DomainMatrix DomainMatrix::operator*(const DomainMatrix& o) const {
    if (cols_ != o.rows_) throw std::invalid_argument("DomainMatrix: shape mismatch in *");
    DomainMatrix r(rows_, o.cols_, domain_ == Domain::ZZ && o.domain_ == Domain::ZZ
                                       ? Domain::ZZ : Domain::QQ);
    for (std::size_t i = 0; i < rows_; ++i) {
        for (std::size_t j = 0; j < o.cols_; ++j) {
            mpq_class s = 0;
            for (std::size_t k = 0; k < cols_; ++k) s += at(i, k) * o.at(k, j);
            r.data_[i * o.cols_ + j] = std::move(s);
        }
    }
    return r;
}

Expr DomainMatrix::det() const {
    if (rows_ != cols_) throw std::invalid_argument("DomainMatrix: det needs a square matrix");
    std::size_t n = rows_;
    std::vector<mpq_class> M = data_;
    auto at_ = [&](std::size_t i, std::size_t j) -> mpq_class& { return M[i * n + j]; };
    mpq_class prev = 1;
    int sign = 1;
    for (std::size_t k = 0; k < n; ++k) {
        if (at_(k, k) == 0) {  // pivot
            std::size_t piv = k + 1;
            while (piv < n && at_(piv, k) == 0) ++piv;
            if (piv == n) return make<Integer>(0);
            for (std::size_t j = 0; j < n; ++j) std::swap(at_(k, j), at_(piv, j));
            sign = -sign;
        }
        for (std::size_t i = k + 1; i < n; ++i) {
            for (std::size_t j = k + 1; j < n; ++j) {
                at_(i, j) = (at_(i, j) * at_(k, k) - at_(i, k) * at_(k, j)) / prev;
            }
            at_(i, k) = 0;
        }
        prev = at_(k, k);
    }
    mpq_class d = at_(n - 1, n - 1);
    if (sign < 0) d = -d;
    return mpq_to_expr(d);
}

DomainMatrix DomainMatrix::rref() const {
    DomainMatrix r = *this;
    auto at_ = [&](std::size_t i, std::size_t j) -> mpq_class& { return r.data_[i * cols_ + j]; };
    std::size_t pivot_row = 0;
    for (std::size_t col = 0; col < cols_ && pivot_row < rows_; ++col) {
        std::size_t sel = pivot_row;
        while (sel < rows_ && at_(sel, col) == 0) ++sel;
        if (sel == rows_) continue;
        for (std::size_t j = 0; j < cols_; ++j) std::swap(at_(sel, j), at_(pivot_row, j));
        mpq_class p = at_(pivot_row, col);
        for (std::size_t j = 0; j < cols_; ++j) at_(pivot_row, j) /= p;
        for (std::size_t i = 0; i < rows_; ++i) {
            if (i == pivot_row || at_(i, col) == 0) continue;
            mpq_class f = at_(i, col);
            for (std::size_t j = 0; j < cols_; ++j) at_(i, j) -= f * at_(pivot_row, j);
        }
        ++pivot_row;
    }
    r.domain_ = Domain::QQ;  // row reduction introduces rationals
    return r;
}

std::size_t DomainMatrix::rank() const {
    DomainMatrix r = rref();
    std::size_t rank = 0;
    for (std::size_t i = 0; i < rows_; ++i) {
        for (std::size_t j = 0; j < cols_; ++j) {
            if (r.at(i, j) != 0) { ++rank; break; }
        }
    }
    return rank;
}

}  // namespace sympp
