#pragma once

// DomainMatrix — a matrix over a concrete computational domain (ℤ or ℚ) with
// exact, fraction-free algorithms. Operating directly on GMP integers/rationals
// (rather than the general Expr tree) makes determinant / rank / rref on numeric
// matrices much faster and keeps ℤ computations fraction-free (Bareiss).
//
//   auto dm = DomainMatrix::from_matrix(M);   // M of Integer/Rational entries
//   dm.det();                                 // exact, fraction-free
//
// Reference: sympy/polys/matrices/domainmatrix.py.

#include <cstddef>
#include <vector>

#include <gmpxx.h>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>
#include <sympp/matrices/matrix.hpp>

namespace sympp {

class SYMPP_EXPORT DomainMatrix {
public:
    enum class Domain { ZZ, QQ };

    DomainMatrix(std::size_t rows, std::size_t cols, Domain domain);

    // Build from a Matrix of Integer/Rational entries (domain inferred: ℤ if all
    // entries are integers, else ℚ). Throws on a non-numeric entry.
    [[nodiscard]] static DomainMatrix from_matrix(const Matrix& m);
    [[nodiscard]] Matrix to_matrix() const;

    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t cols() const noexcept { return cols_; }
    [[nodiscard]] Domain domain() const noexcept { return domain_; }
    [[nodiscard]] const mpq_class& at(std::size_t i, std::size_t j) const;
    void set(std::size_t i, std::size_t j, mpq_class v);

    [[nodiscard]] DomainMatrix operator+(const DomainMatrix& o) const;
    [[nodiscard]] DomainMatrix operator*(const DomainMatrix& o) const;

    // Determinant via Bareiss fraction-free elimination (exact; integer-only
    // intermediates over ℤ). Returns an Integer or Rational Expr.
    [[nodiscard]] Expr det() const;
    [[nodiscard]] std::size_t rank() const;
    [[nodiscard]] DomainMatrix rref() const;

private:
    std::size_t rows_, cols_;
    Domain domain_;
    std::vector<mpq_class> data_;  // row-major
};

}  // namespace sympp
