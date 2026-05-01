#pragma once

// RootOf — lazy/unevaluated form for the k-th root of a univariate polynomial.
//
// Used when a polynomial's root cannot be expressed cleanly in radicals
// (degree ≥ 5 with no rational roots, or degree 3/4 where the radical form
// is unwieldy). Carries (poly_expr, var, index) and stays unevaluated until
// numerically materialized via evalf.
//
// Indexing: 0-based, real roots first in ascending order, then complex
// roots paired by conjugates. The minimal evalf currently handles real
// roots only via bisection on isolating intervals.
//
// Reference: sympy/polys/rootoftools.py::CRootOf

#include <array>
#include <cstddef>
#include <span>
#include <string>

#include <sympp/core/api.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class SYMPP_EXPORT RootOf final : public Basic {
public:
    RootOf(Expr poly_expr, Expr var, std::size_t index);

    [[nodiscard]] TypeId type_id() const noexcept override {
        return TypeId::RootOf;
    }
    [[nodiscard]] std::size_t hash() const noexcept override { return hash_; }
    [[nodiscard]] std::span<const Expr> args() const noexcept override {
        return {args_.data(), args_.size()};
    }
    [[nodiscard]] bool equals(const Basic& other) const noexcept override;
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;

    [[nodiscard]] const Expr& poly_expr() const { return args_[0]; }
    [[nodiscard]] const Expr& var() const { return args_[1]; }
    [[nodiscard]] std::size_t index() const { return index_; }

    // Numeric evaluation of the index-th real root via Cauchy-bound sampling
    // + MPFR bisection. Returns std::nullopt when the polynomial has
    // non-numeric coefficients or fewer real roots than requested.
    [[nodiscard]] std::optional<Expr> try_evalf(int dps) const;

private:
    std::array<Expr, 2> args_;  // [poly_expr, var]
    std::size_t index_;
    std::size_t hash_;
};

// Factory: returns a RootOf Expr. If the polynomial has fewer real roots than
// the requested index, throws std::out_of_range.
[[nodiscard]] SYMPP_EXPORT Expr root_of(const Expr& poly_expr, const Expr& var,
                                          std::size_t index);

}  // namespace sympp
