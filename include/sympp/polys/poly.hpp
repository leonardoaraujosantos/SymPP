#pragma once

// Minimal univariate polynomial.
//
// Poly stores coefficients in increasing degree order (coeffs[0] is the
// constant term). `var` is the symbol the polynomial is in.
//
// Aggressive Phase 4 cut per docs/04-roadmap.md: ship the API surface
// callers need (degree, coefficients, leading_coeff, evaluate, +/-/*) and
// closed-form roots up to degree 2; full Berlekamp-Zassenhaus, multivariate
// factor, Groebner stay deferred.
//
// Reference: sympy/polys/polytools.py::Poly

#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class SYMPP_EXPORT Poly {
public:
    // Build a Poly from an existing Expr in `var`. Collects like-degree
    // terms; coefficients with non-numeric / non-polynomial structure stay
    // symbolic in their slot.
    Poly(const Expr& expr, const Expr& var);
    Poly(std::vector<Expr> coeffs, Expr var);

    [[nodiscard]] std::size_t degree() const noexcept;
    [[nodiscard]] const std::vector<Expr>& coeffs() const noexcept { return coeffs_; }
    [[nodiscard]] const Expr& leading_coeff() const;
    [[nodiscard]] const Expr& constant_term() const;
    [[nodiscard]] const Expr& var() const noexcept { return var_; }

    [[nodiscard]] Poly operator+(const Poly& other) const;
    [[nodiscard]] Poly operator-(const Poly& other) const;
    [[nodiscard]] Poly operator*(const Poly& other) const;

    // Evaluate at numeric / symbolic value.
    [[nodiscard]] Expr eval(const Expr& at) const;

    // Reconstitute as an Expr in var.
    [[nodiscard]] Expr as_expr() const;

    [[nodiscard]] std::string str() const;

    // Closed-form roots for degree 1 and 2; returns std::nullopt for higher.
    // For degree 2: roots of a*x^2 + b*x + c.
    [[nodiscard]] std::vector<Expr> roots() const;

private:
    std::vector<Expr> coeffs_;  // increasing degree, [0]=const, [d]=leading
    Expr var_;

    void trim_leading_zeros();
};

}  // namespace sympp
