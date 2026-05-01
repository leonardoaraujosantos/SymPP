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
#include <tuple>

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

    // Polynomial long division. Returns (quotient, remainder) such that
    // self == q * other + r and degree(r) < degree(other) (or r == 0).
    // Throws on zero divisor or variable mismatch.
    [[nodiscard]] std::pair<Poly, Poly> divmod(const Poly& other) const;
    [[nodiscard]] Poly operator/(const Poly& other) const;
    [[nodiscard]] Poly operator%(const Poly& other) const;

    [[nodiscard]] bool is_zero() const;
    // Divides through by the leading coefficient — returns a monic Poly with
    // the same roots. Caller must ensure the leading coefficient is invertible
    // (numeric, or symbolic-but-nonzero).
    [[nodiscard]] Poly monic() const;
    // Differentiate w.r.t. var, returning a Poly of degree one less.
    [[nodiscard]] Poly diff() const;

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

// Univariate polynomial GCD via Euclidean algorithm. Result is monic
// (leading coefficient 1) for non-zero inputs. Behaves correctly over ℚ.
//
// Reference: sympy/polys/euclidtools.py::dup_gcd
[[nodiscard]] SYMPP_EXPORT Poly gcd(const Poly& a, const Poly& b);

// Rational root finder. Returns roots of the polynomial that are exactly
// rational (i.e., p/q for integers p, q), repeated by multiplicity. Works
// only when every coefficient of `f` is an Integer or Rational; otherwise
// returns an empty vector.
//
// Reference: sympy/polys/polyroots.py::roots_quintic / rational roots theorem
[[nodiscard]] SYMPP_EXPORT std::vector<Expr> rational_roots(const Poly& f);

// Square-free factorization via Yun's algorithm.
//
// Returns (content, [(factor_i, multiplicity_i), ...]) such that
//   f = content * Π factor_i^multiplicity_i
// where each factor_i is monic and square-free, and any two distinct
// factor_i are coprime. The trivial factor 1 is omitted.
//
// Reference: sympy/polys/sqfreetools.py::dup_sqf_list
struct SqfList {
    Expr content;
    std::vector<std::pair<Poly, std::size_t>> factors;
};
[[nodiscard]] SYMPP_EXPORT SqfList sqf_list(const Poly& f);

}  // namespace sympp
