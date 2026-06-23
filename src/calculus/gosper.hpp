#pragma once

// Gosper's algorithm — indefinite/definite summation of hypergeometric terms.
// Internal helper used by summation(); not part of the public API.
//
// A term t(k) is hypergeometric when t(k+1)/t(k) is a rational function of k
// (true for products of factorials, gammas, binomials, rising factorials,
// powers cᵏ, and polynomials). Gosper decides whether such a term is summable
// in closed form and, if so, produces the antidifference, giving
//   Σ_{k=lo}^{hi} t(k) = z(hi+1) − z(lo),  z(k) = y(k)·t(k), y rational.
//
// Reference: Petkovšek, Wilf & Zeilberger, "A = B", Ch. 5.

#include <optional>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp::detail {

// Returns the closed-form Σ_{var=lo}^{hi} t when t is Gosper-summable, else
// nullopt. Finite ranges only (hi must not be ∞). The result is verified by the
// rational Gosper identity before being returned, so it is never wrong.
[[nodiscard]] std::optional<Expr> gosper_summation(const Expr& t,
                                                   const Expr& var,
                                                   const Expr& lo,
                                                   const Expr& hi);

}  // namespace sympp::detail
