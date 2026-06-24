#pragma once

// Zeilberger's algorithm (creative telescoping) — closed forms for definite
// sums of a parametric hypergeometric term F(n, k). Internal helper used by
// summation(); not part of the public API.
//
// For S(n) = Σ_{k=lo}^{hi(n)} F(n, k) the algorithm finds the P-recursive
// recurrence Σ_{j=0}^{J} a_j(n)·S(n+j) = 0 that S satisfies (the Zeilberger
// recurrence). When the recurrence is first order it is solved to a closed
// form. Discovery is done by exact-rational fitting over many integer n and
// the result is verified on held-out points before being returned, so it never
// produces a wrong value.
//
// Reference: Petkovšek, Wilf & Zeilberger, "A = B", Ch. 6.

#include <optional>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp::detail {

// Returns the closed form Σ_{var=lo}^{hi} t when t is a parametric
// hypergeometric term (exactly one extra parameter besides var) whose sum
// satisfies a first-order recurrence, else nullopt.
[[nodiscard]] std::optional<Expr> zeilberger_summation(const Expr& t,
                                                       const Expr& var,
                                                       const Expr& lo,
                                                       const Expr& hi);

}  // namespace sympp::detail
