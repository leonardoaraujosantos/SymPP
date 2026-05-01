#pragma once

// Convenience free-function queries — thin wrappers over Basic::ask.
// Use these in user code; they read like SymPy's `expr.is_real` / `is_positive`.
//
//     if (is_positive(e) == true) { ... }    // definitely positive
//     if (is_real(e) == false) { ... }       // definitely not real
//     if (!is_integer(e).has_value()) { ... } // unknown
//
// Reference: sympy/core/assumptions.py — `Basic.is_<property>` accessors.

#include <optional>

#include <sympp/core/api.hpp>
#include <sympp/core/assumption_key.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

// Direct query that also applies the standard implications between predicates,
// so callers don't have to. For example:
//
//   * an Expr whose direct ask(Positive) == true also answers ask(Nonzero),
//     ask(Nonnegative), and ask(Real) as true;
//   * an Expr whose direct ask(Integer) == true also answers ask(Rational)
//     and ask(Real) as true.
//
// Each node's ask() implementation is responsible only for the facts it
// can determine *directly* from its structure; this function fills in the
// trivial logical consequences.
[[nodiscard]] SYMPP_EXPORT std::optional<bool> ask(const Expr& e, AssumptionKey k) noexcept;

[[nodiscard]] inline std::optional<bool> is_real(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Real);
}
[[nodiscard]] inline std::optional<bool> is_rational(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Rational);
}
[[nodiscard]] inline std::optional<bool> is_integer(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Integer);
}
[[nodiscard]] inline std::optional<bool> is_positive(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Positive);
}
[[nodiscard]] inline std::optional<bool> is_negative(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Negative);
}
[[nodiscard]] inline std::optional<bool> is_zero(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Zero);
}
[[nodiscard]] inline std::optional<bool> is_nonzero(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Nonzero);
}
[[nodiscard]] inline std::optional<bool> is_nonnegative(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Nonnegative);
}
[[nodiscard]] inline std::optional<bool> is_nonpositive(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Nonpositive);
}
[[nodiscard]] inline std::optional<bool> is_finite(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Finite);
}

}  // namespace sympp
