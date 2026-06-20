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

// Like ask() but WITHOUT implication-chasing: the node's own fact for `k`, with
// any active `assuming` scope taking precedence. Compound nodes (Mul, Add, …)
// query their children through this so scoped facts (e.g. x > 0) propagate
// through products and sums. With no `assuming` scope active it is exactly
// `e->ask(k)`, so existing behavior is unchanged.
[[nodiscard]] SYMPP_EXPORT std::optional<bool> direct_ask(const Expr& e, AssumptionKey k) noexcept;

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
[[nodiscard]] inline std::optional<bool> is_even(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Even);
}
[[nodiscard]] inline std::optional<bool> is_odd(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Odd);
}
[[nodiscard]] inline std::optional<bool> is_prime(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Prime);
}
[[nodiscard]] inline std::optional<bool> is_composite(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Composite);
}
[[nodiscard]] inline std::optional<bool> is_irrational(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Irrational);
}
[[nodiscard]] inline std::optional<bool> is_algebraic(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Algebraic);
}
[[nodiscard]] inline std::optional<bool> is_transcendental(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Transcendental);
}
[[nodiscard]] inline std::optional<bool> is_extended_real(const Expr& e) noexcept {
    return ask(e, AssumptionKey::ExtendedReal);
}
[[nodiscard]] inline std::optional<bool> is_infinite(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Infinite);
}
[[nodiscard]] inline std::optional<bool> is_imaginary(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Imaginary);
}
[[nodiscard]] inline std::optional<bool> is_complex(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Complex);
}
[[nodiscard]] inline std::optional<bool> is_commutative(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Commutative);
}
[[nodiscard]] inline std::optional<bool> is_extended_positive(const Expr& e) noexcept {
    return ask(e, AssumptionKey::ExtendedPositive);
}
[[nodiscard]] inline std::optional<bool> is_extended_negative(const Expr& e) noexcept {
    return ask(e, AssumptionKey::ExtendedNegative);
}
[[nodiscard]] inline std::optional<bool> is_extended_nonnegative(const Expr& e) noexcept {
    return ask(e, AssumptionKey::ExtendedNonnegative);
}
[[nodiscard]] inline std::optional<bool> is_extended_nonpositive(const Expr& e) noexcept {
    return ask(e, AssumptionKey::ExtendedNonpositive);
}
[[nodiscard]] inline std::optional<bool> is_hermitian(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Hermitian);
}
[[nodiscard]] inline std::optional<bool> is_antihermitian(const Expr& e) noexcept {
    return ask(e, AssumptionKey::Antihermitian);
}

// Structural parity of an integer-valued expression. Conservative: reports
// even/odd only when provable from the form — `2·n` (n a known integer) is even,
// `2·n+1` is odd, while a bare integer symbol has unknown parity. (These are
// pure-form deductions; there is no Even/Odd assumption key yet.)
[[nodiscard]] SYMPP_EXPORT bool is_provably_even(const Expr& e) noexcept;
[[nodiscard]] SYMPP_EXPORT bool is_provably_odd(const Expr& e) noexcept;

}  // namespace sympp
