#include <sympp/core/queries.hpp>

#include <optional>

#include <gmpxx.h>

#include <sympp/core/assumption_context.hpp>
#include <sympp/core/assumption_key.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

// Convenience: ask the underlying node directly (no implication-chasing) so
// we don't loop through this very function while computing derivations. An
// active `assuming` scope takes precedence, so a scoped fact (e.g. x > 0)
// overrides the symbol's own mask and feeds the implication chain below.
[[nodiscard]] inline std::optional<bool> direct(const Expr& e, AssumptionKey k) noexcept {
    if (assumption_context_active()) {
        if (auto v = assumption_context_get(e, k); v.has_value()) return v;
    }
    return e ? e->ask(k) : std::nullopt;
}

}  // namespace

std::optional<bool> direct_ask(const Expr& e, AssumptionKey k) noexcept {
    return direct(e, k);
}

std::optional<bool> ask(const Expr& e, AssumptionKey k) noexcept {
    if (!e) return std::nullopt;

    if (auto v = direct(e, k); v.has_value()) {
        return v;
    }

    // Derive from related predicates when the direct query is Unknown.
    switch (k) {
        case AssumptionKey::Real: {
            if (direct(e, AssumptionKey::Integer) == true) return true;
            if (direct(e, AssumptionKey::Rational) == true) return true;
            if (direct(e, AssumptionKey::Positive) == true) return true;
            if (direct(e, AssumptionKey::Negative) == true) return true;
            if (direct(e, AssumptionKey::Zero) == true) return true;
            if (direct(e, AssumptionKey::Even) == true) return true;
            if (direct(e, AssumptionKey::Odd) == true) return true;
            if (direct(e, AssumptionKey::Prime) == true) return true;
            if (direct(e, AssumptionKey::Composite) == true) return true;
            if (direct(e, AssumptionKey::Irrational) == true) return true;
            // A (nonzero) imaginary value is not real.
            if (direct(e, AssumptionKey::Imaginary) == true) return false;
            return std::nullopt;
        }
        case AssumptionKey::Rational: {
            if (direct(e, AssumptionKey::Integer) == true) return true;
            if (direct(e, AssumptionKey::Irrational) == true) return false;
            return std::nullopt;
        }
        case AssumptionKey::Integer: {
            // even / odd / prime / composite ⇒ integer
            if (direct(e, AssumptionKey::Even) == true) return true;
            if (direct(e, AssumptionKey::Odd) == true) return true;
            if (direct(e, AssumptionKey::Prime) == true) return true;
            if (direct(e, AssumptionKey::Composite) == true) return true;
            return std::nullopt;
        }
        case AssumptionKey::Nonzero: {
            if (direct(e, AssumptionKey::Positive) == true) return true;
            if (direct(e, AssumptionKey::Negative) == true) return true;
            if (direct(e, AssumptionKey::Zero) == false) return true;
            if (direct(e, AssumptionKey::Odd) == true) return true;  // odd ⇒ ≠ 0
            if (direct(e, AssumptionKey::Prime) == true) return true;  // prime ≥ 2
            if (direct(e, AssumptionKey::Composite) == true) return true;  // ≥ 4
            return std::nullopt;
        }
        case AssumptionKey::Prime: {
            // A non-integer is never prime; a composite is not prime; otherwise
            // primality is a direct fact.
            if (direct(e, AssumptionKey::Integer) == false) return false;
            if (direct(e, AssumptionKey::Composite) == true) return false;
            return std::nullopt;
        }
        case AssumptionKey::Composite: {
            // A non-integer is never composite; a prime is not composite.
            if (direct(e, AssumptionKey::Integer) == false) return false;
            if (direct(e, AssumptionKey::Prime) == true) return false;
            return std::nullopt;
        }
        case AssumptionKey::Irrational: {
            // irrational ⟺ real ∧ ¬rational. A rational or non-real value is not
            // irrational; a known real non-rational is irrational.
            if (direct(e, AssumptionKey::Rational) == true) return false;
            if (direct(e, AssumptionKey::Real) == false) return false;
            if (direct(e, AssumptionKey::Real) == true
                && direct(e, AssumptionKey::Rational) == false) {
                return true;
            }
            return std::nullopt;
        }
        case AssumptionKey::Algebraic: {
            // rational ⇒ algebraic; transcendental ⇒ ¬algebraic; within ℂ a
            // non-transcendental value is algebraic.
            if (direct(e, AssumptionKey::Rational) == true) return true;
            if (direct(e, AssumptionKey::Transcendental) == true) return false;
            if (direct(e, AssumptionKey::Complex) == true
                && direct(e, AssumptionKey::Transcendental) == false) {
                return true;
            }
            return std::nullopt;
        }
        case AssumptionKey::Transcendental: {
            // transcendental ⟺ complex ∧ ¬algebraic. A rational / algebraic value
            // is not transcendental.
            if (direct(e, AssumptionKey::Algebraic) == true) return false;
            if (direct(e, AssumptionKey::Rational) == true) return false;
            if (direct(e, AssumptionKey::Complex) == true
                && direct(e, AssumptionKey::Algebraic) == false) {
                return true;
            }
            return std::nullopt;
        }
        case AssumptionKey::Nonnegative: {
            if (direct(e, AssumptionKey::Positive) == true) return true;
            if (direct(e, AssumptionKey::Zero) == true) return true;
            if (direct(e, AssumptionKey::Negative) == true) return false;
            return std::nullopt;
        }
        case AssumptionKey::Nonpositive: {
            if (direct(e, AssumptionKey::Negative) == true) return true;
            if (direct(e, AssumptionKey::Zero) == true) return true;
            if (direct(e, AssumptionKey::Positive) == true) return false;
            return std::nullopt;
        }
        case AssumptionKey::Finite: {
            // infinite ⟺ ¬finite; otherwise numeric / structural overrides
            // handle the definite cases.
            if (direct(e, AssumptionKey::Infinite) == true) return false;
            return std::nullopt;
        }
        case AssumptionKey::ExtendedReal: {
            // real ⇒ extended_real; a nonzero pure imaginary is off the line.
            if (direct(e, AssumptionKey::Real) == true) return true;
            if (direct(e, AssumptionKey::Imaginary) == true) return false;
            return std::nullopt;
        }
        case AssumptionKey::Infinite: {
            // infinite ⟺ ¬finite.
            if (direct(e, AssumptionKey::Finite) == true) return false;
            if (direct(e, AssumptionKey::Finite) == false) return true;
            return std::nullopt;
        }
        case AssumptionKey::Complex: {
            // real ∨ imaginary ⇒ complex. Covers every real number and symbol
            // without each node having to answer Complex directly.
            if (direct(e, AssumptionKey::Real) == true) return true;
            if (direct(e, AssumptionKey::Imaginary) == true) return true;
            if (direct(e, AssumptionKey::Integer) == true) return true;
            if (direct(e, AssumptionKey::Rational) == true) return true;
            if (direct(e, AssumptionKey::Algebraic) == true) return true;
            if (direct(e, AssumptionKey::Transcendental) == true) return true;
            return std::nullopt;
        }
        case AssumptionKey::Imaginary: {
            // A real value is not imaginary; 0 is real, not imaginary.
            if (direct(e, AssumptionKey::Real) == true) return false;
            if (direct(e, AssumptionKey::Zero) == true) return false;
            return std::nullopt;
        }
        default:
            return std::nullopt;
    }
}

bool is_provably_odd(const Expr& e) noexcept {
    if (!e) return false;
    if (ask(e, AssumptionKey::Odd) == true) return true;  // declared/derived odd
    if (e->type_id() == TypeId::Integer) {
        return mpz_odd_p(static_cast<const Integer&>(*e).value().get_mpz_t()) != 0;
    }
    if (e->type_id() == TypeId::Mul) {
        for (const auto& f : e->args()) {
            if (!is_provably_odd(f)) return false;  // product of odds is odd
        }
        return true;
    }
    if (e->type_id() == TypeId::Add) {
        int odd_terms = 0;
        for (const auto& f : e->args()) {
            if (is_provably_even(f)) continue;
            if (is_provably_odd(f)) { ++odd_terms; continue; }
            return false;  // a term of unknown parity
        }
        return (odd_terms % 2) == 1;
    }
    return false;
}

bool is_provably_even(const Expr& e) noexcept {
    if (!e) return false;
    if (ask(e, AssumptionKey::Even) == true) return true;  // declared/derived even
    if (e->type_id() == TypeId::Integer) {
        return mpz_even_p(static_cast<const Integer&>(*e).value().get_mpz_t()) != 0;
    }
    if (e->type_id() == TypeId::Mul) {
        // Product of integers with at least one even factor → even.
        bool all_integer = true;
        bool has_even = false;
        for (const auto& f : e->args()) {
            if (is_integer(f) != true) all_integer = false;
            if (is_provably_even(f)) has_even = true;
        }
        return all_integer && has_even;
    }
    if (e->type_id() == TypeId::Add) {
        for (const auto& f : e->args()) {
            if (!is_provably_even(f)) return false;  // sum of evens is even
        }
        return true;
    }
    return false;
}

}  // namespace sympp
