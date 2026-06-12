#include <sympp/core/queries.hpp>

#include <optional>

#include <gmpxx.h>

#include <sympp/core/assumption_key.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

// Convenience: ask the underlying node directly (no implication-chasing) so
// we don't loop through this very function while computing derivations.
[[nodiscard]] inline std::optional<bool> direct(const Expr& e, AssumptionKey k) noexcept {
    return e ? e->ask(k) : std::nullopt;
}

}  // namespace

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
            return std::nullopt;
        }
        case AssumptionKey::Rational: {
            if (direct(e, AssumptionKey::Integer) == true) return true;
            return std::nullopt;
        }
        case AssumptionKey::Nonzero: {
            if (direct(e, AssumptionKey::Positive) == true) return true;
            if (direct(e, AssumptionKey::Negative) == true) return true;
            if (direct(e, AssumptionKey::Zero) == false) return true;
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
            // No simple universal implications here; numeric / structural
            // overrides handle the definite cases.
            return std::nullopt;
        }
        default:
            return std::nullopt;
    }
}

bool is_provably_odd(const Expr& e) noexcept {
    if (!e) return false;
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
