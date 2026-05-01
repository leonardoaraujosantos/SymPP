#include <sympp/core/queries.hpp>

#include <optional>

#include <sympp/core/assumption_key.hpp>
#include <sympp/core/basic.hpp>

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

}  // namespace sympp
