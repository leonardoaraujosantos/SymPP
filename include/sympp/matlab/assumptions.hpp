#pragma once

// MATLAB facade — assumptions.
//
// MATLAB lets you write
//
//     syms x positive
//     assume(x, 'positive')
//     assumeAlso(x, 'integer')
//
// and have those facts attached to `x` mutably. SymPP `Symbol` is
// immutable: an assumption mask is fixed at construction. The two
// models can't be reconciled exactly without rewiring the hash-cons
// cache (which would break the immutability invariants
// downstream code relies on).
//
// Honest port: a process-wide registry mapping `name → AssumptionMask`.
// `assume(x, "positive")` registers the fact under x's printed name.
// `refresh(x)` returns a fresh Symbol carrying the registered mask.
// The original `x` is unchanged — the caller must rebind:
//
//     auto x = sym("x");
//     assume(x, "positive");
//     x = refresh(x);              // now x carries the mask
//     simplify(sqrt(x*x));         // uses positivity → x
//
// Property strings recognized: "real", "rational", "integer",
// "positive", "negative", "zero", "nonzero", "nonnegative",
// "nonpositive", "finite". Anything else (e.g. MATLAB's "even",
// "odd", "prime") throws std::runtime_error — SymPP's mask doesn't
// represent those yet (Phase 2 deep-deferred).
//
// Reference: sympy/core/assumptions.py — assumption fact propagation.

#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/core/assumption_key.hpp>
#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/fwd.hpp>

namespace sympp::matlab {

namespace detail {

// Process-wide registry. We keep this inline-defined as a
// function-local static so consumers don't have to link against any
// extra .cpp.
[[nodiscard]] inline std::unordered_map<std::string, AssumptionMask>&
assumption_registry() {
    static std::unordered_map<std::string, AssumptionMask> reg;
    return reg;
}

// Map MATLAB-style property strings onto AssumptionKey.
[[nodiscard]] inline AssumptionKey property_key(std::string_view prop) {
    if (prop == "real")        return AssumptionKey::Real;
    if (prop == "rational")    return AssumptionKey::Rational;
    if (prop == "integer")     return AssumptionKey::Integer;
    if (prop == "positive")    return AssumptionKey::Positive;
    if (prop == "negative")    return AssumptionKey::Negative;
    if (prop == "zero")        return AssumptionKey::Zero;
    if (prop == "nonzero")     return AssumptionKey::Nonzero;
    if (prop == "nonnegative") return AssumptionKey::Nonnegative;
    if (prop == "nonpositive") return AssumptionKey::Nonpositive;
    if (prop == "finite")      return AssumptionKey::Finite;
    throw std::runtime_error(
        std::string("assumption '") + std::string(prop)
        + "' not yet supported in SymPP");
}

[[nodiscard]] inline std::string_view key_name(AssumptionKey k) noexcept {
    switch (k) {
        case AssumptionKey::Real:        return "real";
        case AssumptionKey::Rational:    return "rational";
        case AssumptionKey::Integer:     return "integer";
        case AssumptionKey::Positive:    return "positive";
        case AssumptionKey::Negative:    return "negative";
        case AssumptionKey::Zero:        return "zero";
        case AssumptionKey::Nonzero:     return "nonzero";
        case AssumptionKey::Nonnegative: return "nonnegative";
        case AssumptionKey::Nonpositive: return "nonpositive";
        case AssumptionKey::Finite:      return "finite";
    }
    return "unknown";
}

}  // namespace detail

// assume(x, "positive") — replace any existing assumption mask for
// x's name with one carrying the single property.
inline void assume(const Expr& sym_expr, std::string_view property) {
    AssumptionKey k = detail::property_key(property);
    AssumptionMask m;
    m.set(k, true);
    detail::assumption_registry()[sym_expr->str()] = close_assumptions(m);
}

// assumeAlso(x, "integer") — merge the property into the existing
// mask for x's name (without dropping previously-set facts).
inline void assumeAlso(const Expr& sym_expr, std::string_view property) {
    AssumptionKey k = detail::property_key(property);
    auto& reg = detail::assumption_registry();
    auto& m = reg[sym_expr->str()];     // default-insert if missing
    m.set(k, true);
    m = close_assumptions(m);
}

// assumptions(x) — list of property strings registered for x.
[[nodiscard]] inline std::vector<std::string> assumptions(
    const Expr& sym_expr) {
    auto& reg = detail::assumption_registry();
    auto it = reg.find(sym_expr->str());
    if (it == reg.end()) return {};
    const AssumptionMask& m = it->second;
    std::vector<std::string> out;
    for (auto k : {AssumptionKey::Real, AssumptionKey::Rational,
                    AssumptionKey::Integer, AssumptionKey::Positive,
                    AssumptionKey::Negative, AssumptionKey::Zero,
                    AssumptionKey::Nonzero, AssumptionKey::Nonnegative,
                    AssumptionKey::Nonpositive, AssumptionKey::Finite}) {
        auto v = m.get(k);
        if (v == std::optional<bool>{true}) {
            out.emplace_back(detail::key_name(k));
        }
    }
    return out;
}

// clearAssumptions(x) — remove the registered mask for x's name.
inline void clearAssumptions(const Expr& sym_expr) {
    detail::assumption_registry().erase(sym_expr->str());
}

// refresh(x) — return a fresh Symbol named after x carrying the
// registered mask. If no mask is registered, returns the input
// expression unchanged (only meaningful for Symbol inputs).
[[nodiscard]] inline Expr refresh(const Expr& sym_expr) {
    auto& reg = detail::assumption_registry();
    auto it = reg.find(sym_expr->str());
    if (it == reg.end()) return sym_expr;
    return symbol(sym_expr->str(), it->second);
}

}  // namespace sympp::matlab
