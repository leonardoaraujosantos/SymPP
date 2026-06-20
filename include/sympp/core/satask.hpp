#pragma once

// Boolean / SAT-style `ask` over predicate propositions.
//
// SymPP's primary `ask(expr, key)` answers a *single* predicate three-valued.
// SymPy additionally lets you query arbitrary boolean *combinations* of
// predicates — `ask(Q.positive(x) | Q.negative(x))` — and resolves them with a
// small SAT engine (`satask`). This header provides the same: a `Proposition`
// built from predicate literals with ¬ / ∧ / ∨, and an `ask(Proposition)` that
//
//   1. evaluates it three-valued (Kleene) from the per-literal `ask`, and
//   2. when that is Unknown and every literal concerns the *same* expression,
//      upgrades the answer by refutation — a proposition is True when asserting
//      its negation produces a logically inconsistent assumption mask (and
//      False when asserting the proposition itself does).
//
// So `Q(x, Positive) || Q(x, Negative)` is True for a real, nonzero `x` even
// though neither disjunct is individually decidable — exactly SymPy's result.
//
// Reference: sympy/assumptions/ask.py — `ask`; sympy/assumptions/satask.py.

#include <cstdint>
#include <optional>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/core/assumption_key.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

// A three-valued proposition over predicate literals `ask(expr, key) == value`.
// Build literals with Q()/Qn() and combine with the !, && and || operators.
class SYMPP_EXPORT Proposition {
public:
    enum class Op : std::uint8_t { Lit, Not, And, Or };

    // Literal: asserts `ask(e, k) == expected`.
    Proposition(Expr e, AssumptionKey k, bool expected = true);

    [[nodiscard]] Op op() const noexcept { return op_; }
    [[nodiscard]] const Expr& expr() const noexcept { return expr_; }
    [[nodiscard]] AssumptionKey key() const noexcept { return key_; }
    [[nodiscard]] bool expected() const noexcept { return expected_; }
    [[nodiscard]] const std::vector<Proposition>& args() const noexcept { return args_; }

    // Connective factories (used by the operators below).
    [[nodiscard]] static Proposition make(Op op, std::vector<Proposition> args);

private:
    Proposition() = default;

    Op op_{Op::Lit};
    Expr expr_;
    AssumptionKey key_{};
    bool expected_{true};
    std::vector<Proposition> args_;
};

// Literal builders: Q(e, k) ≡ "ask(e, k) is true"; Qn(e, k) ≡ "... is false".
[[nodiscard]] SYMPP_EXPORT Proposition Q(const Expr& e, AssumptionKey k);
[[nodiscard]] SYMPP_EXPORT Proposition Qn(const Expr& e, AssumptionKey k);

[[nodiscard]] SYMPP_EXPORT Proposition operator!(const Proposition& p);
[[nodiscard]] SYMPP_EXPORT Proposition operator&&(const Proposition& a, const Proposition& b);
[[nodiscard]] SYMPP_EXPORT Proposition operator||(const Proposition& a, const Proposition& b);

// Resolve the proposition to true / false / unknown (std::nullopt).
[[nodiscard]] SYMPP_EXPORT std::optional<bool> ask(const Proposition& p) noexcept;

}  // namespace sympp
