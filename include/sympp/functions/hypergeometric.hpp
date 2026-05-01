#pragma once

// Generalized hypergeometric function pFq and Meijer-G function.
//
// `Hyper((a₁,…,aₚ), (b₁,…,b_q), z)` is the generalized hypergeometric
//
//     pFq(a; b; z) = Σ_{n=0..∞} ((a₁)_n ⋯ (aₚ)_n / (b₁)_n ⋯ (b_q)_n) zⁿ / n!
//
// where (a)_n is the Pochhammer rising factorial. Common cases:
//
//   ₀F₀(z)     = exp(z)
//   ₁F₀(a; z)  = (1 − z)^(−a)
//   ₂F₁(a, b; c; z)   the Gauss hypergeometric
//   ₁F₁(a; b; z)      Kummer confluent
//   ₀F₁( ; b; z)      regular at the origin (Bessel-related)
//
// Encoding: SymPP doesn't have a Tuple node, so the variadic parameter
// lists are encoded position-style:
//
//     args = [Integer(p), Integer(q), a₁, …, aₚ, b₁, …, b_q, z]
//
// The factory `hyper(ap, bq, z)` packs lists into args and applies
// auto-eval rules (numerator/denominator cancellation, z=0, well-known
// closed forms). For deeper rewriting use `hyperexpand`.
//
// Meijer-G is encoded similarly with four parameter lists:
//
//     G^{m,n}_{p,q}((a₁..aₙ), (a_{n+1}..aₚ), (b₁..b_m), (b_{m+1}..b_q), z)
//
//     args = [Integer(n_first), Integer(p_rest), Integer(m_first),
//             Integer(q_rest), <flattened sequences>, z]
//
// SymPP currently keeps Meijer-G mostly opaque — full Slater-theorem
// evaluation and Mellin-Barnes back-transform are deferred-deep.
//
// Reference: sympy/functions/special/hyper.py, sympy/simplify/hyperexpand.py

#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class SYMPP_EXPORT Hyper final : public Function {
public:
    explicit Hyper(std::vector<Expr> args);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::Hyper;
    }
    [[nodiscard]] std::string_view name() const noexcept override { return "hyper"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;

    // Structural accessors. p = |numerator params|, q = |denominator params|.
    [[nodiscard]] std::size_t p() const;
    [[nodiscard]] std::size_t q() const;
    [[nodiscard]] std::vector<Expr> ap() const;
    [[nodiscard]] std::vector<Expr> bq() const;
    [[nodiscard]] Expr z() const;
};

class SYMPP_EXPORT MeijerG final : public Function {
public:
    explicit MeijerG(std::vector<Expr> args);
    [[nodiscard]] FunctionId function_id() const noexcept override {
        return FunctionId::MeijerG;
    }
    [[nodiscard]] std::string_view name() const noexcept override { return "meijerg"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::string str() const override;
};

// hyper(ap, bq, z) — factory with auto-eval. Recognized closed forms:
//   * z == 0                                            → 1
//   * any aᵢ ∈ ap matches any bⱼ ∈ bq                  → cancel both
//     (after cancellation, re-runs the factory)
//   * after cancellation: p == 0 ∧ q == 0              → exp(z)
//   * after cancellation: p == 1 ∧ q == 0              → (1 − z)^(−a₁)
//
// Other identities (₂F₁ → log/arctan/arcsin/etc.) are recognized by
// `hyperexpand`, not at construction time, since they are not always
// the "preferred" form.
[[nodiscard]] SYMPP_EXPORT Expr hyper(const std::vector<Expr>& ap,
                                         const std::vector<Expr>& bq,
                                         const Expr& z);

// hyper(a, b, c, z) — convenience for the ₂F₁ Gauss case. Builds the
// variadic Hyper with ap=[a, b], bq=[c]. Kept as a 4-arg overload so
// the existing dsolve_hypergeometric API continues to work.
[[nodiscard]] SYMPP_EXPORT Expr hyper(const Expr& a, const Expr& b,
                                         const Expr& c, const Expr& z);

// meijerg(an, ap_rest, bm, bq_rest, z) — opaque factory. No auto-eval.
[[nodiscard]] SYMPP_EXPORT Expr meijerg(const std::vector<Expr>& an,
                                           const std::vector<Expr>& ap_rest,
                                           const std::vector<Expr>& bm,
                                           const std::vector<Expr>& bq_rest,
                                           const Expr& z);

}  // namespace sympp
