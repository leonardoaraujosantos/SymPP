#include <sympp/simplify/hyperexpand.hpp>

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/hypergeometric.hpp>
#include <sympp/functions/trigonometric.hpp>

namespace sympp {

namespace {

[[nodiscard]] bool is_hyper(const Expr& e) {
    if (e->type_id() != TypeId::Function) return false;
    const auto& fn = static_cast<const Function&>(*e);
    return fn.function_id() == FunctionId::Hyper;
}

[[nodiscard]] bool both_match(const std::vector<Expr>& v,
                                  const std::vector<Expr>& target) {
    if (v.size() != target.size()) return false;
    // Order-independent match — Pochhammer products commute, so SymPy
    // also normalizes the parameter sets before checking.
    std::vector<bool> used(target.size(), false);
    for (const auto& a : v) {
        bool matched = false;
        for (std::size_t j = 0; j < target.size(); ++j) {
            if (used[j]) continue;
            if (a == target[j]) {
                used[j] = true;
                matched = true;
                break;
            }
        }
        if (!matched) return false;
    }
    return true;
}

// Rewrite a hyper(ap, bq, z) node using the recognized-identity table.
// Returns std::nullopt if no rule fires; the caller leaves the node alone.
[[nodiscard]] std::optional<Expr> rewrite_hyper_node(const Expr& e) {
    if (!is_hyper(e)) return std::nullopt;
    const auto& h = static_cast<const Hyper&>(*e);
    auto a = h.ap();
    auto b = h.bq();
    auto z = h.z();
    auto p = a.size();
    auto q = b.size();

    Expr half = rational(1, 2);
    Expr threehalf = rational(3, 2);
    Expr one = S::One();
    Expr two = integer(2);

    // ₀F₁ family (p == 0, q == 1).
    if (p == 0 && q == 1) {
        // ₀F₁(; 1/2; −z²/4) = cos(z) — match by setting w := -z*4 and
        // checking arg-shape. Simpler: introduce y so z == −y²/4 ↔ y² == −4z.
        // We don't currently invert that, so we recognize the canonical
        // input form ₀F₁(; 1/2; w) = cosh(√w) / 1 and ₀F₁(; 3/2; w) =
        // sinh(√w) / √w. The caller usually supplies the form below.
        // We don't ship those here to avoid Sqrt-with-branch surprises.
    }

    // ₁F₁(a; b; z) family (p == 1, q == 1).
    if (p == 1 && q == 1) {
        // ₁F₁(1; 2; z) = (eᶻ − 1) / z.
        if (a[0] == one && b[0] == two) {
            return (exp(z) - one) / z;
        }
    }

    // ₂F₁ family (p == 2, q == 1).
    if (p == 2 && q == 1) {
        // ₂F₁(1, 1; 2; −z) = log(1 + z)/z      (so input z↦z')
        // ₂F₁(1, 1; 2; z')  = −log(1 − z')/z'.
        if (both_match(a, {one, one}) && b[0] == two) {
            return -log(one - z) / z;
        }
        // ₂F₁(1/2, 1; 3/2; w) = arctan(√w) / √w  on the "w = z²" reading
        // is what users want for compact closed forms. We instead
        // recognize the form most callers produce, which is z² already
        // wrapped: ₂F₁(1/2, 1; 3/2; -z²) = arctan(z)/z, but only when
        // z is recognizable as -y²; that requires structural unwrapping
        // we skip here. A safer rewrite that always fires:
        //
        //   ₂F₁(1/2, 1; 3/2; w) = arctan(√w) / √w
        //
        // is left for a follow-up.
        if (both_match(a, {half, one}) && b[0] == threehalf) {
            // No safe closed form without sqrt-of-arg gymnastics.
        }
        // ₂F₁(1/2, 1/2; 3/2; w) = arcsin(√w) / √w — same caveat.
    }

    return std::nullopt;
}

[[nodiscard]] Expr apply_recursive(const Expr& e) {
    if (!e) return e;
    auto args = e->args();
    if (args.empty()) {
        if (auto r = rewrite_hyper_node(e); r) return *r;
        return e;
    }
    std::vector<Expr> new_args;
    new_args.reserve(args.size());
    for (const auto& a : args) new_args.push_back(apply_recursive(a));
    Expr rebuilt;
    switch (e->type_id()) {
        case TypeId::Add: rebuilt = add(std::move(new_args)); break;
        case TypeId::Mul: rebuilt = mul(std::move(new_args)); break;
        case TypeId::Pow: rebuilt = pow(new_args[0], new_args[1]); break;
        case TypeId::Function:
            // Rebuild the function with rewritten args via its rebuild()
            // hook (which goes back through the auto-eval factories).
            rebuilt = static_cast<const Function&>(*e).rebuild(
                std::move(new_args));
            break;
        default:
            rebuilt = e;
            break;
    }
    if (auto r = rewrite_hyper_node(rebuilt); r) return *r;
    return rebuilt;
}

}  // namespace

Expr hyperexpand(const Expr& e) { return apply_recursive(e); }

}  // namespace sympp
