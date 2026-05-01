#include <sympp/simplify/simplify.hpp>

#include <algorithm>
#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/boolean.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/expr_collections.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/piecewise.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/polys/poly.hpp>

namespace sympp {

namespace {

// Re-run the canonical-form factories on every node — many substitution
// chains leave nominally-canonical Expr that benefit from a fresh sweep.
[[nodiscard]] Expr re_canonicalize(const Expr& e) {
    auto args = e->args();
    if (args.empty()) return e;
    std::vector<Expr> new_args;
    new_args.reserve(args.size());
    for (const auto& a : args) new_args.push_back(re_canonicalize(a));
    switch (e->type_id()) {
        case TypeId::Add: return add(std::move(new_args));
        case TypeId::Mul: return mul(std::move(new_args));
        case TypeId::Pow: return pow(new_args[0], new_args[1]);
        default: return e;
    }
}

}  // namespace

Expr simplify(const Expr& e) {
    if (!e) return e;
    // Sweep 1: canonical form.
    Expr current = re_canonicalize(e);
    // Sweep 2: expand to flush nested products.
    Expr expanded = expand(current);
    // Sweep 3: canonical again to collect like terms after expansion.
    return re_canonicalize(expanded);
}

Expr collect(const Expr& e, const Expr& var) {
    // Build a Poly in var; reconstitute as a sum of (coeff * var^k).
    Poly p(expand(e), var);
    return p.as_expr();
}

// ----- powsimp ---------------------------------------------------------------

namespace {

// Reconstruct an expression of the same shape as `original` with new args.
// Mirrors traversal.cpp::rebuild_with_args; duplicated here to avoid making
// that helper public.
[[nodiscard]] Expr rebuild_node(const Expr& original,
                                std::vector<Expr> new_args) {
    switch (original->type_id()) {
        case TypeId::Add: return add(std::move(new_args));
        case TypeId::Mul: return mul(std::move(new_args));
        case TypeId::Pow: return pow(new_args[0], new_args[1]);
        case TypeId::Function: {
            const auto& f = static_cast<const Function&>(*original);
            return f.rebuild(std::move(new_args));
        }
        default: return original;
    }
}

// True iff combining bases under a shared non-integer exponent is safe.
// Conservative: requires every base to be definitely-positive.
[[nodiscard]] bool can_combine_under_exp(const std::vector<Expr>& bases,
                                          const Expr& exp) {
    if (exp->type_id() == TypeId::Integer) return true;
    for (const auto& b : bases) {
        auto p = is_positive(b);
        if (!p.has_value() || !*p) return false;
    }
    return true;
}

[[nodiscard]] Expr powsimp_node(const Expr& e) {
    if (e->type_id() != TypeId::Mul) return e;

    // Group factors by their exponent (as_base_exp).
    struct Group { Expr exp; std::vector<Expr> bases; };
    std::vector<Group> groups;
    auto find_or_create = [&](const Expr& exp) -> Group& {
        for (auto& g : groups) if (g.exp == exp) return g;
        groups.push_back(Group{exp, {}});
        return groups.back();
    };

    for (const auto& a : e->args()) {
        Expr base, exp;
        if (a->type_id() == TypeId::Pow) {
            base = a->args()[0];
            exp = a->args()[1];
        } else {
            base = a;
            exp = S::One();
        }
        find_or_create(exp).bases.push_back(base);
    }

    // Reassemble — merge groups with ≥ 2 bases when safe.
    std::vector<Expr> out;
    out.reserve(groups.size());
    for (auto& g : groups) {
        if (g.bases.size() == 1) {
            out.push_back(pow(g.bases[0], g.exp));
        } else if (can_combine_under_exp(g.bases, g.exp)) {
            out.push_back(pow(mul(g.bases), g.exp));
        } else {
            for (auto& b : g.bases) out.push_back(pow(b, g.exp));
        }
    }
    return mul(std::move(out));
}

}  // namespace

Expr powsimp(const Expr& e) {
    if (!e) return e;
    auto args = e->args();
    if (args.empty()) return e;
    std::vector<Expr> new_args;
    new_args.reserve(args.size());
    for (const auto& a : args) new_args.push_back(powsimp(a));
    Expr rebuilt = rebuild_node(e, std::move(new_args));
    return powsimp_node(rebuilt);
}

// ----- expand_power_exp / expand_power_base ----------------------------------

namespace {

[[nodiscard]] Expr expand_power_exp_node(const Expr& e) {
    if (e->type_id() != TypeId::Pow) return e;
    const Expr& base = e->args()[0];
    const Expr& exp = e->args()[1];
    if (exp->type_id() != TypeId::Add) return e;
    // x^(a + b + c) → x^a * x^b * x^c
    std::vector<Expr> factors;
    factors.reserve(exp->args().size());
    for (const auto& term : exp->args()) {
        factors.push_back(pow(base, term));
    }
    return mul(std::move(factors));
}

[[nodiscard]] Expr expand_power_base_node(const Expr& e) {
    if (e->type_id() != TypeId::Pow) return e;
    const Expr& base = e->args()[0];
    const Expr& exp = e->args()[1];
    if (base->type_id() != TypeId::Mul) return e;
    // Safe to distribute when exp is integer, or every factor of base is
    // positive.
    bool safe = (exp->type_id() == TypeId::Integer);
    if (!safe) {
        bool all_pos = true;
        for (const auto& f : base->args()) {
            auto p = is_positive(f);
            if (!p.has_value() || !*p) { all_pos = false; break; }
        }
        safe = all_pos;
    }
    if (!safe) return e;
    std::vector<Expr> factors;
    factors.reserve(base->args().size());
    for (const auto& f : base->args()) factors.push_back(pow(f, exp));
    return mul(std::move(factors));
}

template <typename NodeFn>
[[nodiscard]] Expr apply_recursive(const Expr& e, NodeFn fn) {
    if (!e) return e;
    auto args = e->args();
    if (args.empty()) return fn(e);
    std::vector<Expr> new_args;
    new_args.reserve(args.size());
    for (const auto& a : args) new_args.push_back(apply_recursive(a, fn));
    return fn(rebuild_node(e, std::move(new_args)));
}

}  // namespace

Expr expand_power_exp(const Expr& e) {
    return apply_recursive(e, expand_power_exp_node);
}

Expr expand_power_base(const Expr& e) {
    return apply_recursive(e, expand_power_base_node);
}

// ----- trigsimp --------------------------------------------------------------

namespace {

// If `e` is sin(arg)^2 or cos(arg)^2, return (is_sin, arg). Else nullopt.
[[nodiscard]] std::optional<std::pair<bool, Expr>>
detect_trig_pow2(const Expr& e) {
    if (e->type_id() != TypeId::Pow) return std::nullopt;
    if (!(e->args()[1] == integer(2))) return std::nullopt;
    const Expr& base = e->args()[0];
    if (base->type_id() != TypeId::Function) return std::nullopt;
    const auto& f = static_cast<const Function&>(*base);
    if (f.function_id() == FunctionId::Sin) return std::pair{true, base->args()[0]};
    if (f.function_id() == FunctionId::Cos) return std::pair{false, base->args()[0]};
    return std::nullopt;
}

// Decompose a term as (coef * sin(arg)^2) or (coef * cos(arg)^2). Returns
// nullopt when the term has neither, or has more than one such factor (an
// unusual structure we leave alone).
struct TrigSquareTerm {
    bool is_sin;
    Expr arg;
    Expr coef;
};
[[nodiscard]] std::optional<TrigSquareTerm>
as_trig_square_term(const Expr& term) {
    if (auto r = detect_trig_pow2(term); r) {
        return TrigSquareTerm{r->first, r->second, S::One()};
    }
    if (term->type_id() == TypeId::Mul) {
        std::optional<std::pair<bool, Expr>> trig;
        std::vector<Expr> coef_factors;
        for (const auto& f : term->args()) {
            if (auto r = detect_trig_pow2(f); r) {
                if (trig) return std::nullopt;  // multiple — leave alone
                trig = r;
            } else {
                coef_factors.push_back(f);
            }
        }
        if (trig) {
            return TrigSquareTerm{trig->first, trig->second, mul(coef_factors)};
        }
    }
    return std::nullopt;
}

[[nodiscard]] Expr trigsimp_add(const Expr& e) {
    if (e->type_id() != TypeId::Add) return e;

    struct CoefPair { Expr sin_coef = S::Zero(); Expr cos_coef = S::Zero(); };
    std::vector<std::pair<Expr, CoefPair>> by_arg;
    std::vector<Expr> non_trig;

    auto find_or_create = [&](const Expr& arg) -> CoefPair& {
        for (auto& [a, p] : by_arg) if (a == arg) return p;
        by_arg.push_back({arg, CoefPair{}});
        return by_arg.back().second;
    };

    for (const auto& term : e->args()) {
        auto t = as_trig_square_term(term);
        if (!t) { non_trig.push_back(term); continue; }
        auto& cp = find_or_create(t->arg);
        if (t->is_sin) cp.sin_coef = add(cp.sin_coef, t->coef);
        else cp.cos_coef = add(cp.cos_coef, t->coef);
    }

    std::vector<Expr> out = std::move(non_trig);
    for (auto& [arg, cp] : by_arg) {
        // Rewrite a*sin²(x) + b*cos²(x) as b + (a-b)*sin²(x). When a == b
        // this collapses to b — the Pythagorean identity. Otherwise it
        // produces an equivalent form with one fewer trig² (b absorbed).
        out.push_back(cp.cos_coef);
        Expr diff = cp.sin_coef - cp.cos_coef;
        if (!(diff == S::Zero())) {
            out.push_back(mul(diff, pow(sin(arg), integer(2))));
        }
    }
    if (out.empty()) return S::Zero();
    if (out.size() == 1) return out[0];
    return add(std::move(out));
}

}  // namespace

Expr trigsimp(const Expr& e) {
    return apply_recursive(e, trigsimp_add);
}

// ----- radsimp ---------------------------------------------------------------

namespace {

[[nodiscard]] bool is_sqrt_node(const Expr& e) {
    return e->type_id() == TypeId::Pow && e->args()[1] == S::Half();
}

// Decompose an expression `den` into rational_part + Σ coef_i * sqrt(arg_i).
// Walks an Add, splits each term into a sqrt-bearing form or rational form.
struct SqrtDecomp {
    std::vector<std::pair<Expr, Expr>> sqrt_terms;  // (coef, sqrt_arg)
    Expr rational_part = S::Zero();
};

[[nodiscard]] SqrtDecomp decompose_sqrts(const Expr& den) {
    SqrtDecomp result;
    auto handle_term = [&](const Expr& term) {
        Expr sqrt_factor;
        std::vector<Expr> coef_factors;
        if (term->type_id() == TypeId::Mul) {
            for (const auto& f : term->args()) {
                if (is_sqrt_node(f) && !sqrt_factor) {
                    sqrt_factor = f;
                } else {
                    coef_factors.push_back(f);
                }
            }
        } else if (is_sqrt_node(term)) {
            sqrt_factor = term;
        }
        if (sqrt_factor) {
            Expr coef = mul(coef_factors);
            result.sqrt_terms.emplace_back(coef, sqrt_factor->args()[0]);
        } else {
            result.rational_part = add(result.rational_part, term);
        }
    };
    if (den->type_id() == TypeId::Add) {
        for (const auto& t : den->args()) handle_term(t);
    } else {
        handle_term(den);
    }
    return result;
}

}  // namespace

Expr radsimp(const Expr& e) {
    // Combine into a single fraction so we have a clear num/denom.
    Expr t = together(e);
    if (t->type_id() != TypeId::Mul) return t;

    // Locate Pow(_, -1) factor — the denominator.
    Expr den;
    std::vector<Expr> num_factors;
    for (const auto& f : t->args()) {
        if (f->type_id() == TypeId::Pow
            && f->args()[1] == S::NegativeOne()
            && !den) {
            den = f->args()[0];
        } else {
            num_factors.push_back(f);
        }
    }
    if (!den) return t;
    Expr num = mul(num_factors);

    auto decomp = decompose_sqrts(den);
    if (decomp.sqrt_terms.size() != 1) return t;
    // Binomial: den = a + b * sqrt(c). Conjugate: a - b * sqrt(c).
    const auto& [b_coef, sqrt_arg] = decomp.sqrt_terms[0];
    const Expr& a = decomp.rational_part;
    Expr b_sqrt = mul(b_coef, sqrt(sqrt_arg));
    Expr conjugate = a - b_sqrt;
    Expr new_num = num * conjugate;
    // (a + b√c)(a - b√c) = a² - b²c
    Expr new_den = a * a - b_coef * b_coef * sqrt_arg;
    if (new_den == S::Zero()) return t;
    return new_num / new_den;
}

// ----- sqrtdenest ------------------------------------------------------------
//
// sqrt(a + b*sqrt(c)) denests to sqrt(d) + sqrt(e) when (a² - b²c) is a
// non-negative perfect square (call it r). Then:
//   d = (a + r) / 2
//   e = (a - r) / 2
// And sqrt(a + b*sqrt(c)) = sqrt(d) + sqrt(e).
//
// We only attempt this when a, b, c are numeric Integer/Rational so we can
// compute r exactly via Integer::is_square / mpz_rootrem.

namespace {

[[nodiscard]] std::optional<mpq_class> as_mpq(const Expr& e) {
    if (e->type_id() == TypeId::Integer) {
        return mpq_class(static_cast<const Integer&>(*e).value());
    }
    if (e->type_id() == TypeId::Rational) {
        return static_cast<const Rational&>(*e).value();
    }
    return std::nullopt;
}

[[nodiscard]] Expr mpq_to_sympp(const mpq_class& q) {
    mpq_class r = q;
    r.canonicalize();
    if (r.get_den() == 1) return make<Integer>(r.get_num());
    return make<Rational>(r);
}

[[nodiscard]] std::optional<mpq_class> rational_sqrt(const mpq_class& q) {
    if (q < 0) return std::nullopt;
    mpz_class num = q.get_num();
    mpz_class den = q.get_den();
    mpz_class num_root, num_rem, den_root, den_rem;
    mpz_rootrem(num_root.get_mpz_t(), num_rem.get_mpz_t(),
                num.get_mpz_t(), 2);
    mpz_rootrem(den_root.get_mpz_t(), den_rem.get_mpz_t(),
                den.get_mpz_t(), 2);
    if (num_rem != 0 || den_rem != 0) return std::nullopt;
    return mpq_class(num_root, den_root);
}

[[nodiscard]] Expr sqrtdenest_node(const Expr& e) {
    if (!is_sqrt_node(e)) return e;
    const Expr& inner = e->args()[0];
    if (inner->type_id() != TypeId::Add) return e;
    auto decomp = decompose_sqrts(inner);
    if (decomp.sqrt_terms.size() != 1) return e;
    const auto& [b_expr, c_expr] = decomp.sqrt_terms[0];
    auto a_q = as_mpq(decomp.rational_part);
    auto b_q = as_mpq(b_expr);
    auto c_q = as_mpq(c_expr);
    if (!a_q || !b_q || !c_q) return e;

    mpq_class disc = (*a_q) * (*a_q) - (*b_q) * (*b_q) * (*c_q);
    auto r_q = rational_sqrt(disc);
    if (!r_q) return e;
    // r >= 0 holds since rational_sqrt requires non-negative input.
    mpq_class d_q = ((*a_q) + (*r_q)) / mpq_class(2);
    mpq_class f_q = ((*a_q) - (*r_q)) / mpq_class(2);
    if (d_q < 0 || f_q < 0) return e;
    Expr d_e = mpq_to_sympp(d_q);
    Expr f_e = mpq_to_sympp(f_q);
    // Sign of b_q decides whether result is sqrt(d) + sqrt(f) or sqrt(d) - sqrt(f).
    Expr result_pos = sqrt(d_e) + sqrt(f_e);
    if (*b_q >= 0) return result_pos;
    return sqrt(d_e) - sqrt(f_e);
}

}  // namespace

Expr sqrtdenest(const Expr& e) {
    return apply_recursive(e, sqrtdenest_node);
}

}  // namespace sympp
