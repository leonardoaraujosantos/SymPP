#include <sympp/simplify/simplify.hpp>

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/boolean.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/expr_collections.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/piecewise.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/combinatorial.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/polys/poly.hpp>
#include <sympp/simplify/hyperexpand.hpp>

namespace sympp {

namespace {

// Total node count — used as a cheap "is this simpler?" measure, mirroring
// SymPy's count_ops-based form selection inside simplify().
[[nodiscard]] int node_count(const Expr& e) {
    int n = 0;
    auto rec = [&](auto&& self, const Expr& x) -> void {
        if (!x) return;
        ++n;
        for (const auto& a : x->args()) self(self, a);
    };
    rec(rec, e);
    return n;
}

// Does `e` contain a √ (a Pow with exponent 1/2) anywhere?
[[nodiscard]] bool contains_sqrt(const Expr& e) {
    if (e->type_id() == TypeId::Pow && e->args()[1] == S::Half()) return true;
    for (const auto& a : e->args()) {
        if (contains_sqrt(a)) return true;
    }
    return false;
}

// True iff `e` has an un-rationalized surd denominator — a Pow with a negative
// exponent whose base carries a √ (e.g. the 1/(1+√2) in Pow(1+√2, −1)).
// radsimp removes these; a rationalized result is preferred even when its node
// count is larger, so simplify()'s anti-bloat guard exempts this case.
[[nodiscard]] bool has_surd_denominator(const Expr& e) {
    if (e->type_id() == TypeId::Pow) {
        const Expr& ex = e->args()[1];
        const bool neg =
            (ex->type_id() == TypeId::Integer
             && static_cast<const Integer&>(*ex).is_negative())
            || (ex->type_id() == TypeId::Rational
                && static_cast<const Rational&>(*ex).value() < 0);
        if (neg && contains_sqrt(e->args()[0])) return true;
    }
    for (const auto& a : e->args()) {
        if (has_surd_denominator(a)) return true;
    }
    return false;
}

// True if e contains a complex denominator — Pow(d, negative) whose base d
// carries the imaginary unit. A rationalized result (no such denominator) is
// preferred even when larger, so the anti-bloat guard exempts it.
[[nodiscard]] bool has_complex_denominator(const Expr& e) {
    if (e->type_id() == TypeId::Pow
        && e->args()[1]->type_id() == TypeId::Integer
        && static_cast<const Integer&>(*e->args()[1]).is_negative()
        && has(e->args()[0], S::I())) {
        return true;
    }
    for (const auto& a : e->args()) {
        if (has_complex_denominator(a)) return true;
    }
    return false;
}

// True iff `e` is a rational function in its symbols — built only from
// numbers, symbols, +, *, and integer powers. cancel()/Poly() can loop
// forever on transcendental subexpressions (e.g. sin(x)), so simplify()
// must not hand those to the cancel step.
[[nodiscard]] bool is_rational_function(const Expr& e) {
    switch (e->type_id()) {
        case TypeId::Integer:
        case TypeId::Rational:
        case TypeId::Symbol:
            return true;
        case TypeId::Add:
        case TypeId::Mul:
            for (const auto& a : e->args()) {
                if (!is_rational_function(a)) return false;
            }
            return true;
        case TypeId::Pow:
            // Polynomial/rational only when the exponent is a literal integer
            // (covers denominators via negative exponents like (x-1)**(-1)).
            return is_rational_function(e->args()[0])
                   && e->args()[1]->type_id() == TypeId::Integer;
        default:
            return false;
    }
}

// True iff `e` contains a genuine symbol-dependent denominator — a Pow with
// a negative integer exponent whose base involves a free symbol (e.g. the
// (x-1)**(-1) in (x**2-1)/(x-1)). A bare rational coefficient like 1/2 is a
// Rational node, not such a Pow, so it does not count. cancel() is only worth
// running — and only safe from looping — on expressions that actually have
// one of these to reduce.
[[nodiscard]] bool has_symbolic_denominator(const Expr& e) {
    if (e->type_id() == TypeId::Pow) {
        const Expr& exp = e->args()[1];
        if (exp->type_id() == TypeId::Integer
            && static_cast<const Integer&>(*exp).is_negative()
            && !free_symbols(e->args()[0]).empty()) {
            return true;
        }
    }
    for (const auto& a : e->args()) {
        if (has_symbolic_denominator(a)) return true;
    }
    return false;
}

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

// Forward declarations; defined below, used by the simplify pipeline.
[[nodiscard]] Expr pow_of_pow_node(const Expr& e);
[[nodiscard]] Expr pow_of_pow(const Expr& e);
[[nodiscard]] Expr combine_exp(const Expr& e);
[[nodiscard]] Expr exp_to_hyp_add(const Expr& e);
[[nodiscard]] Expr exp_log_sum(const Expr& e);
[[nodiscard]] Expr radical_coeff(const Expr& e);

}  // namespace

Expr simplify(const Expr& e) {
    if (!e) return e;
    // 1. Canonical form.
    Expr canon = re_canonicalize(e);
    // 2. Expand to flush nested products, then rationalize complex denominators
    //    ((1+I)/(1-I) → I, 1/(1+I) → 1/2 − I/2).
    Expr current = expand(rationalize_complex(expand(canon)));
    // 3. Apply pattern-based simplifiers. Each is a no-op when the input
    //    doesn't match its pattern, so chaining is safe; ordering only
    //    affects which form we land on when multiple rules could apply.
    current = trigsimp(current);
    current = powsimp(current);
    current = combine_exp(current);
    current = exp_to_hyp_add(current);
    current = exp_log_sum(current);
    current = pow_of_pow(current);
    current = combsimp(current);
    current = gammasimp(current);
    current = radsimp(current);
    current = sqrtdenest(current);
    current = hyperexpand(current);
    // 4. Canonical sweep again to collect any like terms that surfaced.
    current = re_canonicalize(current);

    // 5. Rational-function cancellation. The pattern pipeline does not reduce
    //    n/d by polynomial GCD, so simplify() could leave (or even inflate)
    //    cancellable fractions like (x**2-1)/(x-1). cancel() runs its own
    //    together()+GCD; adopt its result only when strictly simpler so
    //    simplify() never returns something larger than the pipeline already
    //    produced.
    //
    //    Restricted to *univariate* rational functions with a symbol-dependent
    //    denominator. cancel() no longer hangs on multivariate input
    //    (CANCEL-1 fixed: divmod expands coefficient subtractions), and is
    //    available directly via cancel(); but auto-applying it here for the
    //    multivariate case regressed a downstream ODE form, so simplify keeps
    //    the conservative univariate scope. (Multivariate cancel-in-simplify
    //    is a separate quality task.)
    auto syms = free_symbols(canon);
    if (syms.size() == 1 && is_rational_function(canon)
        && has_symbolic_denominator(canon)) {
        try {
            Expr cancelled = cancel(canon, *syms.begin());
            cancelled = re_canonicalize(cancelled);
            if (node_count(cancelled) < node_count(current)) {
                current = cancelled;
            }
        } catch (const std::exception&) {
            // cancel() rejected the form; keep the pipeline result rather
            // than propagating out of simplify().
        }
    }

    // 6. Global anti-bloat guard. simplify() must never return something
    //    structurally larger than its (canonical) input — SymPy's simplify
    //    makes the same guarantee via a complexity measure. The pattern
    //    pipeline expands eagerly, which inflates already-compact forms
    //    ((x+1)**3, (exp(x)-1)/(exp(x/2)+1), a/b + c/b); when the pipeline
    //    result ends up bigger, fall back to the canonical input.
    if (node_count(current) > node_count(canon)) {
        // Exception: a rationalized surd denominator is the preferred canonical
        // form even when larger. Only fall back to canon when current did not
        // remove a surd denominator that canon still carries.
        const bool removed_surd =
            has_surd_denominator(canon) && !has_surd_denominator(current);
        const bool removed_cx =
            has_complex_denominator(canon) && !has_complex_denominator(current);
        if (!removed_surd && !removed_cx) {
            // Pulling a perfect-power factor out of a radical (√(4a²) → 2√(a²)) is a
            // simplification even when it raises the node count, so apply it after
            // the anti-bloat guard rather than inside the pipeline (where the guard
            // would revert it).
            return radical_coeff(canon);
        }
    }
    return radical_coeff(current);
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

// √(4·a²) → 2·√(a²), √(8·a²) → 2·√(2·a²), (8·x)^(1/3) → 2·x^(1/3): pull the
// perfect-power factor of a positive numeric coefficient out of a radical (a Pow
// with a non-integer rational exponent over a Mul base). The coefficient's
// non-perfect part stays under the radical with the symbolic factors (matching
// SymPy up to the √c·√X ↔ √(c·X) regrouping). Valid since the coefficient is
// positive: (c·X)^e = c^e·X^e on the principal branch regardless of X's sign.
[[nodiscard]] Expr radical_coeff_node(const Expr& e) {
    if (e->type_id() != TypeId::Pow) return e;
    const Expr& base = e->args()[0];
    const Expr& exp = e->args()[1];
    if (exp->type_id() != TypeId::Rational) return e;  // non-integer rational only
    if (base->type_id() != TypeId::Mul) return e;
    Expr coeff = S::One();
    std::vector<Expr> rest;
    for (const auto& f : base->args()) {
        if (f->type_id() == TypeId::Integer || f->type_id() == TypeId::Rational) {
            coeff = mul(coeff, f);
        } else {
            rest.push_back(f);
        }
    }
    if (rest.empty() || coeff == S::One()) return e;       // pure number / no coeff
    if (is_positive(coeff) != std::optional<bool>{true}) return e;  // need c > 0

    // coeff^exp auto-factors to m·(c')^exp (e.g. 8^(1/2) = 2·2^(1/2), 4^(1/2) = 2).
    // Pull out the rational m and keep c' under the radical with `rest`.
    Expr ce = pow(coeff, exp);
    Expr m;
    Expr cprime = S::One();
    if (ce->type_id() == TypeId::Integer || ce->type_id() == TypeId::Rational) {
        m = ce;  // coeff was a perfect power
    } else if (ce->type_id() == TypeId::Mul) {
        std::vector<Expr> nums;
        std::vector<Expr> others;
        for (const auto& f : ce->args()) {
            if (f->type_id() == TypeId::Integer
                || f->type_id() == TypeId::Rational) {
                nums.push_back(f);
            } else {
                others.push_back(f);
            }
        }
        if (nums.empty() || others.size() != 1
            || others[0]->type_id() != TypeId::Pow
            || !(others[0]->args()[1] == exp)) {
            return e;
        }
        m = mul(nums);
        cprime = others[0]->args()[0];
    } else {
        return e;  // coeff^exp stayed a bare radical — nothing to pull
    }
    Expr inside = (cprime == S::One()) ? mul(rest) : mul(cprime, mul(rest));
    return mul(m, pow(inside, exp));
}

// exp(… + c·log(p) + …) → p^c · exp(rest) for positive p (any addend that is a
// log of a positive base, optionally scaled by a constant, is pulled out as a
// power). Mirrors SymPy's expand/simplify of exp over a log-bearing sum.
[[nodiscard]] Expr exp_log_sum_node(const Expr& e) {
    if (e->type_id() != TypeId::Function) return e;
    const auto& fn = static_cast<const Function&>(*e);
    if (fn.function_id() != FunctionId::Exp || fn.args().size() != 1) return e;
    const Expr& arg = fn.args()[0];
    if (arg->type_id() != TypeId::Add) return e;

    // From a term, extract (p, c) when the term is c·log(p) with p positive.
    auto as_log_power = [](const Expr& term) -> std::optional<std::pair<Expr, Expr>> {
        auto from_log = [](const Expr& g) -> std::optional<Expr> {
            if (g->type_id() == TypeId::Function) {
                const auto& gf = static_cast<const Function&>(*g);
                if (gf.function_id() == FunctionId::Log && gf.args().size() == 1
                    && is_positive(gf.args()[0]) == true) {
                    return gf.args()[0];
                }
            }
            return std::nullopt;
        };
        if (auto p = from_log(term)) return std::make_pair(*p, S::One());
        if (term->type_id() == TypeId::Mul) {
            Expr base;
            std::vector<Expr> coeff;
            for (const auto& f : term->args()) {
                if (!base) {
                    if (auto p = from_log(f)) { base = *p; continue; }
                }
                coeff.push_back(f);
            }
            if (base && !coeff.empty()) return std::make_pair(base, mul(coeff));
        }
        return std::nullopt;
    };

    std::vector<Expr> factors;
    std::vector<Expr> rest;
    for (const auto& term : arg->args()) {
        if (auto lp = as_log_power(term)) {
            factors.push_back(pow(lp->first, lp->second));
        } else {
            rest.push_back(term);
        }
    }
    if (factors.empty()) return e;
    Expr prod = mul(std::move(factors));
    Expr rest_exp = rest.empty() ? S::One() : exp(add(std::move(rest)));
    return mul(prod, rest_exp);
}

// Assumptions-gated power-of-power: (bᵖ)^q. SymPP's canonical Pow deliberately
// leaves these alone (branch cuts), but symbol assumptions make specific cases
// safe — matching SymPy's `sqrt(x**2) == Abs(x)` (real) / `== x` (nonnegative).
[[nodiscard]] Expr pow_of_pow_node(const Expr& e) {
    if (e->type_id() != TypeId::Pow) return e;
    const Expr& inner = e->args()[0];
    const Expr& q = e->args()[1];

    // |y|^(2m) = y^(2m) for real y (the only case where |y|² equals y²; for
    // complex y, |y|² = y·ȳ ≠ y²). q must be a positive even integer.
    if (inner->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*inner);
        if (fn.function_id() == FunctionId::Abs && fn.args().size() == 1
            && q->type_id() == TypeId::Integer
            && static_cast<const Integer&>(*q).is_positive()
            && static_cast<const Integer&>(*q).fits_long()
            && static_cast<const Integer&>(*q).to_long() % 2 == 0
            && is_real(fn.args()[0]) == true) {
            return pow(fn.args()[0], q);
        }
    }

    if (inner->type_id() != TypeId::Pow) return e;
    const Expr& base = inner->args()[0];
    const Expr& p = inner->args()[1];

    // b ≥ 0: (bᵖ)^q = b^(p·q) holds for all real p, q.
    if (is_nonnegative(base) == true) {
        return pow(base, mul(p, q));
    }
    // b real: √(b²) = |b|. (b not known nonnegative here.)
    if (is_real(base) == true && p == integer(2) && q == rational(1, 2)) {
        return abs(base);
    }
    return e;
}

// Combine exponential factors in a product: e^a · e^b → e^(a+b), and e^a · (e^b)^k
// → e^(a + k·b). The canonical Mul keeps exp(a)·exp(b) separate (exp is a
// Function, not Pow(E, ·)), so this matches SymPy's `simplify`/`powsimp`. e^0
// folds to 1, so e.g. e^x · e^(−x) → 1.
[[nodiscard]] Expr combine_exp_node(const Expr& e) {
    // (exp(g))^k → exp(k·g) for an INTEGER exponent k. A symbolic or fractional
    // exponent is left as a Pow (matching SymPy, which keeps exp(x)**n and
    // sqrt(exp(x)) for branch-cut safety). The Mul case below merges products;
    // this closes the standalone power.
    if (e->type_id() == TypeId::Pow) {
        const Expr& base = e->args()[0];
        const Expr& ex = e->args()[1];
        if (ex->type_id() == TypeId::Integer
            && base->type_id() == TypeId::Function) {
            const auto& bfn = static_cast<const Function&>(*base);
            if (bfn.function_id() == FunctionId::Exp && bfn.args().size() == 1) {
                return exp(expand(mul(ex, bfn.args()[0])));
            }
        }
        return e;
    }
    if (e->type_id() != TypeId::Mul) return e;
    auto exp_arg_of = [](const Expr& f) -> std::optional<Expr> {
        // exp(a) → a; (exp(a))^k → k·a.
        if (f->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*f);
            if (fn.function_id() == FunctionId::Exp && fn.args().size() == 1) {
                return fn.args()[0];
            }
        }
        if (f->type_id() == TypeId::Pow) {
            const Expr& b = f->args()[0];
            if (b->type_id() == TypeId::Function) {
                const auto& bfn = static_cast<const Function&>(*b);
                if (bfn.function_id() == FunctionId::Exp
                    && bfn.args().size() == 1) {
                    return mul(f->args()[1], bfn.args()[0]);
                }
            }
        }
        return std::nullopt;
    };
    Expr exp_arg = S::Zero();
    int count = 0;
    std::vector<Expr> others;
    for (const auto& f : e->args()) {
        if (auto a = exp_arg_of(f); a.has_value()) {
            exp_arg = exp_arg + *a;
            ++count;
        } else {
            others.push_back(f);
        }
    }
    if (count < 2) return e;  // nothing to merge
    others.push_back(exp(exp_arg));
    return mul(std::move(others));
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

Expr pow_of_pow(const Expr& e) { return apply_recursive(e, pow_of_pow_node); }

Expr radical_coeff(const Expr& e) {
    return apply_recursive(e, radical_coeff_node);
}

Expr combine_exp(const Expr& e) { return apply_recursive(e, combine_exp_node); }

// a·exp(g) + a·exp(−g) → 2a·cosh(g),  a·exp(g) − a·exp(−g) → 2a·sinh(g).
// The mirror of hyp_to_exp_add (cosh±sinh → exp): collect, per argument g, the
// coefficients of exp(g) and exp(−g), and fold the equal / opposite-coefficient
// pairs into cosh / sinh. Unequal-coefficient sums (exp(x)+2·exp(−x)) are left.
[[nodiscard]] Expr exp_to_hyp_add_node(const Expr& e) {
    if (e->type_id() != TypeId::Add) return e;
    auto is_exp = [](const Expr& f) -> std::optional<Expr> {
        if (f->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*f);
            if (fn.function_id() == FunctionId::Exp && fn.args().size() == 1) {
                return fn.args()[0];
            }
        }
        return std::nullopt;
    };
    auto as_exp_term = [&](const Expr& term)
        -> std::optional<std::pair<Expr, Expr>> {  // (arg, coef)
        if (auto a = is_exp(term)) return std::pair{*a, Expr{S::One()}};
        if (term->type_id() == TypeId::Mul) {
            std::optional<Expr> arg;
            std::vector<Expr> coef;
            for (const auto& f : term->args()) {
                if (auto a = is_exp(f)) {
                    if (arg) return std::nullopt;  // exp·exp — leave for combine_exp
                    arg = *a;
                } else {
                    coef.push_back(f);
                }
            }
            if (arg) return std::pair{*arg, mul(std::move(coef))};
        }
        return std::nullopt;
    };

    struct PosNeg { Expr pos = S::Zero(); Expr neg = S::Zero(); };
    std::vector<std::pair<Expr, PosNeg>> buckets;  // keyed by the +g representative
    std::vector<Expr> rest;
    for (const auto& term : e->args()) {
        auto t = as_exp_term(term);
        if (!t) { rest.push_back(term); continue; }
        const Expr& arg = t->first;
        const Expr& coef = t->second;
        Expr narg = mul(S::NegativeOne(), arg);
        bool placed = false;
        for (auto& [g, pn] : buckets) {
            if (g == arg) { pn.pos = add(pn.pos, coef); placed = true; break; }
            if (g == narg) { pn.neg = add(pn.neg, coef); placed = true; break; }
        }
        if (!placed) buckets.push_back({arg, PosNeg{coef, S::Zero()}});
    }

    // Prefer the non-negated argument: cosh is even (cosh(−g)=cosh g), sinh odd
    // (sinh(−g)=−sinh g), so normalising g to its positive representative matches
    // SymPy's printed form (2·cosh(2x), not 2·cosh(−2x)).
    auto neg_leading = [](const Expr& g) -> bool {
        if (is_number(g)) return is_negative(g) == std::optional<bool>{true};
        if (g->type_id() == TypeId::Mul && !g->args().empty()
            && is_number(g->args()[0])) {
            return is_negative(g->args()[0]) == std::optional<bool>{true};
        }
        return false;
    };

    std::vector<Expr> out = rest;
    bool changed = false;
    for (auto& [g, pn] : buckets) {
        const bool has_pos = !(pn.pos == S::Zero());
        const bool has_neg = !(pn.neg == S::Zero());
        if (has_pos && pn.pos == pn.neg) {  // a·e^g + a·e^−g = 2a·cosh g
            Expr gg = neg_leading(g) ? mul(S::NegativeOne(), g) : g;
            out.push_back(mul(integer(2), mul(pn.pos, cosh(gg))));
            changed = true;
        } else if (has_pos && pn.pos == mul(S::NegativeOne(), pn.neg)) {
            // a·e^g − a·e^−g = 2a·sinh g
            Expr gg = g;
            Expr co = pn.pos;
            if (neg_leading(g)) {  // sinh(−g) = −sinh(g)
                gg = mul(S::NegativeOne(), g);
                co = mul(S::NegativeOne(), co);
            }
            out.push_back(mul(integer(2), mul(co, sinh(gg))));
            changed = true;
        } else {
            if (has_pos) out.push_back(mul(pn.pos, exp(g)));
            if (has_neg) out.push_back(mul(pn.neg, exp(mul(S::NegativeOne(), g))));
        }
    }
    if (!changed) return e;
    if (out.empty()) return S::Zero();
    return out.size() == 1 ? out[0] : add(std::move(out));
}

Expr exp_to_hyp_add(const Expr& e) {
    return apply_recursive(e, exp_to_hyp_add_node);
}

Expr exp_log_sum(const Expr& e) { return apply_recursive(e, exp_log_sum_node); }

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

[[nodiscard]] std::size_t count_leaves(const Expr& e) {
    if (!e) return 0;
    auto args = e->args();
    if (args.empty()) return 1;
    std::size_t total = 0;
    for (const auto& a : args) total += count_leaves(a);
    return total;
}

// If `f` is cos(2·g) — a Cos whose argument is a Mul carrying a literal
// integer-2 factor — return g; otherwise nullopt. Used to fold a double-angle
// cosine back into the sin²/cos² buckets via cos(2g) = cos²(g) − sin²(g), which
// lets trigsimp recover the power-reduction identities (1 − cos 2x)/2 = sin²x.
[[nodiscard]] std::optional<Expr> cos_double_arg(const Expr& f) {
    if (f->type_id() != TypeId::Function) return std::nullopt;
    const auto& fn = static_cast<const Function&>(*f);
    if (fn.function_id() != FunctionId::Cos) return std::nullopt;
    const Expr& a = f->args()[0];
    if (a->type_id() != TypeId::Mul) return std::nullopt;
    std::vector<Expr> rest;
    bool found_two = false;
    for (const auto& g : a->args()) {
        if (!found_two && g == integer(2)) { found_two = true; continue; }
        rest.push_back(g);
    }
    if (!found_two || rest.empty()) return std::nullopt;
    return mul(std::move(rest));
}

// Decompose a term as coef·cos(2·arg). Mirrors as_trig_square_term but for the
// double-angle cosine. Returns nullopt when there is no such cosine, or more
// than one (an unusual structure we leave alone).
struct CosDoubleTerm {
    Expr arg;
    Expr coef;
};
[[nodiscard]] std::optional<CosDoubleTerm> as_cos_double_term(const Expr& term) {
    if (auto g = cos_double_arg(term); g) {
        return CosDoubleTerm{*g, S::One()};
    }
    if (term->type_id() == TypeId::Mul) {
        std::optional<Expr> arg;
        std::vector<Expr> coef_factors;
        for (const auto& f : term->args()) {
            if (auto g = cos_double_arg(f); g) {
                if (arg) return std::nullopt;  // multiple — leave alone
                arg = *g;
            } else {
                coef_factors.push_back(f);
            }
        }
        if (arg) {
            return CosDoubleTerm{*arg, mul(std::move(coef_factors))};
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
        if (auto t = as_trig_square_term(term)) {
            auto& cp = find_or_create(t->arg);
            if (t->is_sin) cp.sin_coef = add(cp.sin_coef, t->coef);
            else cp.cos_coef = add(cp.cos_coef, t->coef);
            continue;
        }
        // q·cos(2x) = q·(cos²x − sin²x): fold into the same buckets so the
        // candidate comparison below can pick the power-reduced sin²/cos² form.
        if (auto d = as_cos_double_term(term)) {
            auto& cp = find_or_create(d->arg);
            cp.cos_coef = add(cp.cos_coef, d->coef);
            cp.sin_coef = add(cp.sin_coef, S::NegativeOne() * d->coef);
            continue;
        }
        non_trig.push_back(term);
    }

    if (by_arg.empty()) return e;

    // Sin Pythagorean candidate: a·sin²(x) + b·cos²(x) → b + (a − b)·sin²(x).
    auto sin_pythagorean_form = [&] {
        std::vector<Expr> out = non_trig;
        for (auto& [arg, cp] : by_arg) {
            out.push_back(cp.cos_coef);
            Expr diff = cp.sin_coef - cp.cos_coef;
            if (!(diff == S::Zero())) {
                out.push_back(diff * pow(sin(arg), integer(2)));
            }
        }
        if (out.empty()) return Expr{S::Zero()};
        if (out.size() == 1) return out[0];
        return add(std::move(out));
    }();

    // Cos Pythagorean candidate: a·sin²(x) + b·cos²(x) → a + (b − a)·cos²(x).
    // Needed so e.g. (1 + cos 2x)/2 reduces to cos²(x) rather than 1 − sin²(x).
    auto cos_pythagorean_form = [&] {
        std::vector<Expr> out = non_trig;
        for (auto& [arg, cp] : by_arg) {
            out.push_back(cp.sin_coef);
            Expr diff = cp.cos_coef - cp.sin_coef;
            if (!(diff == S::Zero())) {
                out.push_back(diff * pow(cos(arg), integer(2)));
            }
        }
        if (out.empty()) return Expr{S::Zero()};
        if (out.size() == 1) return out[0];
        return add(std::move(out));
    }();

    // Double-angle candidate: same Add rewritten with
    //   sin²(x) = (1 − cos(2x))/2,  cos²(x) = (1 + cos(2x))/2.
    // Per arg:  a·sin² + b·cos² = (a + b)/2 + ((b − a)/2)·cos(2x).
    // This automatically handles the constant-mixed shapes
    // (1 − 2·sin²(x), 2·cos²(x) − 1) by absorbing them into the Add.
    auto double_angle_form = [&] {
        std::vector<Expr> out = non_trig;
        Expr half = rational(1, 2);
        for (auto& [arg, cp] : by_arg) {
            out.push_back(half * (cp.sin_coef + cp.cos_coef));
            Expr c_2x = half * (cp.cos_coef - cp.sin_coef);
            if (!(c_2x == S::Zero())) {
                out.push_back(c_2x * cos(integer(2) * arg));
            }
        }
        if (out.empty()) return Expr{S::Zero()};
        if (out.size() == 1) return out[0];
        return add(std::move(out));
    }();

    // Pick whichever of the three equivalent forms has the fewest leaves.
    Expr best = sin_pythagorean_form;
    if (count_leaves(cos_pythagorean_form) < count_leaves(best)) {
        best = cos_pythagorean_form;
    }
    if (count_leaves(double_angle_form) < count_leaves(best)) {
        best = double_angle_form;
    }
    return best;
}

// 2·sin(x)·cos(x) = sin(2x).  More generally, k·sin(x)·cos(x) → (k/2)·sin(2x).
// Walks each Mul looking for at least one sin and one cos with the same
// argument; folds the pair into a single sin(2·arg) and leaves everything
// else untouched.
[[nodiscard]] Expr trigsimp_mul(const Expr& e) {
    if (e->type_id() != TypeId::Mul) return e;
    auto args = e->args();
    Expr sin_arg = nullptr;
    std::size_t sin_idx = 0;
    for (std::size_t i = 0; i < args.size(); ++i) {
        const Expr& f = args[i];
        if (f->type_id() != TypeId::Function) continue;
        const auto& fn = static_cast<const Function&>(*f);
        if (fn.function_id() == FunctionId::Sin) {
            sin_arg = f->args()[0];
            sin_idx = i;
            break;
        }
    }
    if (!sin_arg) return e;
    std::size_t cos_idx = args.size();
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (i == sin_idx) continue;
        const Expr& f = args[i];
        if (f->type_id() != TypeId::Function) continue;
        const auto& fn = static_cast<const Function&>(*f);
        if (fn.function_id() == FunctionId::Cos
            && f->args()[0] == sin_arg) {
            cos_idx = i;
            break;
        }
    }
    if (cos_idx == args.size()) return e;
    std::vector<Expr> rest;
    rest.reserve(args.size() - 1);
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (i == sin_idx || i == cos_idx) continue;
        rest.push_back(args[i]);
    }
    rest.push_back(rational(1, 2));
    rest.push_back(sin(integer(2) * sin_arg));
    return mul(std::move(rest));
}

// Classify an Add term as a product of exactly two first-power sin/cos factors,
// returning the coefficient and each factor's (is_sin, arg). Anything else — a
// bare or squared trig, a third trig factor, or a leftover function in the
// coefficient — returns nullopt so the term is left untouched.
struct TwoTrigTerm {
    Expr coef;
    bool sin1;
    Expr arg1;
    bool sin2;
    Expr arg2;
};
[[nodiscard]] std::optional<TwoTrigTerm> as_two_trig_term(const Expr& term) {
    if (term->type_id() != TypeId::Mul) return std::nullopt;
    std::vector<std::pair<bool, Expr>> trigs;  // (is_sin, arg)
    std::vector<Expr> coef_factors;
    for (const auto& f : term->args()) {
        if (f->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*f);
            if (fn.function_id() == FunctionId::Sin) {
                trigs.push_back({true, f->args()[0]});
                continue;
            }
            if (fn.function_id() == FunctionId::Cos) {
                trigs.push_back({false, f->args()[0]});
                continue;
            }
            return std::nullopt;  // another function (tan, …) in the term
        }
        // A power of a trig (sin²) or other trig-bearing factor is not a clean
        // coefficient — bail rather than risk a wrong fold.
        if (f->type_id() == TypeId::Pow
            && f->args()[0]->type_id() == TypeId::Function) {
            return std::nullopt;
        }
        coef_factors.push_back(f);
    }
    if (trigs.size() != 2) return std::nullopt;
    return TwoTrigTerm{coef_factors.empty() ? Expr{S::One()} : mul(coef_factors),
                       trigs[0].first, trigs[0].second, trigs[1].first,
                       trigs[1].second};
}

// Angle-addition identities on a two-term Add (each term a product of two
// trig factors), matching SymPy's simplify:
//   sin(a)cos(b) ± cos(a)sin(b) → sin(a ± b)
//   cos(a)cos(b) ∓ sin(a)sin(b) → cos(a ± b)
// A shared coefficient g carries through (g·… → g·sin/cos(…)).
[[nodiscard]] Expr trigsimp_angle_addition(const Expr& e) {
    if (e->type_id() != TypeId::Add || e->args().size() != 2) return e;
    auto t1 = as_two_trig_term(e->args()[0]);
    auto t2 = as_two_trig_term(e->args()[1]);
    if (!t1 || !t2) return e;
    auto neg = [](const Expr& c) { return mul(S::NegativeOne(), c); };

    const bool sc1 = t1->sin1 != t1->sin2;  // one sin, one cos
    const bool sc2 = t2->sin1 != t2->sin2;
    // sin(a±b): both terms are sin·cos.
    if (sc1 && sc2) {
        Expr a = t1->sin1 ? t1->arg1 : t1->arg2;   // sin argument of term 1
        Expr b = t1->sin1 ? t1->arg2 : t1->arg1;   // cos argument of term 1
        Expr s2 = t2->sin1 ? t2->arg1 : t2->arg2;  // sin argument of term 2
        Expr c2 = t2->sin1 ? t2->arg2 : t2->arg1;  // cos argument of term 2
        if (!(s2 == b && c2 == a)) return e;       // must cross-pair
        if (t1->coef == t2->coef) {
            return mul(t1->coef, sin(add({a, b})));
        }
        if (t1->coef == neg(t2->coef)) {
            return mul(t1->coef, sin(add({a, neg(b)})));
        }
        return e;
    }
    // cos(a±b): one term is cos·cos, the other sin·sin, over the same arg pair.
    const bool cc1 = !t1->sin1 && !t1->sin2, ss1 = t1->sin1 && t1->sin2;
    const bool cc2 = !t2->sin1 && !t2->sin2, ss2 = t2->sin1 && t2->sin2;
    const TwoTrigTerm* cc = nullptr;
    const TwoTrigTerm* ss = nullptr;
    if (cc1 && ss2) { cc = &*t1; ss = &*t2; }
    else if (ss1 && cc2) { cc = &*t2; ss = &*t1; }
    else return e;
    const bool same_args =
        (cc->arg1 == ss->arg1 && cc->arg2 == ss->arg2)
        || (cc->arg1 == ss->arg2 && cc->arg2 == ss->arg1);
    if (!same_args) return e;
    Expr a = cc->arg1, b = cc->arg2;
    if (cc->coef == ss->coef) {
        return mul(cc->coef, cos(add({a, neg(b)})));  // cos(a−b)
    }
    if (cc->coef == neg(ss->coef)) {
        return mul(cc->coef, cos(add({a, b})));  // cos(a+b)
    }
    return e;
}

// Hyperbolic Pythagorean: a·sinh²(x) + b·cosh²(x) using cosh² − sinh² = 1.
// Two equivalent rewrites — pick whichever (with the rest of the Add) is
// smallest, matching SymPy: cosh²−sinh² → 1, and 1+sinh² → cosh².
//   sinh form: b + (a+b)·sinh²(x)        (cosh² = 1 + sinh²)
//   cosh form: −a + (a+b)·cosh²(x)       (sinh² = cosh² − 1)
[[nodiscard]] Expr hypsimp_add(const Expr& e) {
    if (e->type_id() != TypeId::Add) return e;
    auto detect = [](const Expr& f) -> std::optional<std::pair<bool, Expr>> {
        if (f->type_id() != TypeId::Pow || !(f->args()[1] == integer(2))) {
            return std::nullopt;
        }
        const Expr& base = f->args()[0];
        if (base->type_id() != TypeId::Function) return std::nullopt;
        const auto& fn = static_cast<const Function&>(*base);
        if (fn.function_id() == FunctionId::Sinh) {
            return std::pair{true, base->args()[0]};
        }
        if (fn.function_id() == FunctionId::Cosh) {
            return std::pair{false, base->args()[0]};
        }
        return std::nullopt;
    };
    auto as_term = [&](const Expr& term)
        -> std::optional<std::tuple<bool, Expr, Expr>> {
        if (auto r = detect(term)) return std::make_tuple(r->first, r->second, S::One());
        if (term->type_id() == TypeId::Mul) {
            std::optional<std::pair<bool, Expr>> hyp;
            std::vector<Expr> coef;
            for (const auto& f : term->args()) {
                if (auto r = detect(f)) {
                    if (hyp) return std::nullopt;
                    hyp = r;
                } else {
                    coef.push_back(f);
                }
            }
            if (hyp) return std::make_tuple(hyp->first, hyp->second, mul(coef));
        }
        return std::nullopt;
    };

    struct CP { Expr sinh_c = S::Zero(); Expr cosh_c = S::Zero(); };
    std::vector<std::pair<Expr, CP>> by_arg;
    std::vector<Expr> non_hyp;
    auto find_or_create = [&](const Expr& a) -> CP& {
        for (auto& [x, p] : by_arg) if (x == a) return p;
        by_arg.push_back({a, CP{}});
        return by_arg.back().second;
    };
    for (const auto& term : e->args()) {
        if (auto t = as_term(term)) {
            auto& cp = find_or_create(std::get<1>(*t));
            if (std::get<0>(*t)) cp.sinh_c = add(cp.sinh_c, std::get<2>(*t));
            else cp.cosh_c = add(cp.cosh_c, std::get<2>(*t));
        } else {
            non_hyp.push_back(term);
        }
    }
    if (by_arg.empty()) return e;

    auto build = [&](bool sinh_form) {
        std::vector<Expr> out = non_hyp;
        for (auto& [arg, cp] : by_arg) {
            Expr sum = cp.sinh_c + cp.cosh_c;
            if (sinh_form) {
                out.push_back(cp.cosh_c);
                if (!(sum == S::Zero())) out.push_back(sum * pow(sinh(arg), integer(2)));
            } else {
                out.push_back(mul(S::NegativeOne(), cp.sinh_c));
                if (!(sum == S::Zero())) out.push_back(sum * pow(cosh(arg), integer(2)));
            }
        }
        if (out.empty()) return Expr{S::Zero()};
        if (out.size() == 1) return out[0];
        return add(std::move(out));
    };
    Expr sinh_form = build(true);
    Expr cosh_form = build(false);
    Expr best = count_leaves(cosh_form) < count_leaves(sinh_form)
                    ? cosh_form : sinh_form;
    return count_leaves(best) < count_leaves(e) ? best : e;
}

// Additive tanh/coth Pythagorean identities, using tanh² = 1 − sech² and
// coth² = 1 + csch² (with sech² = cosh⁻², csch² = sinh⁻²). Each squared
// tanh/coth/sech/csch term is rewritten into the cosh⁻²/sinh⁻² basis plus a
// constant; the rewrite is kept only when it has strictly fewer leaves — which
// happens exactly when the loose constant cancels:
//   1 − tanh²x → cosh⁻²x,  coth²x − 1 → sinh⁻²x,
//   sech²x + tanh²x → 1,   csch²x − coth²x → −1,  3 − 3tanh²x → 3cosh⁻²x.
// A bare tanh²x (or 2 − tanh²x, where the constant survives) is left untouched.
[[nodiscard]] Expr tanh_coth_pyth_add(const Expr& e) {
    if (e->type_id() != TypeId::Add) return e;
    // Detect a squared tanh/coth/sech/csch factor → (FunctionId, argument).
    auto detect = [](const Expr& f)
        -> std::optional<std::pair<FunctionId, Expr>> {
        if (f->type_id() != TypeId::Pow || !(f->args()[1] == integer(2))) {
            return std::nullopt;
        }
        const Expr& base = f->args()[0];
        if (base->type_id() != TypeId::Function) return std::nullopt;
        const auto& fn = static_cast<const Function&>(*base);
        switch (fn.function_id()) {
            case FunctionId::Tanh:
            case FunctionId::Coth:
            case FunctionId::Sech:
            case FunctionId::Csch:
                return std::pair{fn.function_id(), base->args()[0]};
            default:
                return std::nullopt;
        }
    };
    auto as_term = [&](const Expr& term)
        -> std::optional<std::tuple<FunctionId, Expr, Expr>> {
        if (auto r = detect(term)) {
            return std::make_tuple(r->first, r->second, S::One());
        }
        if (term->type_id() == TypeId::Mul) {
            std::optional<std::pair<FunctionId, Expr>> hyp;
            std::vector<Expr> coef;
            for (const auto& f : term->args()) {
                if (auto r = detect(f)) {
                    if (hyp) return std::nullopt;
                    hyp = r;
                } else {
                    coef.push_back(f);
                }
            }
            if (hyp) {
                return std::make_tuple(hyp->first, hyp->second, mul(coef));
            }
        }
        return std::nullopt;
    };

    struct Coefs {
        Expr tanh2 = S::Zero();
        Expr coth2 = S::Zero();
        Expr sech2 = S::Zero();
        Expr csch2 = S::Zero();
    };
    std::vector<std::pair<Expr, Coefs>> by_arg;
    std::vector<Expr> rest;
    auto find_or_create = [&](const Expr& a) -> Coefs& {
        for (auto& [x, c] : by_arg) if (x == a) return c;
        by_arg.push_back({a, Coefs{}});
        return by_arg.back().second;
    };
    bool saw_convertible = false;
    for (const auto& term : e->args()) {
        if (auto t = as_term(term)) {
            auto& c = find_or_create(std::get<1>(*t));
            const Expr& k = std::get<2>(*t);
            switch (std::get<0>(*t)) {
                case FunctionId::Tanh: c.tanh2 = add(c.tanh2, k); break;
                case FunctionId::Coth: c.coth2 = add(c.coth2, k); break;
                case FunctionId::Sech: c.sech2 = add(c.sech2, k); break;
                case FunctionId::Csch: c.csch2 = add(c.csch2, k); break;
                default: break;
            }
            if (std::get<0>(*t) == FunctionId::Tanh
                || std::get<0>(*t) == FunctionId::Coth) {
                saw_convertible = true;
            }
        } else {
            rest.push_back(term);
        }
    }
    // Only the tanh²/coth² → constant rewrites can make a constant cancel; a
    // pure sech²/csch² sum has nothing to gain, so skip it.
    if (!saw_convertible) return e;

    std::vector<Expr> out = rest;
    for (auto& [arg, c] : by_arg) {
        // tanh² = 1 − cosh⁻²,  coth² = 1 + sinh⁻²,  sech² = cosh⁻²,  csch² = sinh⁻².
        out.push_back(c.tanh2 + c.coth2);  // loose constant from the conversions
        Expr cosh_inv2 = c.sech2 - c.tanh2;
        Expr sinh_inv2 = c.csch2 + c.coth2;
        if (!(cosh_inv2 == S::Zero())) {
            out.push_back(cosh_inv2 * pow(cosh(arg), integer(-2)));
        }
        if (!(sinh_inv2 == S::Zero())) {
            out.push_back(sinh_inv2 * pow(sinh(arg), integer(-2)));
        }
    }
    Expr rebuilt = out.empty() ? Expr{S::Zero()}
                   : out.size() == 1 ? out[0]
                                     : add(std::move(out));
    // These identities are wins only when they shrink the sum; otherwise they
    // just trade a tanh²/coth² for an equally complex 1 + sech² form (which
    // SymPy leaves alone). Compare the number of additive terms.
    auto term_count = [](const Expr& x) -> std::size_t {
        return x->type_id() == TypeId::Add ? x->args().size() : 1;
    };
    return term_count(rebuilt) < e->args().size() ? rebuilt : e;
}

// Additive trig Pythagorean identities (the analogue of tanh_coth_pyth_add, with
// the opposite sign since sec² − tan² = 1, csc² − cot² = 1). Each squared
// tan/cot/sec/csc term is rewritten into the cos⁻²/sin⁻² basis via
// tan² = cos⁻² − 1, cot² = sin⁻² − 1, sec² = cos⁻², csc² = sin⁻², accumulating a
// constant; kept only when it shrinks the number of additive terms:
//   1 + tan²x → cos⁻²x,  sec²x − tan²x → 1,  csc²x − cot²x → 1,
//   1 + cot²x → sin⁻²x,  tan²x − sec²x → −1.
// A bare tan²x (or 2 + tan²x, where the constant survives) is left untouched.
[[nodiscard]] Expr trig_pyth_add(const Expr& e) {
    if (e->type_id() != TypeId::Add) return e;
    auto detect = [](const Expr& f)
        -> std::optional<std::pair<FunctionId, Expr>> {
        if (f->type_id() != TypeId::Pow || !(f->args()[1] == integer(2))) {
            return std::nullopt;
        }
        const Expr& base = f->args()[0];
        if (base->type_id() != TypeId::Function) return std::nullopt;
        const auto& fn = static_cast<const Function&>(*base);
        switch (fn.function_id()) {
            case FunctionId::Tan:
            case FunctionId::Cot:
            case FunctionId::Sec:
            case FunctionId::Csc:
                return std::pair{fn.function_id(), base->args()[0]};
            default:
                return std::nullopt;
        }
    };
    auto as_term = [&](const Expr& term)
        -> std::optional<std::tuple<FunctionId, Expr, Expr>> {
        if (auto r = detect(term)) {
            return std::make_tuple(r->first, r->second, S::One());
        }
        if (term->type_id() == TypeId::Mul) {
            std::optional<std::pair<FunctionId, Expr>> trig;
            std::vector<Expr> coef;
            for (const auto& f : term->args()) {
                if (auto r = detect(f)) {
                    if (trig) return std::nullopt;
                    trig = r;
                } else {
                    coef.push_back(f);
                }
            }
            if (trig) {
                return std::make_tuple(trig->first, trig->second, mul(coef));
            }
        }
        return std::nullopt;
    };

    struct Coefs {
        Expr tan2 = S::Zero();
        Expr cot2 = S::Zero();
        Expr sec2 = S::Zero();
        Expr csc2 = S::Zero();
    };
    std::vector<std::pair<Expr, Coefs>> by_arg;
    std::vector<Expr> rest;
    auto find_or_create = [&](const Expr& a) -> Coefs& {
        for (auto& [x, c] : by_arg) if (x == a) return c;
        by_arg.push_back({a, Coefs{}});
        return by_arg.back().second;
    };
    bool saw_convertible = false;
    for (const auto& term : e->args()) {
        if (auto t = as_term(term)) {
            auto& c = find_or_create(std::get<1>(*t));
            const Expr& k = std::get<2>(*t);
            switch (std::get<0>(*t)) {
                case FunctionId::Tan: c.tan2 = add(c.tan2, k); break;
                case FunctionId::Cot: c.cot2 = add(c.cot2, k); break;
                case FunctionId::Sec: c.sec2 = add(c.sec2, k); break;
                case FunctionId::Csc: c.csc2 = add(c.csc2, k); break;
                default: break;
            }
            if (std::get<0>(*t) == FunctionId::Tan
                || std::get<0>(*t) == FunctionId::Cot) {
                saw_convertible = true;
            }
        } else {
            rest.push_back(term);
        }
    }
    if (!saw_convertible) return e;

    std::vector<Expr> out = rest;
    for (auto& [arg, c] : by_arg) {
        // tan² = cos⁻² − 1,  cot² = sin⁻² − 1,  sec² = cos⁻²,  csc² = sin⁻².
        out.push_back(mul(S::NegativeOne(), c.tan2 + c.cot2));  // loose constant
        Expr cos_inv2 = c.tan2 + c.sec2;
        Expr sin_inv2 = c.cot2 + c.csc2;
        if (!(cos_inv2 == S::Zero())) {
            out.push_back(cos_inv2 * pow(cos(arg), integer(-2)));
        }
        if (!(sin_inv2 == S::Zero())) {
            out.push_back(sin_inv2 * pow(sin(arg), integer(-2)));
        }
    }
    Expr rebuilt = out.empty() ? Expr{S::Zero()}
                   : out.size() == 1 ? out[0]
                                     : add(std::move(out));
    auto term_count = [](const Expr& x) -> std::size_t {
        return x->type_id() == TypeId::Add ? x->args().size() : 1;
    };
    return term_count(rebuilt) < e->args().size() ? rebuilt : e;
}

// cosh(x) + sinh(x) → eˣ and cosh(x) − sinh(x) → e⁻ˣ (and scaled forms, when the
// cosh and sinh coefficients match up to a sign). Collapses a 2-term sum into a
// single exponential, matching SymPy's simplify.
[[nodiscard]] Expr hyp_to_exp_add(const Expr& e) {
    if (e->type_id() != TypeId::Add) return e;
    auto detect = [](const Expr& f) -> std::optional<std::pair<bool, Expr>> {
        // bare cosh(a) / sinh(a) (degree 1, not squared)
        if (f->type_id() != TypeId::Function) return std::nullopt;
        const auto& fn = static_cast<const Function&>(*f);
        if (fn.args().size() != 1) return std::nullopt;
        if (fn.function_id() == FunctionId::Cosh) return std::pair{false, f->args()[0]};
        if (fn.function_id() == FunctionId::Sinh) return std::pair{true, f->args()[0]};
        return std::nullopt;
    };
    auto as_term = [&](const Expr& term)
        -> std::optional<std::tuple<bool, Expr, Expr>> {
        if (auto r = detect(term)) return std::make_tuple(r->first, r->second, S::One());
        if (term->type_id() == TypeId::Mul) {
            std::optional<std::pair<bool, Expr>> hyp;
            std::vector<Expr> coef;
            for (const auto& f : term->args()) {
                if (!hyp) { if (auto r = detect(f)) { hyp = r; continue; } }
                coef.push_back(f);
            }
            if (hyp && !coef.empty()) {
                return std::make_tuple(hyp->first, hyp->second, mul(coef));
            }
        }
        return std::nullopt;
    };

    struct CP { Expr sinh_c = S::Zero(); Expr cosh_c = S::Zero(); };
    std::vector<std::pair<Expr, CP>> by_arg;
    std::vector<Expr> rest;
    auto find_or_create = [&](const Expr& a) -> CP& {
        for (auto& [x, p] : by_arg) if (x == a) return p;
        by_arg.push_back({a, CP{}});
        return by_arg.back().second;
    };
    for (const auto& term : e->args()) {
        if (auto t = as_term(term)) {
            auto& cp = find_or_create(std::get<1>(*t));
            if (std::get<0>(*t)) cp.sinh_c = add(cp.sinh_c, std::get<2>(*t));
            else cp.cosh_c = add(cp.cosh_c, std::get<2>(*t));
        } else {
            rest.push_back(term);
        }
    }

    std::vector<Expr> out = rest;
    bool changed = false;
    for (auto& [arg, cp] : by_arg) {
        if (!(cp.cosh_c == S::Zero()) && cp.sinh_c == cp.cosh_c) {
            out.push_back(cp.cosh_c * exp(arg));            // c(cosh+sinh) = c·eˣ
            changed = true;
        } else if (!(cp.cosh_c == S::Zero())
                   && cp.sinh_c == mul(S::NegativeOne(), cp.cosh_c)) {
            out.push_back(cp.cosh_c * exp(mul(S::NegativeOne(), arg)));
            changed = true;
        } else {
            if (!(cp.cosh_c == S::Zero())) out.push_back(cp.cosh_c * cosh(arg));
            if (!(cp.sinh_c == S::Zero())) out.push_back(cp.sinh_c * sinh(arg));
        }
    }
    if (!changed) return e;
    if (out.empty()) return S::Zero();
    if (out.size() == 1) return out[0];
    return add(std::move(out));
}

// Cancel hyperbolic ratio products inside a Mul by rewriting tanh/coth/sech/csch
// as sinh/cosh ratios and letting Mul recombine same-base powers:
//   tanh·cosh → sinh,  coth·sinh → cosh,  sech·cosh → 1,  csch·sinh → 1,
//   coth·tanh → 1, …   Each factor f^n is replaced by the sinh/cosh power(s)
// it stands for; the rewrite is kept only when it reduces the leaf count, so a
// lone tanh (or 2·tanh) — which would only grow — is left untouched.
[[nodiscard]] Expr hyp_ratio_mul(const Expr& e) {
    if (e->type_id() != TypeId::Mul) return e;
    bool changed = false;
    std::vector<Expr> factors;
    for (const auto& f : e->args()) {
        Expr base = f;
        Expr expn = S::One();
        if (f->type_id() == TypeId::Pow) {
            base = f->args()[0];
            expn = f->args()[1];
        }
        if (base->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*base);
            if (fn.args().size() == 1) {
                const Expr& a = base->args()[0];
                Expr neg = mul(S::NegativeOne(), expn);
                bool handled = true;
                switch (fn.function_id()) {
                    case FunctionId::Tanh:
                        factors.push_back(pow(sinh(a), expn));
                        factors.push_back(pow(cosh(a), neg));
                        break;
                    case FunctionId::Coth:
                        factors.push_back(pow(cosh(a), expn));
                        factors.push_back(pow(sinh(a), neg));
                        break;
                    case FunctionId::Sech:
                        factors.push_back(pow(cosh(a), neg));
                        break;
                    case FunctionId::Csch:
                        factors.push_back(pow(sinh(a), neg));
                        break;
                    default:
                        handled = false;
                        break;
                }
                if (handled) {
                    changed = true;
                    continue;
                }
            }
        }
        factors.push_back(f);
    }
    if (!changed) return e;
    Expr rebuilt = mul(std::move(factors));
    return count_leaves(rebuilt) < count_leaves(e) ? rebuilt : e;
}

// Trig analogue of hyp_ratio_mul: cancel ratio products inside a Mul by
// rewriting tan/cot/sec/csc as sin/cos ratios and letting Mul recombine
// same-base powers:
//   tan·cos → sin,  cot·sin → cos,  sec·cos → 1,  csc·sin → 1,  cot·tan → 1, …
// Kept only when it reduces the leaf count, so a lone tan (or 2·tan) is left
// untouched.
[[nodiscard]] Expr trig_ratio_mul(const Expr& e) {
    if (e->type_id() != TypeId::Mul) return e;
    bool changed = false;
    std::vector<Expr> factors;
    for (const auto& f : e->args()) {
        Expr base = f;
        Expr expn = S::One();
        if (f->type_id() == TypeId::Pow) {
            base = f->args()[0];
            expn = f->args()[1];
        }
        if (base->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*base);
            if (fn.args().size() == 1) {
                const Expr& a = base->args()[0];
                Expr neg = mul(S::NegativeOne(), expn);
                bool handled = true;
                switch (fn.function_id()) {
                    case FunctionId::Tan:
                        factors.push_back(pow(sin(a), expn));
                        factors.push_back(pow(cos(a), neg));
                        break;
                    case FunctionId::Cot:
                        factors.push_back(pow(cos(a), expn));
                        factors.push_back(pow(sin(a), neg));
                        break;
                    case FunctionId::Sec:
                        factors.push_back(pow(cos(a), neg));
                        break;
                    case FunctionId::Csc:
                        factors.push_back(pow(sin(a), neg));
                        break;
                    default:
                        handled = false;
                        break;
                }
                if (handled) {
                    changed = true;
                    continue;
                }
            }
        }
        factors.push_back(f);
    }
    if (!changed) return e;
    Expr rebuilt = mul(std::move(factors));
    return count_leaves(rebuilt) < count_leaves(e) ? rebuilt : e;
}

[[nodiscard]] Expr trigsimp_node(const Expr& e) {
    Expr cur = trigsimp_add(e);
    cur = trigsimp_angle_addition(cur);
    cur = hypsimp_add(cur);
    cur = tanh_coth_pyth_add(cur);
    cur = trig_pyth_add(cur);
    cur = hyp_to_exp_add(cur);
    cur = trig_ratio_mul(cur);
    cur = trigsimp_mul(cur);
    cur = hyp_ratio_mul(cur);
    return cur;
}

}  // namespace

Expr trigsimp(const Expr& e) {
    return apply_recursive(e, trigsimp_node);
}

// ----- expand_trig + fu ------------------------------------------------------

namespace {

[[nodiscard]] Expr expand_trig_node(const Expr& e) {
    if (e->type_id() != TypeId::Function) return e;
    const auto& f = static_cast<const Function&>(*e);
    auto fid = f.function_id();
    if (fid != FunctionId::Sin && fid != FunctionId::Cos
        && fid != FunctionId::Tan) return e;
    const Expr& arg = e->args()[0];
    // Split the argument into a + b, then apply the angle-addition formula.
    // Two sources of a split:
    //   • an Add argument (sin(x+y) → …): a = first term, b = the rest;
    //   • a multiple angle n·g with integer n ≥ 2 (sin(2x), cos(3x), …):
    //     n·g = g + (n−1)·g, and the recursive pass reduces (n−1)·g.
    Expr a, b;
    if (arg->type_id() == TypeId::Add) {
        auto add_args = arg->args();
        if (add_args.size() < 2) return e;
        a = add_args[0];
        if (add_args.size() == 2) {
            b = add_args[1];
        } else {
            std::vector<Expr> rest(add_args.begin() + 1, add_args.end());
            b = add(std::move(rest));
        }
    } else if (arg->type_id() == TypeId::Mul) {
        long n = 1;
        bool found = false;
        std::vector<Expr> rest;
        for (const auto& fac : arg->args()) {
            if (!found && fac->type_id() == TypeId::Integer) {
                const auto& z = static_cast<const Integer&>(*fac);
                if (z.is_positive() && z.fits_long() && z.to_long() >= 2) {
                    n = z.to_long();
                    found = true;
                    continue;
                }
            }
            rest.push_back(fac);
        }
        if (!found || rest.empty()) return e;
        Expr g = rest.size() == 1 ? rest[0] : mul(std::move(rest));
        a = g;
        b = (n == 2) ? g : mul(integer(n - 1), g);
    } else {
        return e;
    }
    if (fid == FunctionId::Sin) {
        return sin(a) * cos(b) + cos(a) * sin(b);
    }
    if (fid == FunctionId::Cos) {
        return cos(a) * cos(b) - sin(a) * sin(b);
    }
    // Tan(a+b) = (tan a + tan b) / (1 - tan a · tan b).
    Expr ta = tan(a);
    Expr tb = tan(b);
    return (ta + tb) / (integer(1) - ta * tb);
}

}  // namespace

Expr expand_trig(const Expr& e) {
    // Iterate to a fixpoint: each pass peels one level of angle addition, so a
    // multiple angle sin(n·x) needs ~n passes to fully reduce to sin(x)/cos(x).
    Expr cur = apply_recursive(e, expand_trig_node);
    for (int i = 0; i < 32; ++i) {
        Expr next = apply_recursive(cur, expand_trig_node);
        if (next == cur) break;
        cur = next;
    }
    return cur;
}

Expr fu(const Expr& e) {
    Expr c1 = e;
    Expr c2 = trigsimp(e);
    Expr c3 = trigsimp(expand_trig(e));
    Expr best = c1;
    std::size_t best_n = count_leaves(c1);
    if (auto n = count_leaves(c2); n < best_n) { best = c2; best_n = n; }
    if (auto n = count_leaves(c3); n < best_n) { best = c3; best_n = n; }
    return best;
}

// ----- rewrite ---------------------------------------------------------------

namespace {

// sin/cos/tan and sinh/cosh/tanh → exponentials (Euler / hyperbolic identities).
[[nodiscard]] Expr rewrite_exp_node(const Expr& e) {
    if (e->type_id() != TypeId::Function) return e;
    const auto& fn = static_cast<const Function&>(*e);
    if (fn.args().size() != 1) return e;
    const Expr& g = fn.args()[0];
    const Expr neg = S::NegativeOne();
    const Expr& I = S::I();
    const Expr half = rational(1, 2);
    Expr Ig = mul(I, g);
    Expr nIg = mul(mul(neg, I), g);
    Expr ng = mul(neg, g);
    switch (fn.function_id()) {
        case FunctionId::Sin:  // −i·(e^{ig} − e^{−ig})/2
            return mul({mul(neg, I), half, exp(Ig) - exp(nIg)});
        case FunctionId::Cos:  // (e^{ig} + e^{−ig})/2
            return mul(half, exp(Ig) + exp(nIg));
        case FunctionId::Tan:  // −i·(e^{ig} − e^{−ig})/(e^{ig} + e^{−ig})
            return mul(mul(neg, I), (exp(Ig) - exp(nIg)) / (exp(Ig) + exp(nIg)));
        case FunctionId::Sinh:  // (e^g − e^{−g})/2
            return mul(half, exp(g) - exp(ng));
        case FunctionId::Cosh:  // (e^g + e^{−g})/2
            return mul(half, exp(g) + exp(ng));
        case FunctionId::Tanh:  // (e^g − e^{−g})/(e^g + e^{−g})
            return (exp(g) - exp(ng)) / (exp(g) + exp(ng));
        default:
            return e;
    }
}

}  // namespace

Expr rewrite(const Expr& e, std::string_view target) {
    if (target == "exp") return apply_recursive(e, rewrite_exp_node);
    return e;  // unknown target — no-op, matching SymPy's "leave as is"
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

    Expr den;
    Expr num;
    if (t->type_id() == TypeId::Pow && t->args()[1] == S::NegativeOne()) {
        // A bare reciprocal 1/den (e.g. 1/(1+√2)) — no numerator factors.
        den = t->args()[0];
        num = S::One();
    } else if (t->type_id() == TypeId::Mul) {
        // Locate the Pow(_, -1) factor — the denominator.
        std::vector<Expr> num_factors;
        for (const auto& f : t->args()) {
            if (f->type_id() == TypeId::Pow && f->args()[1] == S::NegativeOne()
                && !den) {
                den = f->args()[0];
            } else {
                num_factors.push_back(f);
            }
        }
        if (!den) return t;
        num = mul(num_factors);
    } else {
        return t;
    }

    auto decomp = decompose_sqrts(den);
    Expr conjugate;
    Expr new_den;
    if (decomp.sqrt_terms.size() == 1) {
        // Binomial: den = a + b·√c. Conjugate a − b·√c; (a+b√c)(a−b√c) = a²−b²c.
        const auto& [b_coef, sqrt_arg] = decomp.sqrt_terms[0];
        const Expr& a = decomp.rational_part;
        conjugate = a - mul(b_coef, sqrt(sqrt_arg));
        new_den = a * a - b_coef * b_coef * sqrt_arg;
    } else if (decomp.sqrt_terms.size() == 2
               && decomp.rational_part == S::Zero()) {
        // Two-surd: den = b₁·√c₁ + b₂·√c₂. Conjugate b₁·√c₁ − b₂·√c₂;
        // the product is the rational b₁²c₁ − b₂²c₂.
        const auto& [b1, c1] = decomp.sqrt_terms[0];
        const auto& [b2, c2] = decomp.sqrt_terms[1];
        conjugate = mul(b1, sqrt(c1)) - mul(b2, sqrt(c2));
        new_den = b1 * b1 * c1 - b2 * b2 * c2;
    } else {
        return t;
    }
    if (new_den == S::Zero()) return t;
    Expr new_num = num * conjugate;
    // Distribute so the result is the compact rationalized form (√2−1, not
    // −(−√2+1)); this also keeps it small enough to survive simplify()'s
    // anti-bloat guard, which otherwise reverts to the reciprocal.
    return expand(new_num / new_den);
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

// ----- combsimp / gammasimp --------------------------------------------------

namespace {

[[nodiscard]] std::optional<Expr> as_factorial_arg(const Expr& e) {
    if (!e || e->type_id() != TypeId::Function) return std::nullopt;
    const auto& f = static_cast<const Function&>(*e);
    if (f.function_id() != FunctionId::Factorial) return std::nullopt;
    return e->args()[0];
}

[[nodiscard]] std::optional<Expr> as_gamma_arg(const Expr& e) {
    if (!e || e->type_id() != TypeId::Function) return std::nullopt;
    const auto& f = static_cast<const Function&>(*e);
    if (f.function_id() != FunctionId::Gamma) return std::nullopt;
    return e->args()[0];
}

// Returns the falling product (a)(a-1)(a-2)...(a-k+1) — k factors total.
[[nodiscard]] Expr falling_factorial(const Expr& a, long k) {
    if (k <= 0) return S::One();
    Expr current = a;
    Expr product = current;
    for (long i = 1; i < k; ++i) {
        current = current - integer(static_cast<long>(1));
        product = product * current;
    }
    return product;
}

// Generic ratio simplifier: given a Mul with positive and negative factors
// of `func(arg)`, pair them up where (pos.arg - neg.arg) is a small non-
// negative integer k, replacing the pair with the falling product
// (pos.arg)*(pos.arg-1)*...*(neg.arg+1).
[[nodiscard]] Expr simplify_func_ratio(
    const Expr& e,
    std::optional<Expr> (*matcher)(const Expr&)) {
    if (e->type_id() != TypeId::Mul) return e;

    std::vector<Expr> pos_args;          // arg of positive func() factors
    std::vector<Expr> neg_args;          // arg of (func())^-1 factors
    std::vector<Expr> rest;
    for (const auto& f : e->args()) {
        if (auto a = matcher(f); a) {
            pos_args.push_back(*a);
            continue;
        }
        if (f->type_id() == TypeId::Pow
            && f->args()[1] == S::NegativeOne()) {
            if (auto a = matcher(f->args()[0]); a) {
                neg_args.push_back(*a);
                continue;
            }
        }
        rest.push_back(f);
    }

    // Pair up by integer-difference matching.
    std::vector<bool> pos_used(pos_args.size(), false);
    std::vector<bool> neg_used(neg_args.size(), false);
    std::vector<Expr> pairings;
    for (std::size_t i = 0; i < pos_args.size(); ++i) {
        for (std::size_t j = 0; j < neg_args.size(); ++j) {
            if (neg_used[j]) continue;
            // diff = pos_args[i] - neg_args[j]; if it's a small positive
            // integer k, factorial cancels into a falling product of k terms.
            Expr diff = simplify(pos_args[i] - neg_args[j]);
            if (diff->type_id() == TypeId::Integer) {
                const auto& z = static_cast<const Integer&>(*diff);
                if (z.fits_long()) {
                    long k = z.to_long();
                    if (k >= 0 && k <= 50) {
                        // For factorial(a)/factorial(b), product is
                        //   a*(a-1)*...*(b+1)  → falling_factorial(a, a-b)
                        // For gamma(a)/gamma(b), product is
                        //   (a-1)*(a-2)*...*b  → falling_factorial(a-1, a-b)
                        // Both are k = a - b terms.
                        Expr top = (matcher == as_factorial_arg)
                                       ? pos_args[i]
                                       : (pos_args[i] - integer(1));
                        pairings.push_back(falling_factorial(top, k));
                        pos_used[i] = true;
                        neg_used[j] = true;
                        break;
                    }
                    if (k < 0 && k >= -50) {
                        // Denominator is the larger factorial: the ratio is the
                        // reciprocal of the rising product. m = b - a terms.
                        //   factorial(a)/factorial(b) → 1/falling_factorial(b, m)
                        //   gamma(a)/gamma(b)         → 1/falling_factorial(b-1, m)
                        long m = -k;
                        Expr bot = (matcher == as_factorial_arg)
                                       ? neg_args[j]
                                       : (neg_args[j] - integer(1));
                        pairings.push_back(
                            pow(falling_factorial(bot, m), integer(-1)));
                        pos_used[i] = true;
                        neg_used[j] = true;
                        break;
                    }
                }
            }
        }
    }

    bool any_paired = std::any_of(pos_used.begin(), pos_used.end(),
                                  [](bool b) { return b; });
    if (!any_paired) return e;

    // Reassemble.
    std::vector<Expr> out = std::move(rest);
    for (auto& p : pairings) out.push_back(std::move(p));
    auto rebuild_func = [matcher](const Expr& arg) -> Expr {
        // matcher == as_factorial_arg → factorial; as_gamma_arg → gamma.
        if (matcher == as_factorial_arg) return factorial(arg);
        return gamma(arg);
    };
    for (std::size_t i = 0; i < pos_args.size(); ++i) {
        if (!pos_used[i]) out.push_back(rebuild_func(pos_args[i]));
    }
    for (std::size_t j = 0; j < neg_args.size(); ++j) {
        if (!neg_used[j]) out.push_back(pow(rebuild_func(neg_args[j]),
                                            integer(-1)));
    }
    return mul(std::move(out));
}

// Euler reflection: gamma(z)*gamma(1-z) -> pi / sin(pi*z). Scans the Mul for
// two numerator gamma factors whose arguments sum to 1 and folds each such
// pair. The surviving argument is chosen as the one free of a leading additive
// constant (prefer `z` over `1-z`) so the output matches SymPy's `pi/sin(pi*z)`.
[[nodiscard]] Expr gamma_reflection(const Expr& e) {
    if (e->type_id() != TypeId::Mul) return e;
    const auto& args = e->args();

    std::vector<Expr> g_args;            // arg of each numerator gamma factor
    std::vector<std::size_t> g_idx;      // its position in `args`
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (auto a = as_gamma_arg(args[i]); a) {
            g_args.push_back(*a);
            g_idx.push_back(i);
        }
    }

    std::vector<bool> consumed(args.size(), false);
    std::vector<bool> used(g_args.size(), false);
    std::vector<Expr> replacements;
    for (std::size_t i = 0; i < g_args.size(); ++i) {
        if (used[i]) continue;
        for (std::size_t j = i + 1; j < g_args.size(); ++j) {
            if (used[j]) continue;
            if (simplify(g_args[i] + g_args[j]) != S::One()) continue;
            // Prefer the argument without an additive integer term so the
            // result reads pi/sin(pi*z) rather than pi/sin(pi*(1-z)).
            const Expr& z =
                (g_args[j]->type_id() == TypeId::Add
                 && g_args[i]->type_id() != TypeId::Add)
                    ? g_args[i]
                    : (g_args[i]->type_id() == TypeId::Add ? g_args[j]
                                                           : g_args[i]);
            replacements.push_back(
                mul({S::Pi(), pow(sin(mul({S::Pi(), z})), integer(-1))}));
            used[i] = used[j] = true;
            consumed[g_idx[i]] = consumed[g_idx[j]] = true;
            break;
        }
    }
    if (replacements.empty()) return e;

    std::vector<Expr> out;
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (!consumed[i]) out.push_back(args[i]);
    }
    for (auto& r : replacements) out.push_back(std::move(r));
    return mul(std::move(out));
}

[[nodiscard]] std::optional<std::pair<Expr, Expr>> as_binomial(const Expr& e) {
    if (!e || e->type_id() != TypeId::Function) return std::nullopt;
    const auto& f = static_cast<const Function&>(*e);
    if (f.function_id() != FunctionId::Binomial) return std::nullopt;
    return std::make_pair(e->args()[0], e->args()[1]);
}

// combsimp on a binomial: collapse binomial(n, k) to its polynomial form when
// either k or n-k is a small non-negative integer m, via the gamma identity
//   binomial(n, k) = falling_factorial(n, m) / m!   (m = n-k, or m = k)
// which is valid for symbolic n. Examples: binomial(n,n)→1,
// binomial(n,n-1)→n, binomial(n+1,n)→n+1, binomial(n,2)→n*(n-1)/2.
[[nodiscard]] Expr combsimp_binomial(const Expr& e) {
    auto b = as_binomial(e);
    if (!b) return e;
    auto small_nonneg = [](const Expr& x) -> std::optional<long> {
        if (x->type_id() != TypeId::Integer) return std::nullopt;
        const auto& z = static_cast<const Integer&>(*x);
        if (!z.fits_long()) return std::nullopt;
        long v = z.to_long();
        if (v < 0 || v > 50) return std::nullopt;
        return v;
    };
    // m! as a single folded Integer (Mul folds numeric factors).
    auto m_factorial = [](long m) -> Expr {
        Expr r = S::One();
        for (long i = 2; i <= m; ++i) r = r * integer(i);
        return r;
    };
    const Expr& n = b->first;
    const Expr& k = b->second;
    // Prefer m = n-k (handles the symmetric tail: binomial(n,n), n-1, ...).
    std::optional<long> m = small_nonneg(simplify(n - k));
    if (!m) m = small_nonneg(k);   // otherwise m = k (small head: 0,1,2,...).
    if (!m) return e;
    return mul({falling_factorial(n, *m), pow(m_factorial(*m), integer(-1))});
}

// Gamma recurrence z*gamma(z) = gamma(z+1). Within a Mul, repeatedly absorb a
// numerator factor equal to a gamma argument z (raising it to gamma(z+1)) or a
// denominator factor equal to z-1 (lowering gamma(z) to gamma(z-1)). Iterates
// to a fixpoint so chains collapse: x*(x+1)*gamma(x) -> gamma(x+2).
[[nodiscard]] Expr gamma_recurrence(const Expr& e) {
    if (e->type_id() != TypeId::Mul) return e;
    std::vector<Expr> factors(e->args().begin(), e->args().end());
    bool changed = true;
    while (changed) {
        changed = false;
        for (std::size_t g = 0; g < factors.size() && !changed; ++g) {
            auto za = as_gamma_arg(factors[g]);
            if (!za) continue;
            const Expr z = *za;
            // Upward: a numerator factor equal to z.
            for (std::size_t i = 0; i < factors.size(); ++i) {
                if (i == g) continue;
                if (simplify(factors[i] - z) == S::Zero()) {
                    factors[g] = gamma(z + S::One());
                    factors.erase(factors.begin() + static_cast<long>(i));
                    changed = true;
                    break;
                }
            }
            if (changed) break;
            // Downward: a denominator factor f^-1 with f == z-1.
            for (std::size_t i = 0; i < factors.size(); ++i) {
                if (i == g) continue;
                if (factors[i]->type_id() == TypeId::Pow
                    && factors[i]->args()[1] == S::NegativeOne()
                    && simplify(factors[i]->args()[0] - (z - S::One()))
                           == S::Zero()) {
                    factors[g] = gamma(z - S::One());
                    factors.erase(factors.begin() + static_cast<long>(i));
                    changed = true;
                    break;
                }
            }
        }
    }
    return mul(std::move(factors));
}

// Does any subexpression apply the function `id`?
[[nodiscard]] bool contains_function_id(const Expr& e, FunctionId id) {
    if (!e) return false;
    if (e->type_id() == TypeId::Function
        && static_cast<const Function&>(*e).function_id() == id) {
        return true;
    }
    for (const auto& a : e->args()) {
        if (contains_function_id(a, id)) return true;
    }
    return false;
}

// combsimp on a ratio/product of binomials: rewrite binomial(a,b) =
// a!/(b!·(a−b)!) so the factorial-ratio cancellation closes shifts like
// binomial(n,k)/binomial(n,k−1) → (n−k+1)/k. Only adopted when the rewrite
// fully resolves (no factorial or binomial remains); otherwise binomial(n,k)
// would expand into an uglier factorial form, so the input is kept. Mirrors
// SymPy's combsimp, which collapses such ratios.
[[nodiscard]] Expr combsimp_binomial_ratio(const Expr& e) {
    if (e->type_id() != TypeId::Mul) return e;
    bool has_binom = false;
    std::vector<Expr> expanded;
    for (const auto& f : e->args()) {
        Expr base = f;
        Expr exp_e = S::One();
        if (f->type_id() == TypeId::Pow) {
            base = f->args()[0];
            exp_e = f->args()[1];
        }
        if (auto b = as_binomial(base)) {
            has_binom = true;
            const Expr& a = b->first;
            const Expr& k = b->second;
            Expr neg_exp = mul(S::NegativeOne(), exp_e);
            // binomial(a,k)^p = a!^p · k!^(−p) · (a−k)!^(−p).
            expanded.push_back(pow(factorial(a), exp_e));
            expanded.push_back(pow(factorial(k), neg_exp));
            expanded.push_back(pow(factorial(simplify(a - k)), neg_exp));
        } else {
            expanded.push_back(f);
        }
    }
    if (!has_binom) return e;
    Expr canceled = gamma_recurrence(
        simplify_func_ratio(mul(std::move(expanded)), as_factorial_arg));
    if (contains_function_id(canceled, FunctionId::Factorial)
        || contains_function_id(canceled, FunctionId::Binomial)) {
        return e;  // didn't fully resolve — keep the cleaner binomial form
    }
    return canceled;
}

[[nodiscard]] Expr combsimp_node(const Expr& e) {
    if (Expr ratio = combsimp_binomial_ratio(e); !(ratio == e)) return ratio;
    return gamma_recurrence(
        combsimp_binomial(simplify_func_ratio(e, as_factorial_arg)));
}

[[nodiscard]] Expr gammasimp_node(const Expr& e) {
    return gamma_recurrence(
        gamma_reflection(simplify_func_ratio(e, as_gamma_arg)));
}

}  // namespace

Expr combsimp(const Expr& e) { return apply_recursive(e, combsimp_node); }
Expr gammasimp(const Expr& e) { return apply_recursive(e, gammasimp_node); }

// ----- cse -------------------------------------------------------------------

namespace {

[[nodiscard]] int subtree_size(const Expr& e) {
    if (!e) return 0;
    int s = 1;
    for (const auto& a : e->args()) s += subtree_size(a);
    return s;
}

void count_subtrees(const Expr& e, ExprMap<int>& counts) {
    if (!e || e->args().empty()) return;
    counts[e]++;
    for (const auto& a : e->args()) count_subtrees(a, counts);
}

}  // namespace

CSEResult cse(const Expr& e) {
    ExprMap<int> counts;
    count_subtrees(e, counts);

    // Candidates appear ≥ 2 times.
    std::vector<Expr> candidates;
    for (const auto& [expr, cnt] : counts) {
        if (cnt >= 2) candidates.push_back(expr);
    }

    // Order by descending size — outermost subtrees first. After
    // substituting a large pattern, any smaller pattern entirely
    // contained in it that doesn't appear elsewhere will lose its
    // multi-count and stop being a candidate; we'll skip it implicitly
    // by re-checking counts on the rewritten expression at each step.
    std::sort(candidates.begin(), candidates.end(),
              [](const Expr& a, const Expr& b) {
                  int sa = subtree_size(a);
                  int sb = subtree_size(b);
                  if (sa != sb) return sa > sb;
                  return a->str() < b->str();  // tie-break stable
              });

    Expr current = e;
    std::vector<std::pair<Expr, Expr>> subs_list;
    int next_id = 0;
    for (auto& cand : candidates) {
        // Re-count on the current rewritten expression — earlier
        // substitutions may have removed occurrences.
        ExprMap<int> live_counts;
        count_subtrees(current, live_counts);
        auto it = live_counts.find(cand);
        if (it == live_counts.end() || it->second < 2) continue;

        Expr temp = symbol("_cse_" + std::to_string(next_id++));
        // Replace `cand` with `temp` everywhere in `current` and in
        // earlier substitution definitions.
        current = subs(current, cand, temp);
        for (auto& [t, def] : subs_list) {
            def = subs(def, cand, temp);
        }
        subs_list.emplace_back(temp, cand);
    }
    return CSEResult{std::move(subs_list), std::move(current)};
}

// ----- nsimplify -------------------------------------------------------------

namespace {

[[nodiscard]] std::optional<double> as_double(const Expr& e) {
    if (e->type_id() == TypeId::Integer) {
        return static_cast<const Integer&>(*e).value().get_d();
    }
    if (e->type_id() == TypeId::Rational) {
        return static_cast<const Rational&>(*e).value().get_d();
    }
    if (e->type_id() == TypeId::Float) {
        // We don't have a clean Float -> double helper; round-trip via str.
        return std::stod(e->str());
    }
    return std::nullopt;
}

[[nodiscard]] bool close_enough(double a, double b, double tol) {
    return std::fabs(a - b) < tol * std::max(1.0, std::fabs(a));
}

// Approximate value of a SymPP Expr as a double. Returns nullopt if the
// expression isn't a numeric/named-symbol leaf the table covers.
[[nodiscard]] std::optional<double> evalf_to_double(const Expr& e, int dps) {
    if (auto d = as_double(e); d) return d;
    if (e == S::Pi()) return 3.14159265358979323846;
    if (e == S::E()) return 2.71828182845904523536;
    if (e == S::EulerGamma()) return 0.577215664901532860606;
    if (e == S::Catalan()) return 0.915965594177219015055;
    Expr v = evalf(e, dps);
    if (v->type_id() == TypeId::Float) return std::stod(v->str());
    return std::nullopt;
}

}  // namespace

Expr nsimplify(const Expr& e, int dps) {
    // Only operates on numeric-valued expressions; pass non-numeric through.
    auto v_opt = evalf_to_double(e, dps);
    if (!v_opt) return e;
    double v = *v_opt;

    const double tol = 1e-12;

    // Try simple Rational p/q with small q (Stern-Brocot-style).
    for (long q = 1; q <= 1000; ++q) {
        long p = static_cast<long>(std::lround(v * static_cast<double>(q)));
        if (close_enough(v,
                          static_cast<double>(p) / static_cast<double>(q),
                          tol)) {
            return rational(p, q);
        }
    }

    // pi, e, EulerGamma, Catalan multiples (k/q) for small q, k.
    auto try_constant = [&](double cval, const Expr& sym) -> std::optional<Expr> {
        for (long q = 1; q <= 60; ++q) {
            long p = static_cast<long>(std::lround(v / cval * static_cast<double>(q)));
            if (close_enough(v,
                              cval * static_cast<double>(p) / static_cast<double>(q),
                              tol)) {
                if (q == 1 && p == 1) return sym;
                if (q == 1) return mul(integer(p), sym);
                return mul(rational(p, q), sym);
            }
        }
        return std::nullopt;
    };
    if (auto r = try_constant(3.14159265358979323846, S::Pi())) return *r;
    if (auto r = try_constant(2.71828182845904523536, S::E())) return *r;

    // sqrt of small rationals.
    for (long q = 1; q <= 50; ++q) {
        for (long p = 1; p <= 50; ++p) {
            double candidate = std::sqrt(static_cast<double>(p)
                                          / static_cast<double>(q));
            if (close_enough(v, candidate, tol)) {
                if (q == 1) return sqrt(integer(p));
                return sqrt(rational(p, q));
            }
        }
    }

    // No clean form found; return original.
    return e;
}

}  // namespace sympp
