#include <sympp/integrals/integrate.hpp>

#include <optional>
#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/core/undefined_function.hpp>
#include <sympp/calculus/diff.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/simplify/simplify.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/special.hpp>
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/polys/poly.hpp>

namespace sympp {

namespace {

// Helper: does `e` depend on `var`? Uses has() from traversal.
[[nodiscard]] bool depends_on(const Expr& e, const Expr& var) noexcept {
    return has(e, var);
}

// Decompose `e` as (a, b) such that e == a*var + b, with a, b independent
// of var. Returns nullopt if e isn't affine in var.
[[nodiscard]] std::optional<std::pair<Expr, Expr>>
as_affine(const Expr& e, const Expr& var) {
    if (!e) return std::nullopt;
    if (!has(e, var)) return std::pair{S::Zero(), e};
    if (e == var) return std::pair{S::One(), S::Zero()};
    if (e->type_id() == TypeId::Mul) {
        // Split into constant and var-dependent factors. Any factor of
        // var must itself be affine, and exactly one such factor is
        // permitted (multiplying two affines in var produces a quadratic).
        std::vector<Expr> const_factors;
        std::vector<Expr> var_factors;
        for (const auto& a : e->args()) {
            if (has(a, var)) var_factors.push_back(a);
            else const_factors.push_back(a);
        }
        if (var_factors.size() != 1) return std::nullopt;
        auto sub = as_affine(var_factors[0], var);
        if (!sub) return std::nullopt;
        Expr scale = mul(const_factors);
        return std::pair{mul(scale, sub->first), mul(scale, sub->second)};
    }
    if (e->type_id() == TypeId::Add) {
        Expr a_acc = S::Zero();
        Expr b_acc = S::Zero();
        for (const auto& term : e->args()) {
            auto sub = as_affine(term, var);
            if (!sub) return std::nullopt;
            a_acc = add(a_acc, sub->first);
            b_acc = add(b_acc, sub->second);
        }
        return std::pair{a_acc, b_acc};
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<Expr> integrate_term(const Expr& term, const Expr& var);

// Forward-decl of the public entry — integrate_term needs to recurse
// through the full dispatch (Add linearity, trig, rational) when it
// pulls a constant out of a Mul whose remainder is an Add.
}  // namespace
Expr integrate(const Expr& expr, const Expr& var);
namespace {

// Try the elementary table on a single term. Returns std::nullopt if the
// term is outside the table.
[[nodiscard]] std::optional<Expr> integrate_term(const Expr& term, const Expr& var) {
    // Constant w.r.t. var: ∫ c dx = c*x
    if (!depends_on(term, var)) {
        return mul(term, var);
    }

    // var itself: ∫ x dx = x^2/2
    if (term == var) {
        return mul(rational(1, 2), pow(var, integer(2)));
    }

    // (a*x + b)^n with n constant: ∫(ax+b)^n dx = (ax+b)^(n+1) / (a*(n+1))
    // and n == -1 → log(ax+b)/a.
    if (term->type_id() == TypeId::Pow) {
        const auto& base = term->args()[0];
        const auto& exp = term->args()[1];
        if (!depends_on(exp, var)) {
            auto aff = as_affine(base, var);
            if (aff && !(aff->first == S::Zero())) {
                if (exp == S::NegativeOne()) {
                    return log(base) / aff->first;
                }
                Expr n_plus_1 = add(exp, S::One());
                return pow(base, n_plus_1) / (aff->first * n_plus_1);
            }
        }
        // Reciprocal-square trig of an affine argument:
        //   ∫ 1/cos²(ax+b) dx =  sin(ax+b)/(a·cos(ax+b))   (= tan/a)
        //   ∫ 1/sin²(ax+b) dx = -cos(ax+b)/(a·sin(ax+b))   (= -cot/a)
        if (exp == integer(-2) && base->type_id() == TypeId::Function) {
            const auto& bfn = static_cast<const Function&>(*base);
            if (bfn.args().size() == 1) {
                const auto& inner = bfn.args()[0];
                auto inner_aff = as_affine(inner, var);
                if (inner_aff && !(inner_aff->first == S::Zero())) {
                    const Expr& a = inner_aff->first;
                    if (bfn.function_id() == FunctionId::Cos) {
                        return sin(inner) / mul(a, cos(inner));
                    }
                    if (bfn.function_id() == FunctionId::Sin) {
                        return mul(S::NegativeOne(), cos(inner)) / mul(a, sin(inner));
                    }
                    // Hyperbolic counterparts:
                    //   ∫ 1/cosh²(ax+b) dx =  tanh(ax+b)/a   (sech² → tanh)
                    //   ∫ 1/sinh²(ax+b) dx = -cosh(ax+b)/(a·sinh(ax+b))  (−coth/a)
                    if (bfn.function_id() == FunctionId::Cosh) {
                        return tanh(inner) / a;
                    }
                    if (bfn.function_id() == FunctionId::Sinh) {
                        return mul(S::NegativeOne(), cosh(inner)) / mul(a, sinh(inner));
                    }
                }
            }
        }
    }

    // 1/x: ∫ 1/x dx = log(x). 1/x is Pow(x, -1) — caught above.

    // sin(ax+b), cos(ax+b), exp(ax+b) where inner is affine in var.
    if (term->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*term);
        if (fn.args().size() == 1) {
            const auto& inner = fn.args()[0];
            auto aff = as_affine(inner, var);
            if (aff && !(aff->first == S::Zero())) {
                const Expr& a = aff->first;
                switch (fn.function_id()) {
                    case FunctionId::Sin:
                        // ∫sin(ax+b) dx = -cos(ax+b)/a
                        return mul(S::NegativeOne(), cos(inner)) / a;
                    case FunctionId::Cos:
                        return sin(inner) / a;
                    case FunctionId::Tan:
                        // ∫tan(ax+b) dx = -log(cos(ax+b))/a
                        return mul(S::NegativeOne(), log(cos(inner))) / a;
                    case FunctionId::Cot:
                        // ∫cot(ax+b) dx = log(sin(ax+b))/a
                        return log(sin(inner)) / a;
                    case FunctionId::Sec: {
                        // ∫sec(ax+b) dx = (log(sin+1) − log(sin−1)) / (2a),
                        // the half-angle log form SymPy prints (≡ log|sec+tan|).
                        Expr s = sin(inner);
                        return (log(s + S::One()) - log(s - S::One()))
                               / (integer(2) * a);
                    }
                    case FunctionId::Csc: {
                        // ∫csc(ax+b) dx = (log(cos−1) − log(cos+1)) / (2a).
                        Expr c = cos(inner);
                        return (log(c - S::One()) - log(c + S::One()))
                               / (integer(2) * a);
                    }
                    case FunctionId::Exp:
                        return exp(inner) / a;
                    case FunctionId::Sinh:
                        // ∫sinh(ax+b) dx = cosh(ax+b)/a
                        return cosh(inner) / a;
                    case FunctionId::Cosh:
                        // ∫cosh(ax+b) dx = sinh(ax+b)/a
                        return sinh(inner) / a;
                    case FunctionId::Tanh:
                        // ∫tanh(ax+b) dx = log(cosh(ax+b))/a
                        return log(cosh(inner)) / a;
                    case FunctionId::Coth:
                        // ∫coth(ax+b) dx = log(sinh(ax+b))/a
                        return log(sinh(inner)) / a;
                    case FunctionId::Sech:
                        // ∫sech(ax+b) dx = atan(sinh(ax+b))/a (Gudermannian).
                        return atan(sinh(inner)) / a;
                    case FunctionId::Csch:
                        // ∫csch(ax+b) dx = log(tanh((ax+b)/2))/a.
                        return log(tanh(inner / integer(2))) / a;
                    default:
                        break;
                }
            }
        }
    }

    // Mul: pull out factors that don't depend on var, then recurse on the rest.
    if (term->type_id() == TypeId::Mul) {
        std::vector<Expr> constant_factors;
        std::vector<Expr> var_factors;
        for (const auto& a : term->args()) {
            if (depends_on(a, var)) {
                var_factors.push_back(a);
            } else {
                constant_factors.push_back(a);
            }
        }
        if (var_factors.empty()) {
            // Entire term is constant — handled at the top of this function;
            // shouldn't reach here.
            return mul(mul(std::move(constant_factors)), var);
        }
        // If we didn't pull anything constant out, recursing on
        // mul(var_factors) would just be the same term — bail to avoid
        // infinite recursion. Higher-level strategies (try_rational,
        // integration by parts, etc.) get a chance instead.
        if (constant_factors.empty()) {
            return std::nullopt;
        }
        // Try integrating the remaining var-dependent product. Use the
        // full integrate() dispatch (not just integrate_term) so an Add
        // inner gets linearity, and so trig / rational strategies still
        // apply.
        Expr inner = mul(std::move(var_factors));
        Expr sub = integrate(inner, var);
        // integrate() returns an Integral(_, _) marker on failure (an
        // UndefinedFunction named "Integral"); propagate as nullopt so
        // higher strategies can try other approaches.
        if (sub->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*sub);
            if (fn.name() == "Integral") return std::nullopt;
        }
        return mul(mul(std::move(constant_factors)), sub);
    }

    // Outside the table.
    return std::nullopt;
}

// Heurisch — u-substitution recognition for non-linear inner arguments.
// Detects the chain-rule-reverse pattern c·g'(x)·f(g(x)) where f is one
// of {sin, cos, exp, 1/u}, and emits c·F(g(x)) where F is the
// antiderivative of f. Catches standard textbook integrals like
// ∫2x·exp(x²) → exp(x²) that the table-based path can't recognize.
[[nodiscard]] std::optional<Expr> try_heurisch(const Expr& expr, const Expr& var);

// Integration by parts: ∫u dv = uv - ∫v du. Handles two patterns:
//   * standalone log(affine) → x*log(ax+b) + (b/a)*log(ax+b) - x
//   * Mul where one factor is a single sin/cos/exp of affine and the
//     remaining factor depends on var (typically a polynomial):
//     u = remaining; dv = target dx. Recurses on ∫v du via integrate.
[[nodiscard]] std::optional<Expr> try_integration_by_parts(
    const Expr& expr, const Expr& var);

// Trig identity rewrites:
//   sin²(u) → (1 - cos(2u))/2
//   cos²(u) → (1 + cos(2u))/2
//   sin(p)cos(q) → (sin(p+q) + sin(p-q))/2
// After rewriting, recurse via integrate so the existing affine-argument
// rules pick up the resulting sin(linear) / cos(linear) terms.
[[nodiscard]] std::optional<Expr> try_trig_reduction(const Expr& expr, const Expr& var);

// ∫ sin(g)^m · cos(g)^n dx for an affine argument g and non-negative integer
// powers (at least one ≥ 1). Odd power → u = cos/sin substitution into a
// polynomial; both even → half-angle reduction, recursing via integrate.
[[nodiscard]] std::optional<Expr> try_trig_power(const Expr& expr, const Expr& var);

// ∫ tan(g)^n dx (n ≥ 2 integer, g affine) via the reduction
// ∫tanⁿ = tan^(n-1)/((n-1)·g') − ∫tan^(n-2), recursing through integrate.
[[nodiscard]] std::optional<Expr> try_tan_power(const Expr& expr, const Expr& var);

// ∫ tanh(g)^n / coth(g)^n dx (n ≥ 2 integer, g affine) — the hyperbolic analogue
// of try_tan_power. Both share the reduction ∫fⁿ = ∫f^(n-2) − f^(n-1)/((n-1)·g')
// (tanh from tanh²=1−sech², coth from coth²=1+csch²), recursing through integrate
// to the n=1 table case (∫tanh = log(cosh)/g', ∫coth = log(sinh)/g').
[[nodiscard]] std::optional<Expr> try_tanh_coth_power(const Expr& expr, const Expr& var);

// ∫ sec(g)^n / csc(g)^n dx (n ≥ 3 integer, g affine) via the by-parts reduction
// ∫secⁿ =  sec^(n-2)·tan/((n-1)·g') + (n-2)/(n-1)·∫sec^(n-2)
// ∫cscⁿ = −csc^(n-2)·cot/((n-1)·g') + (n-2)/(n-1)·∫csc^(n-2). n=1 is the table
// case (INT-24) and n=2 the trig-reduction square (INT-25).
[[nodiscard]] std::optional<Expr> try_sec_csc_power(const Expr& expr, const Expr& var);

// ∫ sech(g)^n / csch(g)^n dx (n ≥ 3 integer, g affine) — the hyperbolic analogue
// of try_sec_csc_power. The Pythagorean sign differs (coth²−csch²=1 vs
// csc²−cot²=1), so the csch rest term is subtracted:
// ∫sechⁿ =  sech^(n-2)·tanh/((n-1)·g') + (n-2)/(n-1)·∫sech^(n-2)
// ∫cschⁿ = −csch^(n-2)·coth/((n-1)·g') − (n-2)/(n-1)·∫csch^(n-2). n=1/n=2 are
// the table / square cases (INT-26).
[[nodiscard]] std::optional<Expr> try_sech_csch_power(const Expr& expr, const Expr& var);

// ∫ sinh(g)^m · cosh(g)^n dx — the hyperbolic analogue of try_trig_power,
// using cosh²−sinh²=1 and the half-angle forms of cosh(2g).
[[nodiscard]] std::optional<Expr> try_hyperbolic_power(const Expr& expr, const Expr& var);

// Try to integrate `expr` as a rational function in var. Uses
// polynomial division for the improper part and apart() for the proper
// part. Returns nullopt when:
//   * expr has no top-level denominator structure, or
//   * num / den isn't pure-polynomial in var, or
//   * apart can't decompose the proper remainder (e.g. repeated roots
//     beyond what apart handles in this build).
[[nodiscard]] std::optional<Expr> try_rational(const Expr& expr, const Expr& var);

// Product mixing a hyperbolic sinh/cosh(affine) with a sin/cos/exp(affine)
// factor: rewrite the hyperbolics to exponentials (sinh g = (e^g − e^−g)/2,
// cosh g = (e^g + e^−g)/2), expand, and integrate term by term — each term is a
// c·e^(·)·sin/cos(·) (cyclic closed form) or a pure exponential. Closes
// ∫sin·sinh, ∫cos·cosh, ∫sin·cosh, ∫e^x·sinh, etc.
[[nodiscard]] std::optional<Expr> try_hyperbolic_to_exp(const Expr& expr,
                                                        const Expr& var);

// ∫ 1/(a·x²+b·x+c)ⁿ dx for integer n ≥ 2 (a constant numerator) via the standard
// reduction Iₙ = a·x/(2(n−1)c·Qⁿ⁻¹) + a(2n−3)/(2(n−1)c)·Iₙ₋₁, recursing through
// integrate down to I₁ (try_arctan_quadratic / try_rational). A linear term is
// removed first by completing the square (as in try_sqrt_quadratic).
[[nodiscard]] std::optional<Expr> try_quadratic_power(const Expr& expr, const Expr& var);

// Weierstrass (half-angle) substitution t = tan(var/2): a rational function of
// sin(var)/cos(var)/tan(var)/cot(var)/sec(var)/csc(var) becomes a rational
// function of t, which try_rational closes. Last-resort fallback — its tan(x/2)
// output is uglier than the dedicated trig integrators, so it runs only after
// everything else. Requires the bare argument `var` (affine arguments bail).
[[nodiscard]] std::optional<Expr> try_weierstrass(const Expr& expr, const Expr& var);

// ∫ 1/(a·x² + b·x + c) dx for an irreducible quadratic (discriminant
// 4ac − b² > 0): the arctangent rule  2·atan((2ax+b)/√D)/√D. Reducible
// denominators (real roots) are left to try_rational, which splits them into
// logs. Returns nullopt unless the denominator is a degree-2 polynomial in
// var with numeric coefficients and positive discriminant.
[[nodiscard]] std::optional<Expr> try_arctan_quadratic(const Expr& expr, const Expr& var);

// ∫ (p·x + q)/(a·x² + b·x + c) dx for an irreducible quadratic denominator:
// split the numerator into the part proportional to the denominator's
// derivative (→ log) plus a constant remainder (→ the arctangent rule). E.g.
// ∫(2x+3)/(x²+1) = log(x²+1) + 3·atan(x). Reducible / repeated-root quadratics
// are left to try_rational; constant numerators to try_arctan_quadratic.
[[nodiscard]] std::optional<Expr> try_linear_over_quadratic(const Expr& expr,
                                                            const Expr& var);

// ∫ 1/√(a·x² + c) dx for a pure quadratic under the root (no linear term)
// with c > 0:  a > 0 → asinh(x·√(a/c))/√a ;  a < 0 → asin(x·√(−a/c))/√(−a).
// A linear term, or c ≤ 0 (the acosh / log forms), is out of scope. Returns
// nullopt unless the radicand is a degree-2 polynomial with rational
// coefficients, zero linear term, and positive constant term.
[[nodiscard]] std::optional<Expr> try_sqrt_quadratic(const Expr& expr, const Expr& var);

// ∫ x/√(a·x² + c) dx = √(a·x²+c)/a (pure quadratic radicand, no linear term).
// The numerator must be the bare var (times a constant); a constant numerator is
// try_sqrt_quadratic's job. Closes ∫x/√(1−x²) etc. — the x·f' term that the
// inverse-trig/hyperbolic by-parts integrals reduce to.
[[nodiscard]] std::optional<Expr> try_x_over_sqrt_quadratic(const Expr& expr, const Expr& var);

// ∫ exp(c·x²) dx for a concrete negative c (a real Gaussian) → the error
// function: with a = −c > 0, √π·erf(√a·x)/(2√a). Returns nullopt unless the
// exponent is a pure quadratic (no linear/constant term) with negative
// rational leading coefficient.
[[nodiscard]] std::optional<Expr> try_gaussian(const Expr& expr, const Expr& var);

// ∫ sin(c·x)/x = Si(c·x), ∫ cos(c·x)/x = Ci(c·x), ∫ exp(c·x)/x = Ei(c·x) — the
// special-integral functions, for a monomial argument c·x (no constant term).
[[nodiscard]] std::optional<Expr> try_expint_integral(const Expr& expr, const Expr& var);

// ∫ P(x)·(a·x+b)^r dx for a polynomial P, an affine base a·x+b, and a
// non-integer rational exponent r — e.g. ∫x·√(x+1). Substitute u = a·x+b so the
// integrand becomes Σ cₖ·u^(k+r), which integrates term-by-term. Returns nullopt
// unless exactly one (affine)^(rational) factor is present and the rest is
// polynomial in var.
[[nodiscard]] std::optional<Expr> try_algebraic_linear_sub(const Expr& expr, const Expr& var);

}  // namespace

// Recursion-depth backstop. Integration by parts recurses through
// integrate(); cyclic integrands (e.g. exp(x)*sin(x), where the original
// integral reappears) would otherwise recurse without bound and overflow
// the stack. When the limit is hit we bail to the unevaluated Integral
// marker instead of crashing. The limit is far above any legitimate
// closed-form nesting depth.
constexpr int kMaxIntegrateDepth = 64;
thread_local int g_integrate_depth = 0;

struct IntegrateDepthGuard {
    IntegrateDepthGuard() { ++g_integrate_depth; }
    ~IntegrateDepthGuard() { --g_integrate_depth; }
};

Expr integrate(const Expr& expr, const Expr& var) {
    if (!expr || !var) return S::Zero();

    if (g_integrate_depth >= kMaxIntegrateDepth) {
        return function_symbol("Integral")(expr, var);
    }
    IntegrateDepthGuard depth_guard;

    // Linearity: split Add into per-term integrals.
    if (expr->type_id() == TypeId::Add) {
        std::vector<Expr> terms;
        terms.reserve(expr->args().size());
        for (const auto& a : expr->args()) {
            terms.push_back(integrate(a, var));
        }
        return add(std::move(terms));
    }

    if (auto r = integrate_term(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_trig_reduction(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_trig_power(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_tan_power(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_tanh_coth_power(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_sec_csc_power(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_sech_csch_power(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_hyperbolic_power(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_rational(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_arctan_quadratic(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_linear_over_quadratic(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_quadratic_power(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_sqrt_quadratic(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_x_over_sqrt_quadratic(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_algebraic_linear_sub(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_expint_integral(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_hyperbolic_to_exp(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_heurisch(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_gaussian(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_integration_by_parts(expr, var); r.has_value()) {
        return *r;
    }
    if (auto r = try_weierstrass(expr, var); r.has_value()) {
        return *r;
    }

    // Outside the closed-form table — return an unevaluated marker. Use an
    // UndefinedFunction named "Integral" so callers can detect it.
    return function_symbol("Integral")(expr, var);
}

namespace {

// Detect the public `Integral(_, _)` failure marker.
[[nodiscard]] bool is_integral_marker(const Expr& e) {
    if (!e || e->type_id() != TypeId::Function) return false;
    const auto& fn = static_cast<const Function&>(*e);
    return fn.name() == "Integral";
}

// Try to interpret expr as c * g'(x) * f(g(x)) for some non-trivial
// inner expression g(x) and outer pattern f ∈ {sin, cos, exp, log,
// reciprocal}. Strategy:
//   1. Scan the expression for candidate inner forms — function args
//      and Pow bases that depend on var.
//   2. For each candidate g, compute g' and divide expr by g'. Simplify.
//   3. Check that the quotient is independent of var aside from the g
//      occurrences themselves — i.e. it depends on var only through g.
//      Operationally: substitute g → fresh symbol u in the quotient
//      and verify the result is var-free.
//   4. If the quotient as a function of u admits a closed-form
//      antiderivative (try_integrate against the fresh symbol u via the
//      table), substitute back and return.
std::optional<Expr> try_heurisch(const Expr& expr, const Expr& var) {
    // Collect candidate inner expressions.
    std::vector<Expr> candidates;
    auto add_candidate = [&](const Expr& cand) {
        if (!cand || !depends_on(cand, var)) return;
        if (cand == var) return;  // linear-sub already handled
        for (const auto& c : candidates) if (c == cand) return;
        candidates.push_back(cand);
    };
    auto walk = [&](auto&& self, const Expr& e) -> void {
        if (!e) return;
        if (e->type_id() == TypeId::Function) {
            // The function application itself is a candidate inner g — e.g.
            // g = log(x) for ∫1/(x·log(x)). Without this, a function buried as
            // a Mul factor would never be considered (only its args were).
            add_candidate(e);
            for (const auto& a : e->args()) {
                add_candidate(a);
                self(self, a);
            }
        }
        if (e->type_id() == TypeId::Pow) {
            add_candidate(e->args()[0]);
            self(self, e->args()[0]);
            self(self, e->args()[1]);
        }
        if (e->type_id() == TypeId::Add || e->type_id() == TypeId::Mul) {
            for (const auto& a : e->args()) self(self, a);
        }
    };
    walk(walk, expr);

    // Sort candidates by descending node count — try larger / outer
    // candidates first.
    auto node_count = [](const Expr& e) {
        int n = 0;
        auto rec = [&](auto&& self, const Expr& x) -> void {
            if (!x) return;
            ++n;
            for (const auto& a : x->args()) self(self, a);
        };
        rec(rec, e);
        return n;
    };
    std::sort(candidates.begin(), candidates.end(),
              [&](const Expr& a, const Expr& b) {
                  return node_count(a) > node_count(b);
              });

    auto u = symbol("__heurisch_u");
    for (const auto& g : candidates) {
        Expr gp = simplify(diff(g, var));
        if (gp == S::Zero()) continue;
        // q = expr / g'. Substitute g → u first so the symbolic structure
        // collapses (e.g. exp(x²) → exp(u)). Apply expand_power_base so
        // patterns like (2x)^(-1) split into 2^(-1) * x^(-1), letting
        // canonical Mul base collection cancel x * x^(-1) cleanly.
        Expr q = expr / gp;
        Expr q_sub = subs(q, g, u);
        q_sub = simplify(expand_power_base(q_sub));
        if (depends_on(q_sub, var)) continue;
        // Try to integrate q_sub against u via the table only — calling
        // the full integrate() here would loop back into try_heurisch.
        // Linearity over Add is handled inline; the rest is integrate_term.
        std::optional<Expr> integrated_opt;
        if (q_sub->type_id() == TypeId::Add) {
            std::vector<Expr> parts;
            parts.reserve(q_sub->args().size());
            bool ok = true;
            for (const auto& term : q_sub->args()) {
                auto sub_int = integrate_term(term, u);
                if (!sub_int) { ok = false; break; }
                parts.push_back(*sub_int);
            }
            if (ok) integrated_opt = add(std::move(parts));
        } else {
            integrated_opt = integrate_term(q_sub, u);
        }
        // When the table can't close it, q_sub may still be a rational function
        // of u — e.g. g = exp(x) turns ∫1/(1+exp(x)) into ∫1/(u·(1+u)). Hand
        // that to try_rational (partial fractions), which decomposes into
        // strictly simpler pieces and so terminates; the depth guard backstops
        // any pathological recursion through its internal integrate() calls.
        if (!integrated_opt) integrated_opt = try_rational(q_sub, u);
        // Irreducible-quadratic fallback: a bare or numerically-scaled
        // 1/(a·u²+b·u+c) (or (linear)/quadratic) — e.g. g = exp(x) turns
        // ∫1/(eˣ+e⁻ˣ) into ∫1/(u²+1), and g = x² turns ∫x/(x⁴+1) into
        // ∫1/(2(u²+1)). The table and try_rational don't close these; the
        // dedicated quadratic helpers do, once a leading numeric factor is
        // pulled aside.
        if (!integrated_opt) {
            Expr coeff = S::One();
            Expr core = q_sub;
            if (q_sub->type_id() == TypeId::Mul) {
                std::vector<Expr> num_consts, rest;
                for (const auto& f : q_sub->args()) {
                    if (f->type_id() == TypeId::Integer
                        || f->type_id() == TypeId::Rational) {
                        num_consts.push_back(f);
                    } else {
                        rest.push_back(f);
                    }
                }
                if (!num_consts.empty() && !rest.empty()) {
                    coeff = mul(num_consts);
                    core = mul(rest);
                }
            }
            auto qr = try_arctan_quadratic(core, u);
            if (!qr) qr = try_linear_over_quadratic(core, u);
            if (qr) integrated_opt = mul(coeff, qr.value());
        }
        if (!integrated_opt) continue;
        Expr integrated = *integrated_opt;
        if (is_integral_marker(integrated)) continue;
        // Walk the result to ensure it didn't smuggle an Integral marker.
        bool has_marker = false;
        auto rec_check = [&](auto&& self, const Expr& e) -> void {
            if (!e) return;
            if (is_integral_marker(e)) { has_marker = true; return; }
            for (const auto& a : e->args()) {
                if (has_marker) return;
                self(self, a);
            }
        };
        rec_check(rec_check, integrated);
        if (has_marker) continue;
        // Substitute u back to g.
        Expr result = subs(integrated, u, g);
        return simplify(result);
    }
    return std::nullopt;
}

// True if `e` is a polynomial in `var` (only non-negative integer powers of
// var). Such a `u` differentiates to zero in finitely many steps, so using it
// as the `u` factor of integration by parts terminates. A factor like x^(-1)
// is NOT a polynomial — its derivatives x^(-2), x^(-3), … grow without bound,
// so by parts would recurse forever (e.g. ∫exp(x)/x, which is non-elementary).
[[nodiscard]] bool is_polynomial_in(const Expr& e, const Expr& var) {
    if (!depends_on(e, var)) return true;
    if (e == var) return true;
    switch (e->type_id()) {
        case TypeId::Add:
        case TypeId::Mul:
            for (const auto& a : e->args()) {
                if (!is_polynomial_in(a, var)) return false;
            }
            return true;
        case TypeId::Pow: {
            const Expr& exp = e->args()[1];
            if (exp->type_id() != TypeId::Integer) return false;
            if (static_cast<const Integer&>(*exp).is_negative()) return false;
            return is_polynomial_in(e->args()[0], var);
        }
        default:
            return false;  // a function of var, a fractional power, etc.
    }
}

// True if `f` is log(affine-in-var) or a positive-integer power of one.
// Used to drive integration by parts over log and log^n factors.
[[nodiscard]] bool is_log_or_log_power(const Expr& f, const Expr& var) {
    auto is_log_affine = [&](const Expr& g) {
        if (g->type_id() != TypeId::Function) return false;
        const auto& fn = static_cast<const Function&>(*g);
        return fn.function_id() == FunctionId::Log && fn.args().size() == 1
               && as_affine(fn.args()[0], var).has_value();
    };
    if (is_log_affine(f)) return true;
    if (f->type_id() == TypeId::Pow) {
        const auto& e = f->args()[1];
        return e->type_id() == TypeId::Integer
               && static_cast<const Integer&>(*e).is_positive()
               && is_log_affine(f->args()[0]);
    }
    return false;
}

// Functions f for which ∫f dx = x·f − ∫x·f' closes because f' is elementary and
// x·f' (or, more generally, a polynomial × f') integrates. The special-integral
// functions (erf, Si, …) and the inverse trig / hyperbolic functions qualify; a
// function left with the default 0-derivative would yield a bogus x·f, so the
// set is an explicit whitelist. Drives both the standalone and the
// polynomial × f integration-by-parts branches.
[[nodiscard]] bool is_by_parts_fn(FunctionId id) {
    return id == FunctionId::Erf || id == FunctionId::Erfc
           || id == FunctionId::Erfi || id == FunctionId::Si
           || id == FunctionId::Ci || id == FunctionId::Ei
           || id == FunctionId::Shi || id == FunctionId::Chi
           || id == FunctionId::Asin || id == FunctionId::Acos
           || id == FunctionId::Atan || id == FunctionId::Acot
           || id == FunctionId::Asinh || id == FunctionId::Acosh
           || id == FunctionId::Atanh;
}

std::optional<Expr> try_integration_by_parts(const Expr& expr, const Expr& var) {
    // Standalone log(affine).
    if (expr->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*expr);
        if (fn.function_id() == FunctionId::Log && fn.args().size() == 1) {
            auto aff = as_affine(fn.args()[0], var);
            if (aff && !(aff->first == S::Zero())) {
                const Expr& inner = fn.args()[0];
                const Expr& a = aff->first;
                const Expr& b = aff->second;
                // ∫log(ax+b) dx = x*log(ax+b) + (b/a)*log(ax+b) - x
                return var * log(inner) + (b / a) * log(inner) - var;
            }
        }
    }

    // Standalone special-integral function f(affine): by parts with u = f,
    // dv = dx, v = x → ∫f dx = x·f − ∫x·f'. Closes because x·f' is elementary
    // for these (erf' = 2e^(−x²)/√π so x·erf' integrates, Si' = sin(x)/x so
    // x·Si' = sin(x), etc.). Whitelisted to functions with a correct derivative
    // — a function left with the default 0-derivative would yield a bogus x·f.
    if (expr->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*expr);
        const FunctionId id = fn.function_id();
        if (is_by_parts_fn(id) && fn.args().size() == 1) {
            auto aff = as_affine(fn.args()[0], var);
            if (aff && !(aff->first == S::Zero())) {
                Expr remaining = integrate(var * diff(expr, var), var);
                if (!is_integral_marker(remaining)) {
                    return expand(var * expr - remaining);
                }
            }
        }
    }

    // Standalone log(affine)^n (n ≥ 2): by parts with u = log^n, dv = dx, v = x.
    //   ∫ log^n dx = x·log^n − ∫ x·(log^n)' dx.
    // For log(c·x) the remaining x·(log^n)' collapses to n·log^(n-1), so it
    // recurses down to ∫log = the case above. The marker guard bails on
    // anything that does not reduce, so a stray case never loops or emits a
    // wrong closed form.
    if (expr->type_id() == TypeId::Pow && is_log_or_log_power(expr, var)) {
        Expr remaining = integrate(var * diff(expr, var), var);
        if (!is_integral_marker(remaining)) {
            return expand(var * expr - remaining);
        }
    }

    // Cyclic case: ∫ c·e^(a·x+·)·sin/cos(g·x+·) dx. Generic by-parts would
    // recurse exp·sin → exp·cos → exp·sin … forever, so solve it in closed
    // form. With E = e^(arg_e), S = sin(arg_t), C = cos(arg_t), and a, g the
    // x-coefficients of arg_e, arg_t:
    //   ∫ E·S dx = E·(a·S − g·C)/(a²+g²)
    //   ∫ E·C dx = E·(a·C + g·S)/(a²+g²)
    if (expr->type_id() == TypeId::Mul) {
        Expr exp_factor, trig_factor;
        std::vector<Expr> consts;
        bool ok = true;
        for (const auto& f : expr->args()) {
            if (!exp_factor && f->type_id() == TypeId::Function) {
                const auto& fn = static_cast<const Function&>(*f);
                if (fn.function_id() == FunctionId::Exp && fn.args().size() == 1
                    && as_affine(fn.args()[0], var)) {
                    exp_factor = f;
                    continue;
                }
            }
            if (!trig_factor && f->type_id() == TypeId::Function) {
                const auto& fn = static_cast<const Function&>(*f);
                if ((fn.function_id() == FunctionId::Sin
                     || fn.function_id() == FunctionId::Cos)
                    && fn.args().size() == 1 && as_affine(fn.args()[0], var)) {
                    trig_factor = f;
                    continue;
                }
            }
            // Any leftover factor must be constant w.r.t. var.
            if (depends_on(f, var)) { ok = false; break; }
            consts.push_back(f);
        }
        if (ok && exp_factor && trig_factor) {
            const auto& tfn = static_cast<const Function&>(*trig_factor);
            Expr arg_e = exp_factor->args()[0];
            Expr arg_t = trig_factor->args()[0];
            Expr a = as_affine(arg_e, var)->first;  // x-coefficient of exp arg
            Expr g = as_affine(arg_t, var)->first;  // x-coefficient of trig arg
            Expr denom = a * a + g * g;
            Expr S = sin(arg_t);
            Expr C = cos(arg_t);
            Expr core = (tfn.function_id() == FunctionId::Sin)
                            ? exp_factor * (a * S - g * C) / denom
                            : exp_factor * (a * C + g * S) / denom;
            if (!consts.empty()) core = mul(consts) * core;
            return simplify(core);
        }
    }

    // Polynomial × e^(a·x+·) × sin/cos(g·x+·): by parts with u = the polynomial
    // and dv = e·trig, whose antiderivative is the cyclic closed form above.
    // Differentiating u lowers its degree by one each step, so the recursion
    // bottoms out at the bare cyclic case (constant u). Without this, the cyclic
    // block bails (the polynomial factor makes its `rest` non-constant) and the
    // single-function by-parts below bails too (u = x·sin(x) isn't polynomial).
    if (expr->type_id() == TypeId::Mul) {
        Expr exp_factor, trig_factor;
        std::vector<Expr> poly_factors;
        for (const auto& f : expr->args()) {
            if (!exp_factor && f->type_id() == TypeId::Function) {
                const auto& fn = static_cast<const Function&>(*f);
                if (fn.function_id() == FunctionId::Exp && fn.args().size() == 1
                    && as_affine(fn.args()[0], var)) {
                    exp_factor = f;
                    continue;
                }
            }
            if (!trig_factor && f->type_id() == TypeId::Function) {
                const auto& fn = static_cast<const Function&>(*f);
                if ((fn.function_id() == FunctionId::Sin
                     || fn.function_id() == FunctionId::Cos)
                    && fn.args().size() == 1 && as_affine(fn.args()[0], var)) {
                    trig_factor = f;
                    continue;
                }
            }
            poly_factors.push_back(f);
        }
        Expr u = poly_factors.empty() ? S::One() : mul(poly_factors);
        if (exp_factor && trig_factor && depends_on(u, var)
            && is_polynomial_in(u, var)) {
            Expr v = integrate(exp_factor * trig_factor, var);  // cyclic form
            if (!is_integral_marker(v)) {
                // expand so v·u' = poly·e·sin + poly·e·cos splits over the Add,
                // letting integrate() recurse per term with deg(u) reduced.
                Expr remaining = integrate(expand(v * diff(u, var)), var);
                if (!is_integral_marker(remaining)) {
                    return simplify(u * v - remaining);
                }
            }
        }
    }

    // Polynomial × log(affine), by parts with u = log(ax+b), dv = rest dx:
    //   ∫ rest·log(ax+b) dx = log(ax+b)·∫rest − ∫(∫rest)·a/(ax+b) dx.
    // The remaining integral is rational (∫rest is polynomial, du = a/(ax+b)),
    // so try_rational closes it. The marker/depth guards bail on anything that
    // does not reduce, so a stray log in `rest` can't loop.
    if (expr->type_id() == TypeId::Mul) {
        Expr log_factor;
        std::vector<Expr> rest_factors;
        for (const auto& f : expr->args()) {
            if (!log_factor && is_log_or_log_power(f, var)) {
                log_factor = f;
                continue;
            }
            rest_factors.push_back(f);
        }
        if (log_factor && !rest_factors.empty()) {
            Expr rest = mul(rest_factors);
            if (depends_on(rest, var)) {
                Expr v = integrate(rest, var);
                if (!is_integral_marker(v)) {
                    Expr remaining = integrate(v * diff(log_factor, var), var);
                    if (!is_integral_marker(remaining)) {
                        // expand (not simplify) keeps the polynomial terms
                        // distributed: x²·log(x)/2 − x²/4 rather than a
                        // common-factor-wrapped 1/8·(…).
                        return expand(log_factor * v - remaining);
                    }
                }
            }
        }
    }

    // Polynomial × f(affine) for a whitelisted by-parts function f (inverse
    // trig/hyperbolic, erf, Si, …), by parts with u = f, dv = rest dx:
    //   ∫ rest·f dx = f·∫rest − ∫(∫rest)·f' dx.
    // For atan/acot/atanh, f' is rational so the remaining integral is rational
    // (closed by try_rational); for asin/acos/asinh/acosh it is a polynomial over
    // √(quadratic), closed when low-degree. The marker guard bails otherwise.
    if (expr->type_id() == TypeId::Mul) {
        Expr fn_factor;
        std::vector<Expr> rest_factors;
        for (const auto& f : expr->args()) {
            if (!fn_factor && f->type_id() == TypeId::Function) {
                const auto& fn = static_cast<const Function&>(*f);
                if (is_by_parts_fn(fn.function_id()) && fn.args().size() == 1
                    && as_affine(fn.args()[0], var)) {
                    fn_factor = f;
                    continue;
                }
            }
            rest_factors.push_back(f);
        }
        if (fn_factor && !rest_factors.empty()) {
            Expr rest = mul(rest_factors);
            if (is_polynomial_in(rest, var)) {
                Expr v = integrate(rest, var);
                if (!is_integral_marker(v)) {
                    Expr remaining = integrate(v * diff(fn_factor, var), var);
                    if (!is_integral_marker(remaining)) {
                        return expand(fn_factor * v - remaining);
                    }
                }
            }
        }
    }

    // Mul with one factor being sin/cos/exp/sinh/cosh(affine) — or a positive
    // integer power of sin/cos/sinh/cosh(affine), whose antiderivative integrate
    // already supplies — and the rest forming a non-trivial polynomial u. The
    // polynomial is differentiated down each step, so the recursion terminates
    // (the depth guard backs it up); sinh/cosh don't cycle the way exp·sin/cos
    // does. This closes ∫x·cos²(x), ∫x·sin³(x), ∫x·cosh²(x), etc.
    auto is_byparts_target = [&](const Expr& f) -> bool {
        const Function* fn = nullptr;
        if (f->type_id() == TypeId::Function) {
            fn = &static_cast<const Function&>(*f);
        } else if (f->type_id() == TypeId::Pow
                   && f->args()[1]->type_id() == TypeId::Integer
                   && static_cast<const Integer&>(*f->args()[1]).is_positive()
                   && f->args()[0]->type_id() == TypeId::Function) {
            fn = &static_cast<const Function&>(*f->args()[0]);
        }
        if (!fn || fn->args().size() != 1 || !depends_on(f, var)) return false;
        auto aff = as_affine(fn->args()[0], var);
        if (!aff || aff->first == S::Zero()) return false;
        const FunctionId id = fn->function_id();
        // A bare exp is fine; powers are only safe for the non-cyclic trig/hyp
        // functions whose powers integrate (sin/cos/sinh/cosh).
        if (f->type_id() == TypeId::Function && id == FunctionId::Exp) return true;
        return id == FunctionId::Sin || id == FunctionId::Cos
               || id == FunctionId::Sinh || id == FunctionId::Cosh;
    };

    if (expr->type_id() != TypeId::Mul) return std::nullopt;
    Expr target;
    std::vector<Expr> rest_factors;
    for (const auto& f : expr->args()) {
        if (!target && is_byparts_target(f)) {
            target = f;
            continue;
        }
        rest_factors.push_back(f);
    }
    if (!target) return std::nullopt;

    Expr u = mul(rest_factors);
    if (!depends_on(u, var)) return std::nullopt;
    // u must be a polynomial in var so that differentiating it (du each step)
    // terminates. A non-polynomial u like x^(-1) makes by parts recurse on
    // ∫v·u' with an ever-growing negative power (∫exp(x)/x → ∫exp(x)/x² → …),
    // which is non-elementary — bail to the marker instead of looping.
    if (!is_polynomial_in(u, var)) return std::nullopt;

    Expr v = integrate(target, var);
    if (is_integral_marker(v)) return std::nullopt;

    Expr du = diff(u, var);
    // expand so a product like (x/2 + sin(2x)/4)·2x distributes into a sum the
    // linearity path can integrate term by term (else it stays a Mul·Add that no
    // single strategy matches — e.g. ∫x²·cos²(x)).
    Expr remaining = integrate(expand(v * du), var);
    if (is_integral_marker(remaining)) return std::nullopt;

    return u * v - remaining;
}

std::optional<Expr> try_trig_reduction(const Expr& expr, const Expr& var) {
    // sin²(u) → (1 - cos(2u))/2 ; cos²(u) → (1 + cos(2u))/2
    if (expr->type_id() == TypeId::Pow
        && expr->args()[1] == integer(2)) {
        const Expr& base = expr->args()[0];
        if (base->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*base);
            if (fn.args().size() == 1 && depends_on(fn.args()[0], var)) {
                const Expr& u = fn.args()[0];
                if (fn.function_id() == FunctionId::Sin) {
                    Expr rewritten = (integer(1) - cos(integer(2) * u))
                                     / integer(2);
                    return integrate(rewritten, var);
                }
                if (fn.function_id() == FunctionId::Cos) {
                    Expr rewritten = (integer(1) + cos(integer(2) * u))
                                     / integer(2);
                    return integrate(rewritten, var);
                }
                // tan²(u) → 1/cos²(u) − 1 (Pythagorean). Only for an affine u,
                // where ∫1/cos²(u) is tabulated; otherwise fall through.
                if (fn.function_id() == FunctionId::Tan
                    && as_affine(u, var)) {
                    Expr rewritten = pow(cos(u), integer(-2)) - integer(1);
                    return integrate(rewritten, var);
                }
                // The reciprocal trio, rewritten to the 1/cos², 1/sin² table
                // cases (affine u only): sec²(u) = 1/cos²(u),
                // csc²(u) = 1/sin²(u), cot²(u) = 1/sin²(u) − 1 (Pythagorean).
                if (fn.function_id() == FunctionId::Sec
                    && as_affine(u, var)) {
                    return integrate(pow(cos(u), integer(-2)), var);
                }
                if (fn.function_id() == FunctionId::Csc
                    && as_affine(u, var)) {
                    return integrate(pow(sin(u), integer(-2)), var);
                }
                if (fn.function_id() == FunctionId::Cot
                    && as_affine(u, var)) {
                    Expr rewritten = pow(sin(u), integer(-2)) - integer(1);
                    return integrate(rewritten, var);
                }
                // Hyperbolic reciprocal squares, rewritten to the 1/cosh²,
                // 1/sinh² table cases (affine u only): sech²(u) = 1/cosh²(u),
                // csch²(u) = 1/sinh²(u), coth²(u) = 1/sinh²(u) + 1 (since
                // coth² − csch² = 1).
                if (fn.function_id() == FunctionId::Sech
                    && as_affine(u, var)) {
                    return integrate(pow(cosh(u), integer(-2)), var);
                }
                if (fn.function_id() == FunctionId::Csch
                    && as_affine(u, var)) {
                    return integrate(pow(sinh(u), integer(-2)), var);
                }
                if (fn.function_id() == FunctionId::Coth
                    && as_affine(u, var)) {
                    Expr rewritten = pow(sinh(u), integer(-2)) + integer(1);
                    return integrate(rewritten, var);
                }
            }
        }
    }

    // sin(p)cos(q) → (sin(p+q) + sin(p-q))/2 — product-to-sum.
    if (expr->type_id() == TypeId::Mul) {
        Expr sin_arg, cos_arg;
        std::vector<Expr> rest;
        for (const auto& f : expr->args()) {
            if (!sin_arg && f->type_id() == TypeId::Function) {
                const auto& fn = static_cast<const Function&>(*f);
                if (fn.function_id() == FunctionId::Sin
                    && fn.args().size() == 1) {
                    sin_arg = f->args()[0];
                    continue;
                }
            }
            if (!cos_arg && f->type_id() == TypeId::Function) {
                const auto& fn = static_cast<const Function&>(*f);
                if (fn.function_id() == FunctionId::Cos
                    && fn.args().size() == 1) {
                    cos_arg = f->args()[0];
                    continue;
                }
            }
            rest.push_back(f);
        }
        if (sin_arg && cos_arg
            && depends_on(sin_arg, var) && depends_on(cos_arg, var)) {
            Expr rewritten = (sin(sin_arg + cos_arg)
                              + sin(sin_arg - cos_arg))
                             / integer(2);
            for (const auto& r : rest) rewritten = rewritten * r;
            return integrate(rewritten, var);
        }
    }
    return std::nullopt;
}

namespace {

// Parse `expr` as a product of sin(g)^m · cos(g)^n with a single shared
// argument g; returns (g, m, n) with m, n ≥ 0 integers. nullopt if any factor
// is not a non-negative integer power of sin/cos of the same argument.
struct SinCosPowers { Expr g; long m; long n; };

[[nodiscard]] std::optional<SinCosPowers> parse_sin_cos_powers(const Expr& expr) {
    Expr g;
    long m = 0;
    long n = 0;
    auto take = [&](const Expr& f) -> bool {
        Expr base = f;
        long k = 1;
        if (f->type_id() == TypeId::Pow) {
            const Expr& e = f->args()[1];
            if (e->type_id() != TypeId::Integer) return false;
            const auto& z = static_cast<const Integer&>(*e);
            if (!z.is_positive() || !z.fits_long()) return false;
            k = z.to_long();
            base = f->args()[0];
        }
        if (base->type_id() != TypeId::Function) return false;
        const auto& fn = static_cast<const Function&>(*base);
        if (fn.args().size() != 1) return false;
        const FunctionId id = fn.function_id();
        if (id != FunctionId::Sin && id != FunctionId::Cos) return false;
        if (g && !(fn.args()[0] == g)) return false;  // mismatched arguments
        g = fn.args()[0];
        if (id == FunctionId::Sin) m += k; else n += k;
        return true;
    };
    if (expr->type_id() == TypeId::Mul) {
        for (const auto& f : expr->args()) {
            if (!take(f)) return std::nullopt;
        }
    } else if (!take(expr)) {
        return std::nullopt;
    }
    if (!g || (m + n) == 0 || (m + n) > 24) return std::nullopt;
    return SinCosPowers{g, m, n};
}

}  // namespace

std::optional<Expr> try_trig_power(const Expr& expr, const Expr& var) {
    auto parsed = parse_sin_cos_powers(expr);
    if (!parsed) return std::nullopt;
    const Expr& g = parsed->g;
    const long m = parsed->m;
    const long n = parsed->n;
    auto aff = as_affine(g, var);
    if (!aff || aff->first == S::Zero()) return std::nullopt;
    const Expr& a = aff->first;  // dg/dx

    // Odd cos power: u = sin(g), cosⁿ = cos·(1−u²)^((n−1)/2), cos·dx = du/a.
    //   ∫ sinᵐcosⁿ dx = (1/a) ∫ uᵐ (1−u²)^k du.
    // Odd sin power (symmetric): u = cos(g), gives −(1/a)∫ (1−u²)^k uⁿ du.
    Expr u = symbol("_u_trigpow");
    if (n % 2 == 1) {
        long k = (n - 1) / 2;
        Expr poly = expand(pow(u, integer(m))
                           * pow(integer(1) - pow(u, integer(2)), integer(k)));
        Expr anti = integrate(poly, u);
        if (is_integral_marker(anti)) return std::nullopt;
        return expand(subs(anti, u, sin(g)) / a);
    }
    if (m % 2 == 1) {
        long k = (m - 1) / 2;
        Expr poly = expand(pow(integer(1) - pow(u, integer(2)), integer(k))
                           * pow(u, integer(n)));
        Expr anti = integrate(poly, u);
        if (is_integral_marker(anti)) return std::nullopt;
        return expand(mul(S::NegativeOne(), subs(anti, u, cos(g))) / a);
    }

    // Both even: half-angle reduction, then recurse via integrate(). sin²(g) =
    // (1−cos2g)/2, cos²(g) = (1+cos2g)/2; the result is a polynomial in cos(2g)
    // of strictly lower degree, so the recursion (with the depth guard) ends.
    Expr two_g = mul(integer(2), g);
    Expr s2 = (integer(1) - cos(two_g)) / integer(2);
    Expr c2 = (integer(1) + cos(two_g)) / integer(2);
    Expr rewritten = pow(s2, integer(m / 2)) * pow(c2, integer(n / 2));
    return integrate(expand(rewritten), var);
}

std::optional<Expr> try_tan_power(const Expr& expr, const Expr& var) {
    if (expr->type_id() != TypeId::Pow) return std::nullopt;
    const Expr& base = expr->args()[0];
    const Expr& e = expr->args()[1];
    if (e->type_id() != TypeId::Integer) return std::nullopt;
    const auto& z = static_cast<const Integer&>(*e);
    if (!z.fits_long()) return std::nullopt;
    const long n = z.to_long();
    if (n < 2 || n > 24) return std::nullopt;  // n=1 is the table case
    if (base->type_id() != TypeId::Function) return std::nullopt;
    const auto& fn = static_cast<const Function&>(*base);
    const FunctionId id = fn.function_id();
    if ((id != FunctionId::Tan && id != FunctionId::Cot)
        || fn.args().size() != 1) {
        return std::nullopt;
    }
    const Expr& g = fn.args()[0];
    auto aff = as_affine(g, var);
    if (!aff || aff->first == S::Zero()) return std::nullopt;
    const Expr& a = aff->first;

    // Power-reduction by the Pythagorean identity, recursing down to the n=1
    // table case (∫tan = −log(cos)/a, ∫cot = log(sin)/a):
    //   ∫tanⁿ =  tan^(n-1)/((n-1)·g') − ∫tan^(n-2)
    //   ∫cotⁿ = −cot^(n-1)/((n-1)·g') − ∫cot^(n-2)
    if (id == FunctionId::Tan) {
        Expr first = pow(tan(g), integer(n - 1)) / (integer(n - 1) * a);
        Expr rest = integrate(pow(tan(g), integer(n - 2)), var);
        if (is_integral_marker(rest)) return std::nullopt;
        return first - rest;
    }
    Expr first = mul(S::NegativeOne(), pow(cot(g), integer(n - 1)))
                 / (integer(n - 1) * a);
    Expr rest = integrate(pow(cot(g), integer(n - 2)), var);
    if (is_integral_marker(rest)) return std::nullopt;
    return first - rest;
}

std::optional<Expr> try_tanh_coth_power(const Expr& expr, const Expr& var) {
    if (expr->type_id() != TypeId::Pow) return std::nullopt;
    const Expr& base = expr->args()[0];
    const Expr& e = expr->args()[1];
    if (e->type_id() != TypeId::Integer) return std::nullopt;
    const auto& z = static_cast<const Integer&>(*e);
    if (!z.fits_long()) return std::nullopt;
    const long n = z.to_long();
    if (n < 2 || n > 24) return std::nullopt;  // n=1 is the table case
    if (base->type_id() != TypeId::Function) return std::nullopt;
    const auto& fn = static_cast<const Function&>(*base);
    const FunctionId id = fn.function_id();
    if ((id != FunctionId::Tanh && id != FunctionId::Coth)
        || fn.args().size() != 1) {
        return std::nullopt;
    }
    const Expr& g = fn.args()[0];
    auto aff = as_affine(g, var);
    if (!aff || aff->first == S::Zero()) return std::nullopt;
    const Expr& a = aff->first;

    // Both tanh and coth share the same reduction, recursing to the n=1 table
    // case (∫tanh = log(cosh)/g', ∫coth = log(sinh)/g') and the n=0 case ∫1 = x:
    //   ∫tanhⁿ = ∫tanh^(n-2) − tanh^(n-1)/((n-1)·g')   (tanh² = 1 − sech²)
    //   ∫cothⁿ = ∫coth^(n-2) − coth^(n-1)/((n-1)·g')   (coth² = 1 + csch²)
    Expr first = pow(base, integer(n - 1)) / (integer(n - 1) * a);
    Expr rest = integrate(pow(base, integer(n - 2)), var);
    if (is_integral_marker(rest)) return std::nullopt;
    return rest - first;
}

std::optional<Expr> try_sec_csc_power(const Expr& expr, const Expr& var) {
    if (expr->type_id() != TypeId::Pow) return std::nullopt;
    const Expr& base = expr->args()[0];
    const Expr& e = expr->args()[1];
    if (e->type_id() != TypeId::Integer) return std::nullopt;
    const auto& z = static_cast<const Integer&>(*e);
    if (!z.fits_long()) return std::nullopt;
    const long n = z.to_long();
    // n=1 is the table case (INT-24); n=2 the trig-reduction square (INT-25).
    if (n < 3 || n > 24) return std::nullopt;
    if (base->type_id() != TypeId::Function) return std::nullopt;
    const auto& fn = static_cast<const Function&>(*base);
    const FunctionId id = fn.function_id();
    if ((id != FunctionId::Sec && id != FunctionId::Csc)
        || fn.args().size() != 1) {
        return std::nullopt;
    }
    const Expr& g = fn.args()[0];
    auto aff = as_affine(g, var);
    if (!aff || aff->first == S::Zero()) return std::nullopt;
    const Expr& a = aff->first;

    // By-parts reduction, recursing to the n=1 table / n=2 square cases:
    //   ∫secⁿ =  sec^(n-2)·tan/((n-1)·g') + (n-2)/(n-1)·∫sec^(n-2)
    //   ∫cscⁿ = −csc^(n-2)·cot/((n-1)·g') + (n-2)/(n-1)·∫csc^(n-2)
    Expr rest = integrate(pow(base, integer(n - 2)), var);
    if (is_integral_marker(rest)) return std::nullopt;
    Expr rest_term = rational(n - 2, n - 1) * rest;
    if (id == FunctionId::Sec) {
        Expr boundary = pow(sec(g), integer(n - 2)) * tan(g)
                        / (integer(n - 1) * a);
        return boundary + rest_term;
    }
    Expr boundary = mul(S::NegativeOne(), pow(csc(g), integer(n - 2)) * cot(g))
                    / (integer(n - 1) * a);
    return boundary + rest_term;
}

std::optional<Expr> try_sech_csch_power(const Expr& expr, const Expr& var) {
    if (expr->type_id() != TypeId::Pow) return std::nullopt;
    const Expr& base = expr->args()[0];
    const Expr& e = expr->args()[1];
    if (e->type_id() != TypeId::Integer) return std::nullopt;
    const auto& z = static_cast<const Integer&>(*e);
    if (!z.fits_long()) return std::nullopt;
    const long n = z.to_long();
    // n=1 / n=2 are the table / square cases (INT-26).
    if (n < 3 || n > 24) return std::nullopt;
    if (base->type_id() != TypeId::Function) return std::nullopt;
    const auto& fn = static_cast<const Function&>(*base);
    const FunctionId id = fn.function_id();
    if ((id != FunctionId::Sech && id != FunctionId::Csch)
        || fn.args().size() != 1) {
        return std::nullopt;
    }
    const Expr& g = fn.args()[0];
    auto aff = as_affine(g, var);
    if (!aff || aff->first == S::Zero()) return std::nullopt;
    const Expr& a = aff->first;

    // By-parts reduction, recursing to the n=1 table / n=2 square cases. Unlike
    // the trig case, csch subtracts the rest term (coth² = 1 + csch²):
    //   ∫sechⁿ =  sech^(n-2)·tanh/((n-1)·g') + (n-2)/(n-1)·∫sech^(n-2)
    //   ∫cschⁿ = −csch^(n-2)·coth/((n-1)·g') − (n-2)/(n-1)·∫csch^(n-2)
    Expr rest = integrate(pow(base, integer(n - 2)), var);
    if (is_integral_marker(rest)) return std::nullopt;
    Expr rest_mag = rational(n - 2, n - 1) * rest;
    if (id == FunctionId::Sech) {
        Expr boundary = pow(sech(g), integer(n - 2)) * tanh(g)
                        / (integer(n - 1) * a);
        return boundary + rest_mag;
    }
    Expr boundary = mul(S::NegativeOne(), pow(csch(g), integer(n - 2)) * coth(g))
                    / (integer(n - 1) * a);
    return boundary - rest_mag;
}

namespace {

// sinh(g)^m · cosh(g)^n parser, the hyperbolic analogue of
// parse_sin_cos_powers.
[[nodiscard]] std::optional<SinCosPowers> parse_sinh_cosh_powers(
    const Expr& expr) {
    Expr g;
    long m = 0;
    long n = 0;
    auto take = [&](const Expr& f) -> bool {
        Expr base = f;
        long k = 1;
        if (f->type_id() == TypeId::Pow) {
            const Expr& e = f->args()[1];
            if (e->type_id() != TypeId::Integer) return false;
            const auto& z = static_cast<const Integer&>(*e);
            if (!z.is_positive() || !z.fits_long()) return false;
            k = z.to_long();
            base = f->args()[0];
        }
        if (base->type_id() != TypeId::Function) return false;
        const auto& fn = static_cast<const Function&>(*base);
        if (fn.args().size() != 1) return false;
        const FunctionId id = fn.function_id();
        if (id != FunctionId::Sinh && id != FunctionId::Cosh) return false;
        if (g && !(fn.args()[0] == g)) return false;
        g = fn.args()[0];
        if (id == FunctionId::Sinh) m += k; else n += k;
        return true;
    };
    if (expr->type_id() == TypeId::Mul) {
        for (const auto& f : expr->args()) {
            if (!take(f)) return std::nullopt;
        }
    } else if (!take(expr)) {
        return std::nullopt;
    }
    if (!g || (m + n) == 0 || (m + n) > 24) return std::nullopt;
    return SinCosPowers{g, m, n};
}

}  // namespace

std::optional<Expr> try_hyperbolic_power(const Expr& expr, const Expr& var) {
    auto parsed = parse_sinh_cosh_powers(expr);
    if (!parsed) return std::nullopt;
    const Expr& g = parsed->g;
    const long m = parsed->m;
    const long n = parsed->n;
    auto aff = as_affine(g, var);
    if (!aff || aff->first == S::Zero()) return std::nullopt;
    const Expr& a = aff->first;

    Expr u = symbol("_u_hyppow");
    // Odd cosh power: u = sinh(g), cosh^(n-1) = (1+u²)^k, cosh·dx = du/a.
    if (n % 2 == 1) {
        long k = (n - 1) / 2;
        Expr poly = expand(pow(u, integer(m))
                           * pow(integer(1) + pow(u, integer(2)), integer(k)));
        Expr anti = integrate(poly, u);
        if (is_integral_marker(anti)) return std::nullopt;
        return expand(subs(anti, u, sinh(g)) / a);
    }
    // Odd sinh power: u = cosh(g), sinh^(m-1) = (u²−1)^k, sinh·dx = du/a.
    if (m % 2 == 1) {
        long k = (m - 1) / 2;
        Expr poly = expand(pow(pow(u, integer(2)) - integer(1), integer(k))
                           * pow(u, integer(n)));
        Expr anti = integrate(poly, u);
        if (is_integral_marker(anti)) return std::nullopt;
        return expand(subs(anti, u, cosh(g)) / a);
    }
    // Both even: sinh²(g) = (cosh2g−1)/2, cosh²(g) = (cosh2g+1)/2; recurse.
    Expr two_g = mul(integer(2), g);
    Expr s2 = (cosh(two_g) - integer(1)) / integer(2);
    Expr c2 = (cosh(two_g) + integer(1)) / integer(2);
    Expr rewritten = pow(s2, integer(m / 2)) * pow(c2, integer(n / 2));
    return integrate(expand(rewritten), var);
}

std::optional<Expr> try_rational(const Expr& expr, const Expr& var) {
    // Bring to single fraction.
    Expr t = together(expr);
    Expr num = S::One();
    Expr den;
    if (t->type_id() == TypeId::Pow
        && t->args()[1] == S::NegativeOne()) {
        den = t->args()[0];
    } else if (t->type_id() == TypeId::Mul) {
        std::vector<Expr> num_factors;
        for (const auto& f : t->args()) {
            if (f->type_id() == TypeId::Pow
                && f->args()[1] == S::NegativeOne()) {
                if (den) return std::nullopt;
                den = f->args()[0];
            } else {
                num_factors.push_back(f);
            }
        }
        if (!den) return std::nullopt;
        num = mul(num_factors);
    } else {
        return std::nullopt;
    }

    // Verify num and den are pure polynomials in var (constant
    // coefficients only — reject things like sin(x)/cos(x)).
    Poly num_p(expand(num), var);
    Poly den_p(expand(den), var);
    for (const auto& c : num_p.coeffs()) if (has(c, var)) return std::nullopt;
    for (const auto& c : den_p.coeffs()) if (has(c, var)) return std::nullopt;
    if (den_p.degree() == 0) return std::nullopt;

    // Polynomial division: num = q * den + r, deg(r) < deg(den).
    auto [q, r_poly] = num_p.divmod(den_p);

    if (r_poly.is_zero()) {
        // Pure polynomial result — let the linearity path handle it.
        return integrate(q.as_expr(), var);
    }

    // Decompose the proper remainder via apart().
    Expr proper = r_poly.as_expr() / den_p.as_expr();
    Expr apart_form = apart(proper, var);

    Expr proper_int;
    if (apart_form == proper) {
        // apart couldn't split the proper part — it is a single irreducible
        // quadratic denominator with a constant or linear numerator. For a bare
        // proper fraction (q == 0) defer to the dedicated quadratic helpers
        // downstream (re-integrating here would recurse). For the improper case
        // (q ≠ 0) the polynomial quotient must still be integrated, so close the
        // remainder directly via those helpers instead of dropping the whole
        // result to the marker.
        if (q.is_zero()) return std::nullopt;
        if (den_p.degree() != 2) return std::nullopt;
        if (r_poly.degree() <= 0) {
            auto inv = try_arctan_quadratic(
                pow(den_p.as_expr(), S::NegativeOne()), var);
            if (!inv) return std::nullopt;
            proper_int = mul(r_poly.as_expr(), inv.value());
        } else {
            auto lin = try_linear_over_quadratic(proper, var);
            if (!lin) return std::nullopt;
            proper_int = lin.value();
        }
    } else {
        proper_int = integrate(apart_form, var);
        if (is_integral_marker(proper_int)) return std::nullopt;
    }

    return integrate(q.as_expr(), var) + proper_int;
}

std::optional<Expr> try_hyperbolic_to_exp(const Expr& expr, const Expr& var) {
    if (expr->type_id() != TypeId::Mul) return std::nullopt;

    // Require at least one sinh/cosh(affine) factor and at least one
    // sin/cos/exp(affine) factor — so rewriting the hyperbolics to exponentials
    // pairs them with an exp·trig cyclic form or collapses to a pure exponential.
    // (Pure sinh·cosh products are try_hyperbolic_power's job.)
    bool has_hyp = false, has_partner = false;
    for (const auto& f : expr->args()) {
        if (f->type_id() != TypeId::Function) continue;
        const auto& fn = static_cast<const Function&>(*f);
        if (fn.args().size() != 1 || !as_affine(fn.args()[0], var)) continue;
        const FunctionId id = fn.function_id();
        if (id == FunctionId::Sinh || id == FunctionId::Cosh) has_hyp = true;
        else if (id == FunctionId::Sin || id == FunctionId::Cos
                 || id == FunctionId::Exp) {
            has_partner = true;
        }
    }
    if (!has_hyp || !has_partner) return std::nullopt;

    // Rewrite every sinh/cosh(affine) factor to exponentials.
    std::vector<Expr> factors;
    for (const auto& f : expr->args()) {
        if (f->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*f);
            if (fn.args().size() == 1) {
                const Expr& g = fn.args()[0];
                if (fn.function_id() == FunctionId::Sinh) {
                    factors.push_back((exp(g) - exp(mul(S::NegativeOne(), g)))
                                      / integer(2));
                    continue;
                }
                if (fn.function_id() == FunctionId::Cosh) {
                    factors.push_back((exp(g) + exp(mul(S::NegativeOne(), g)))
                                      / integer(2));
                    continue;
                }
            }
        }
        factors.push_back(f);
    }

    // Merge the exponential factors within a single product term: the canonical
    // Mul does not fold e^a·e^b → e^(a+b), so e^x·e^−x would otherwise stay an
    // unevaluated product and block the pure-exponential (exp × cosh/sinh) case.
    auto combine_exps = [](const Expr& term) -> Expr {
        if (term->type_id() != TypeId::Mul) return term;
        Expr exp_arg = S::Zero();
        std::vector<Expr> others;
        bool found = false;
        for (const auto& f : term->args()) {
            if (f->type_id() == TypeId::Function) {
                const auto& fn = static_cast<const Function&>(*f);
                if (fn.function_id() == FunctionId::Exp && fn.args().size() == 1) {
                    exp_arg = exp_arg + fn.args()[0];
                    found = true;
                    continue;
                }
            }
            others.push_back(f);
        }
        if (!found) return term;
        others.push_back(exp(exp_arg));  // exp(0) → 1 folds in the factory
        return mul(others);
    };

    Expr rewritten = expand(mul(factors));
    if (rewritten->type_id() == TypeId::Add) {
        std::vector<Expr> terms;
        terms.reserve(rewritten->args().size());
        for (const auto& term : rewritten->args()) {
            terms.push_back(combine_exps(term));
        }
        rewritten = add(terms);
    } else {
        rewritten = combine_exps(rewritten);
    }

    Expr result = integrate(rewritten, var);
    if (is_integral_marker(result)) return std::nullopt;
    return result;
}

// True if `e` is a rational function of `t`: built only from t, t-free
// constants, +, ×, and integer powers. A non-integer power or a function of t
// makes it non-rational. Used to gate the Weierstrass result so it only claims
// integrands that become genuinely rational in t (otherwise the t-integral can
// be non-elementary and send `integrate` into a long/unbounded search).
[[nodiscard]] bool is_rational_in(const Expr& e, const Expr& t) {
    if (!depends_on(e, t)) return true;
    switch (e->type_id()) {
        case TypeId::Symbol:
            return true;  // depends on t ⇒ is t
        case TypeId::Add:
        case TypeId::Mul:
            for (const auto& a : e->args()) {
                if (!is_rational_in(a, t)) return false;
            }
            return true;
        case TypeId::Pow:
            return e->args()[1]->type_id() == TypeId::Integer
                   && is_rational_in(e->args()[0], t);
        default:
            return false;  // a function of t, fractional power, etc.
    }
}

// True if `e` contains a trig function of var raised to a power (i.e. a Pow node
// whose base is sin/cos/tan/cot/sec/csc(…var…)). The Weierstrass substitution
// turns such powers into high-degree nested rationals in t whose normalisation
// (cancel) or integration (try_rational's Poly GCD) can run away, so they are
// excluded; trig appearing only to the first power inside a polynomial
// denominator (the classic ∫1/(a+b·cos x) family) is fine.
[[nodiscard]] bool has_trig_power_of(const Expr& e, const Expr& var) {
    if (e->type_id() == TypeId::Pow) {
        const Expr& base = e->args()[0];
        if (base->type_id() == TypeId::Function && depends_on(base, var)) {
            switch (static_cast<const Function&>(*base).function_id()) {
                case FunctionId::Sin:
                case FunctionId::Cos:
                case FunctionId::Tan:
                case FunctionId::Cot:
                case FunctionId::Sec:
                case FunctionId::Csc:
                    return true;
                default:
                    break;
            }
        }
    }
    for (const auto& a : e->args()) {
        if (has_trig_power_of(a, var)) return true;
    }
    return false;
}

std::optional<Expr> try_weierstrass(const Expr& expr, const Expr& var) {
    if (!depends_on(expr, var)) return std::nullopt;
    // Exclude integrands with a trig function raised to a power — their
    // substituted form can send cancel()/try_rational into a runaway. The common
    // rational-trig integrals (denominator linear in sin/cos) have no such power.
    if (has_trig_power_of(expr, var)) return std::nullopt;
    Expr t = symbol("_weierstrass_t");

    // Rewrite the reciprocal/quotient trig functions of var into sin/cos(var),
    // then apply the half-angle substitution. Each subs is a no-op when that
    // function is absent.
    Expr e = expr;
    e = subs(e, tan(var), sin(var) / cos(var));
    e = subs(e, cot(var), cos(var) / sin(var));
    e = subs(e, sec(var), pow(cos(var), S::NegativeOne()));
    e = subs(e, csc(var), pow(sin(var), S::NegativeOne()));

    // sin(var) = 2t/(1+t²), cos(var) = (1−t²)/(1+t²).
    Expr one_pt2 = integer(1) + pow(t, integer(2));
    Expr sin_sub = integer(2) * t / one_pt2;
    Expr cos_sub = (integer(1) - pow(t, integer(2))) / one_pt2;
    e = subs(e, sin(var), sin_sub);
    e = subs(e, cos(var), cos_sub);

    // Any surviving var means the integrand was not a rational function of the
    // trig functions of the bare argument var (e.g. a polynomial factor, an
    // affine trig argument like sin(2x), or exp/log of var) — bail.
    if (depends_on(e, var)) return std::nullopt;

    // dx = 2/(1+t²) dt; bring to a single fraction.
    Expr integrand = together(e * integer(2) / one_pt2);

    // A non-rational integrand (e.g. √(tan x) → √(2t/(1−t²))) would hand
    // `integrate` a non-elementary algebraic integral that can loop — bail.
    if (!is_rational_in(integrand, t)) return std::nullopt;

    Expr antideriv = integrate(integrand, t);
    if (is_integral_marker(antideriv)) return std::nullopt;

    // Back-substitute t = tan(var/2).
    return simplify(subs(antideriv, t, tan(var / integer(2))));
}

std::optional<Expr> try_quadratic_power(const Expr& expr, const Expr& var) {
    if (expr->type_id() != TypeId::Pow) return std::nullopt;
    const Expr& base = expr->args()[0];
    const Expr& e = expr->args()[1];
    if (e->type_id() != TypeId::Integer) return std::nullopt;
    const auto& z = static_cast<const Integer&>(*e);
    if (!z.fits_long()) return std::nullopt;
    const long ei = z.to_long();
    if (ei >= -1) return std::nullopt;       // n=1 (exp −1) is try_arctan's job
    const long n = -ei;
    if (n > 24) return std::nullopt;
    if (!depends_on(base, var)) return std::nullopt;

    Poly p(expand(base), var);
    if (p.degree() != 2) return std::nullopt;
    const Expr& c = p.coeffs()[0];
    const Expr& b = p.coeffs()[1];
    const Expr& a = p.coeffs()[2];
    if (is_rational(a) != true || is_rational(b) != true
        || is_rational(c) != true) {
        return std::nullopt;
    }

    // Complete the square so the linear term vanishes: u = x + b/(2a) gives
    // a·u² + (c − b²/(4a)); du = dx, so back-substitute afterwards.
    if (!(b == S::Zero())) {
        Expr shift = b / (integer(2) * a);
        Expr cprime = simplify(c - b * b / (integer(4) * a));
        Expr shifted = pow(a * pow(var, integer(2)) + cprime, e);
        auto inner = try_quadratic_power(shifted, var);
        if (!inner.has_value()) return std::nullopt;
        return simplify(subs(inner.value(), var, var + shift));
    }
    if (c == S::Zero()) return std::nullopt;  // a·x² alone: not this rule's job

    // Iₙ = x/(2(n−1)c·Qⁿ⁻¹) + (2n−3)/(2(n−1)c)·Iₙ₋₁, Q = a·x²+c. The leading
    // coefficient a cancels in the derivation (x² = (Q−c)/a, 2a·x² = 2(Q−c)), so
    // it does not appear here. Recurse through integrate: Iₙ₋₁ routes back here
    // for n−1 ≥ 2, else to try_arctan_quadratic / try_rational at the
    // I₁ = ∫1/(a·x²+c) base case.
    Expr prev = integrate(pow(base, integer(-(n - 1))), var);
    if (is_integral_marker(prev)) return std::nullopt;
    Expr denom = integer(2) * integer(n - 1) * c;
    Expr term1 = var / (denom * pow(base, integer(n - 1)));
    Expr coeff = integer(2 * n - 3) / denom;
    return simplify(term1 + coeff * prev);
}

std::optional<Expr> try_arctan_quadratic(const Expr& expr, const Expr& var) {
    // Match 1/(quadratic): a Pow with exponent −1 over a var-dependent base.
    if (expr->type_id() != TypeId::Pow) return std::nullopt;
    if (!(expr->args()[1] == S::NegativeOne())) return std::nullopt;
    const Expr& base = expr->args()[0];
    if (!depends_on(base, var)) return std::nullopt;

    Poly p(expand(base), var);
    if (p.degree() != 2) return std::nullopt;
    const Expr& c = p.coeffs()[0];
    const Expr& b = p.coeffs()[1];
    const Expr& a = p.coeffs()[2];

    // Need rational coefficients to decide irreducibility from the sign of the
    // discriminant. Symbolic coefficients (e.g. 1/(k²+x²)) are out of scope.
    if (is_rational(a) != true || is_rational(b) != true
        || is_rational(c) != true) {
        return std::nullopt;
    }

    // D = 4ac − b². D > 0 ⇒ no real roots ⇒ arctangent. D = 0 ⇒ a repeated
    // real root. D < 0 ⇒ distinct real roots, which try_rational splits into
    // logs, so leave those alone.
    Expr disc = integer(4) * a * c - b * b;

    if (disc == S::Zero()) {
        // a·x² + b·x + c = a·(x − r)², r = −b/(2a):
        //   ∫ 1/(a·(x − r)²) dx = −2/(2a·x + b).
        return integer(-2) / (integer(2) * a * var + b);
    }
    if (is_positive(disc) != true) return std::nullopt;

    // ∫ 1/(a·x² + b·x + c) dx = 2·atan((2a·x + b)/√D) / √D.
    Expr root_d = sqrt(disc);
    Expr arg = (integer(2) * a * var + b) / root_d;
    return simplify(integer(2) * atan(arg) / root_d);
}

std::optional<Expr> try_linear_over_quadratic(const Expr& expr,
                                              const Expr& var) {
    // Match num · den^(-1) where den is a quadratic and num is linear in var.
    Expr num;
    Expr den;
    if (expr->type_id() == TypeId::Mul) {
        std::vector<Expr> num_factors;
        for (const auto& f : expr->args()) {
            if (f->type_id() == TypeId::Pow && f->args()[1] == S::NegativeOne()) {
                if (den) return std::nullopt;  // more than one denominator
                den = f->args()[0];
            } else {
                num_factors.push_back(f);
            }
        }
        if (!den) return std::nullopt;
        num = mul(num_factors);
    } else {
        return std::nullopt;  // a constant numerator is try_arctan_quadratic's job
    }
    if (!depends_on(den, var)) return std::nullopt;

    Poly dp(expand(den), var);
    Poly np(expand(num), var);
    if (dp.degree() != 2 || np.degree() != 1) return std::nullopt;
    const Expr& c = dp.coeffs()[0];
    const Expr& b = dp.coeffs()[1];
    const Expr& a = dp.coeffs()[2];
    const Expr& q = np.coeffs()[0];
    const Expr& p = np.coeffs()[1];
    // Numeric coefficients only (irreducibility is decided from the sign of D).
    for (const auto& e : {a, b, c, p, q}) {
        if (is_rational(e) != true) return std::nullopt;
    }
    // Only the irreducible case (D = 4ac − b² > 0). Reducible / repeated-root
    // denominators are split by try_rational into partial fractions.
    Expr disc = integer(4) * a * c - b * b;
    if (is_positive(disc) != true) return std::nullopt;

    // p·x + q = (p/2a)·(2a·x + b) + (q − p·b/2a):
    //   ∫ = (p/2a)·log(a·x²+b·x+c) + (q − p·b/2a)·∫ 1/(a·x²+b·x+c) dx.
    Expr two_a = integer(2) * a;
    Expr log_coeff = p / two_a;
    Expr remainder = q - p * b / two_a;
    Expr result = mul(log_coeff, log(dp.as_expr()));
    if (!(simplify(remainder) == S::Zero())) {
        auto inv = try_arctan_quadratic(pow(den, S::NegativeOne()), var);
        if (!inv) return std::nullopt;  // shouldn't happen given D > 0
        result = result + remainder * (*inv);
    }
    return simplify(result);
}

std::optional<Expr> try_sqrt_quadratic(const Expr& expr, const Expr& var) {
    // Match 1/√(quadratic) [exponent −1/2] or √(quadratic) [exponent +1/2].
    if (expr->type_id() != TypeId::Pow) return std::nullopt;
    const bool reciprocal = (expr->args()[1] == rational(-1, 2));
    const bool numerator = (expr->args()[1] == rational(1, 2));
    if (!reciprocal && !numerator) return std::nullopt;
    const Expr& base = expr->args()[0];
    if (!depends_on(base, var)) return std::nullopt;

    Poly p(expand(base), var);
    if (p.degree() != 2) return std::nullopt;
    const Expr& c = p.coeffs()[0];
    const Expr& b = p.coeffs()[1];
    const Expr& a = p.coeffs()[2];

    // Completing the square: with a linear term, the substitution u = x + b/(2a)
    // shifts the radicand to a·u² + (c − b²/(4a)) and kills the linear term.
    // Since du = dx, the antiderivative is just the pure-quadratic result with
    // x ← x + b/(2a). Reuse this routine recursively for the shifted (b = 0)
    // radicand, for both the reciprocal and numerator exponents.
    if (!(b == S::Zero())) {
        if (is_rational(a) != true || is_rational(b) != true
            || is_rational(c) != true) {
            return std::nullopt;
        }
        Expr shift = b / (integer(2) * a);
        Expr cprime = simplify(c - b * b / (integer(4) * a));
        Expr shifted = pow(a * pow(var, integer(2)) + cprime, expr->args()[1]);
        auto inner = try_sqrt_quadratic(shifted, var);
        if (!inner.has_value()) return std::nullopt;
        return simplify(subs(inner.value(), var, var + shift));
    }
    if (is_rational(a) != true || is_rational(c) != true) return std::nullopt;

    // ∫ √(a·x² + c) dx = (x/2)·√(a·x²+c) + (c/2)·∫ 1/√(a·x²+c) dx (by parts).
    // Reduce to the reciprocal case below and reuse its asin/asinh/log result;
    // a nullopt inner integral (c == 0, or a < 0 with c ≤ 0 — no real region)
    // propagates, so those fall through unevaluated as before.
    if (numerator) {
        auto inv = try_sqrt_quadratic(pow(base, rational(-1, 2)), var);
        if (!inv.has_value()) return std::nullopt;
        return simplify(var * expr / integer(2) + c * inv.value() / integer(2));
    }

    if (is_positive(a) == true && is_positive(c) == true) {
        // ∫ 1/√(a·x² + c) dx = asinh(x·√(a/c)) / √a.
        return simplify(asinh(var * sqrt(a / c)) / sqrt(a));
    }
    if (is_negative(a) == true && is_positive(c) == true) {
        // a < 0: ∫ 1/√(c − |a|·x²) dx = asin(x·√(|a|/c)) / √|a|.
        Expr a_pos = mul(S::NegativeOne(), a);
        return simplify(asin(var * sqrt(a_pos / c)) / sqrt(a_pos));
    }
    if (is_positive(a) == true && is_negative(c) == true) {
        // a > 0, c < 0: ∫ 1/√(a·x² + c) dx = log(√a·x + √(a·x²+c)) / √a.
        // (acosh-equivalent; this is the form SymPy prints.) The radicand
        // √(a·x²+c) is real only where a·x²+c ≥ 0, as usual for this integral.
        Expr sa = sqrt(a);
        return simplify(log(sa * var + pow(base, rational(1, 2))) / sa);
    }
    // a < 0, c ≤ 0 (radicand has no positive region) and c == 0 fall through.
    return std::nullopt;
}

std::optional<Expr> try_x_over_sqrt_quadratic(const Expr& expr, const Expr& var) {
    if (expr->type_id() != TypeId::Mul) return std::nullopt;
    // Split into the √-reciprocal factor (quadratic)^(-1/2) and a numerator that
    // collects everything else (a linear polynomial in var, p·x + q).
    Expr radical;                 // pow(quad, -1/2)
    std::vector<Expr> num_factors;
    for (const auto& f : expr->args()) {
        if (!radical && f->type_id() == TypeId::Pow
            && f->args()[1] == rational(-1, 2) && depends_on(f->args()[0], var)) {
            radical = f;
            continue;
        }
        num_factors.push_back(f);
    }
    if (!radical || num_factors.empty()) return std::nullopt;

    Poly q_poly(expand(radical->args()[0]), var);
    if (q_poly.degree() != 2) return std::nullopt;
    const Expr& c = q_poly.coeffs()[0];
    const Expr& b = q_poly.coeffs()[1];
    const Expr& a = q_poly.coeffs()[2];
    if (is_rational(a) != true || is_rational(b) != true
        || is_rational(c) != true) {
        return std::nullopt;
    }

    // The numerator must be linear in var: N = p·x + q.
    Poly n_poly(expand(mul(num_factors)), var);
    if (n_poly.degree() > 1) return std::nullopt;
    const Expr p = n_poly.degree() >= 1 ? n_poly.coeffs()[1] : S::Zero();
    const Expr& qn = n_poly.coeffs()[0];

    // d/dx √(a·x²+b·x+c) = (2a·x+b)/(2√Q), so
    //   ∫ x/√Q dx = √Q/a − (b/2a)·∫ 1/√Q dx,
    // and for the general linear numerator N = p·x + q:
    //   ∫ N/√Q dx = (p/a)·√Q + (q − p·b/(2a))·∫ 1/√Q dx.
    // The reciprocal term (present whenever b ≠ 0 or q ≠ p·b/(2a)) is handled by
    // try_sqrt_quadratic, including its completing-the-square branch.
    Expr sqrtQ = pow(radical->args()[0], rational(1, 2));
    Expr result = p / a * sqrtQ;
    Expr recip_coeff = simplify(qn - p * b / (integer(2) * a));
    if (!(recip_coeff == S::Zero())) {
        auto recip = try_sqrt_quadratic(pow(radical->args()[0], rational(-1, 2)),
                                        var);
        if (!recip.has_value()) return std::nullopt;
        result = result + recip_coeff * recip.value();
    }
    return simplify(result);
}

std::optional<Expr> try_gaussian(const Expr& expr, const Expr& var) {
    if (expr->type_id() != TypeId::Function) return std::nullopt;
    const auto& fn = static_cast<const Function&>(*expr);
    if (fn.function_id() != FunctionId::Exp || fn.args().size() != 1) {
        return std::nullopt;
    }
    Poly p(expand(fn.args()[0]), var);
    if (p.degree() != 2) return std::nullopt;
    const Expr& c0 = p.coeffs()[0];
    const Expr& c1 = p.coeffs()[1];
    const Expr& c2 = p.coeffs()[2];
    // Pure c·x² only (no linear or constant term — a linear/constant part needs
    // completing the square, out of scope), with c a non-zero rational.
    // Negative c → real erf; positive c → the imaginary error function erfi.
    if (!(c0 == S::Zero()) || !(c1 == S::Zero())) return std::nullopt;
    if (is_rational(c2) != true) return std::nullopt;

    if (is_negative(c2) == true) {
        // ∫ exp(−a·x²) dx = √π·erf(√a·x) / (2·√a),  a = −c₂ > 0.
        Expr a = mul(S::NegativeOne(), c2);
        Expr sa = sqrt(a);
        return simplify(sqrt(S::Pi()) * erf(sa * var) / (integer(2) * sa));
    }
    if (is_positive(c2) == true) {
        // ∫ exp(a·x²) dx = √π·erfi(√a·x) / (2·√a),  a = c₂ > 0.
        Expr sa = sqrt(c2);
        return simplify(sqrt(S::Pi()) * erfi(sa * var) / (integer(2) * sa));
    }
    return std::nullopt;
}

std::optional<Expr> try_expint_integral(const Expr& expr, const Expr& var) {
    if (expr->type_id() != TypeId::Mul) return std::nullopt;
    Expr func;                    // the sin/cos/exp factor
    bool has_recip = false;       // saw a 1/var factor
    std::vector<Expr> consts;     // constant prefactors
    for (const auto& f : expr->args()) {
        if (!has_recip && f->type_id() == TypeId::Pow
            && f->args()[0] == var && f->args()[1] == S::NegativeOne()) {
            has_recip = true;
            continue;
        }
        if (!func && f->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*f);
            const FunctionId id = fn.function_id();
            if ((id == FunctionId::Sin || id == FunctionId::Cos
                 || id == FunctionId::Exp || id == FunctionId::Sinh
                 || id == FunctionId::Cosh) && fn.args().size() == 1) {
                auto aff = as_affine(fn.args()[0], var);
                // Monomial argument c·x only (no constant term b): ∫f(c·x+b)/x
                // is not an elementary special-integral function.
                if (aff && !(aff->first == S::Zero())
                    && aff->second == S::Zero()) {
                    func = f;
                    continue;
                }
            }
        }
        if (depends_on(f, var)) return std::nullopt;  // some other var factor
        consts.push_back(f);
    }
    if (!has_recip || !func) return std::nullopt;

    const auto& fn = static_cast<const Function&>(*func);
    const Expr& g = fn.args()[0];
    Expr result;
    switch (fn.function_id()) {
        case FunctionId::Sin: result = sinint(g); break;     // ∫sin(c·x)/x = Si
        case FunctionId::Cos: result = cosint(g); break;     // ∫cos(c·x)/x = Ci
        case FunctionId::Exp: result = expint_ei(g); break;  // ∫exp(c·x)/x = Ei
        case FunctionId::Sinh: result = sinhint(g); break;   // ∫sinh(c·x)/x = Shi
        case FunctionId::Cosh: result = coshint(g); break;   // ∫cosh(c·x)/x = Chi
        default: return std::nullopt;
    }
    if (!consts.empty()) result = mul(mul(consts), result);
    return result;
}

std::optional<Expr> try_algebraic_linear_sub(const Expr& expr, const Expr& var) {
    if (!depends_on(expr, var)) return std::nullopt;

    // Split the integrand into factors. Build with a push_back loop rather than
    // vector::assign(span): under -O3 the latter inlines a tail-destroy path
    // that trips gcc's -Wnull-dereference (a false positive on the element
    // shared_ptr destructor).
    std::vector<Expr> factors;
    if (expr->type_id() == TypeId::Mul) {
        auto sp = expr->args();
        factors.reserve(sp.size());
        for (const auto& f : sp) factors.push_back(f);
    } else {
        factors.push_back(expr);
    }

    // Find exactly one (affine)^(non-integer rational) factor; the remaining
    // factors must together form a polynomial in var. A second fractional power
    // of var lands in `rest` and fails the polynomial test below, so we bail.
    std::optional<std::pair<Expr, Expr>> lin;  // (a, b) of the affine base a·x+b
    Expr root_base;
    Expr r;  // the fractional exponent
    bool found = false;
    std::vector<Expr> rest;
    for (const auto& f : factors) {
        if (!found && f->type_id() == TypeId::Pow
            && f->args()[1]->type_id() == TypeId::Rational
            && depends_on(f->args()[0], var)) {
            auto aff = as_affine(f->args()[0], var);
            if (aff.has_value() && !(aff->first == S::Zero())) {
                lin = aff;
                root_base = f->args()[0];
                r = f->args()[1];
                found = true;
                continue;
            }
        }
        rest.push_back(f);
    }
    if (!found) return std::nullopt;

    Expr poly_part = rest.empty() ? S::One() : mul(rest);
    if (!is_polynomial_in(poly_part, var)) return std::nullopt;

    // Substitute u = a·x+b ⟹ x = (u−b)/a. The integrand (poly_part)·u^r becomes
    // a polynomial in u times u^r, i.e. Σ cₖ·u^(k+r), which integrates to
    // Σ cₖ·u^(k+r+1)/(k+r+1). r is a non-integer Rational, so k+r+1 ≠ 0 always.
    const Expr& a = lin->first;
    const Expr& b = lin->second;
    auto u = symbol("__alg_lin_u");
    Expr xsub = (u - b) / a;
    Expr poly_u = expand(subs(poly_part, var, xsub));
    if (depends_on(poly_u, var)) return std::nullopt;  // substitution didn't clear var

    Poly P(poly_u, u);
    const auto& cs = P.coeffs();
    std::vector<Expr> terms;
    terms.reserve(cs.size());
    for (std::size_t k = 0; k < cs.size(); ++k) {
        if (cs[k] == S::Zero()) continue;
        Expr exp_new = integer(static_cast<long>(k)) + r + S::One();  // k + r + 1
        terms.push_back(cs[k] * pow(u, exp_new) / exp_new);
    }
    Expr result_u = add(std::move(terms)) / a;
    // Back-substitute u = a·x+b using the original base expression.
    return simplify(subs(result_u, u, root_base));
}

}  // namespace

Expr integrate(const Expr& expr, const Expr& var,
              const Expr& lower, const Expr& upper) {
    auto antider = integrate(expr, var);
    return subs(antider, var, upper) - subs(antider, var, lower);
}

std::optional<Expr> manualintegrate(const Expr& expr, const Expr& var) {
    Expr result = integrate(expr, var);
    // Detect the Integral(_, _) marker integrate() emits on failure.
    if (result->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*result);
        if (fn.name() == "Integral") return std::nullopt;
    }
    // Failure can also be wrapped — Mul(c, Integral(...)). Walk one
    // level: if any args are Integral markers, treat as failure.
    if (result->type_id() == TypeId::Mul
        || result->type_id() == TypeId::Add) {
        for (const auto& a : result->args()) {
            if (a->type_id() == TypeId::Function) {
                const auto& fn = static_cast<const Function&>(*a);
                if (fn.name() == "Integral") return std::nullopt;
            }
        }
    }
    return result;
}

}  // namespace sympp
