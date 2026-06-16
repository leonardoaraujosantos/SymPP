#include <sympp/solvers/solve.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

#include <gmpxx.h>
#include <mpfr.h>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/add.hpp>
#include <sympp/core/boolean.hpp>
#include <sympp/core/expand.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/symbol.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/special.hpp>  // lambertw
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/polys/poly.hpp>
#include <sympp/polys/rootof.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp {

namespace {

// Polynomial-only root finder. This is the original solve() body; it is
// kept separate so solveset()'s internal fallback can use it WITHOUT
// re-entering the transcendental solveset path below (which would recurse).
// Largest polynomial degree for which the RootOf supplement runs factor_list.
// Kronecker factorization is exponential in degree, so cap it to keep solve()
// responsive; beyond this, an unsolvable polynomial still returns its radical
// roots (possibly empty) rather than risking a pathological factorization.
constexpr std::size_t kRootOfFactorDegreeCap = 12;

[[nodiscard]] std::vector<Expr> solve_poly(const Expr& expr, const Expr& var) {
    Poly p(expr, var);
    std::vector<Expr> roots = p.roots();
    // The closed-form solver covers factors of degree <= 4 (Cardano/Ferrari) plus
    // rational roots of higher-degree factors. An irreducible factor of degree
    // >= 5 is not solvable by radicals, so roots() leaves it unrepresented and
    // solve() would return an empty list — implying "no solutions" for e.g.
    // x^5 - x - 1, which has five. When roots are missing, represent each such
    // factor's real roots as RootOf (rendered CRootOf(poly, k), as SymPy does).
    // Complex roots of these factors are not yet representable (RootOf is
    // real-root-only) and are left out — a known partial-parity limitation.
    //
    // Gated on degree to avoid paying for factorization on the common low-degree
    // path (and to bound the exponential Kronecker cost): only attempt this when
    // the polynomial is degree 5..cap and the radical solver left roots missing.
    const std::size_t pdeg = p.degree();
    if (pdeg < 5 || pdeg > kRootOfFactorDegreeCap || roots.size() >= pdeg) {
        return roots;
    }
    try {
        const FactorList fl = factor_list(p);
        for (const auto& [factor, mult] : fl.factors) {
            (void)mult;  // distinct roots only; multiplicity adds no new solution
            const std::size_t deg = factor.degree();
            if (deg < 5) continue;  // already handled by the radical solver
            const Expr fexpr = factor.as_expr();
            for (std::size_t k = 0; k < deg; ++k) {
                Expr r = root_of(fexpr, var, k);
                // try_evalf returns nullopt past the last real root: stop there.
                if (!static_cast<const RootOf&>(*r).try_evalf(20).has_value()) {
                    break;
                }
                if (std::none_of(roots.begin(), roots.end(),
                                 [&](const Expr& u) { return u == r; })) {
                    roots.push_back(std::move(r));
                }
            }
        }
    } catch (...) {
        // Factorization/RootOf failure: fall back to the radical roots alone.
    }
    return roots;
}

// Does `expr` contain a function (log, exp, sin, …) that depends on `var`?
// Used to decide whether the transcendental solveset fallback is worth a try.
[[nodiscard]] bool has_function_of_var(const Expr& e, const Expr& var) {
    if (e->type_id() == TypeId::Function && has(e, var)) return true;
    for (const auto& a : e->args()) {
        if (has_function_of_var(a, var)) return true;
    }
    return false;
}

// A non-integer power of var (e.g. sqrt(x) = x^(1/2)) that the polynomial path
// cannot handle but the invert chain can (g^p = c → g = c^(1/p)).
[[nodiscard]] bool has_radical_of_var(const Expr& e, const Expr& var) {
    if (e->type_id() == TypeId::Pow && has(e->args()[0], var)
        && e->args()[1]->type_id() == TypeId::Rational) {
        return true;
    }
    for (const auto& a : e->args()) {
        if (has_radical_of_var(a, var)) return true;
    }
    return false;
}

// For sin/cos, a real solution of f(x) = c exists only when |c| <= 1. Reject a
// c that is a real number outside [-1, 1] (those have purely complex solutions,
// which the inverse functions would leave as unevaluated asin/acos). A
// non-numeric c (a symbolic parameter, or an in-range irrational like 1/3) is
// accepted — SymPy returns asin(c) there too.
[[nodiscard]] bool trig_value_in_range(const Expr& c) {
    Expr v = evalf(c, 30);
    if (v->type_id() != TypeId::Float) return true;  // symbolic → accept
    try {
        return std::fabs(std::stod(v->str())) <= 1.0 + 1e-9;
    } catch (...) {
        return true;
    }
}

// Strict form of trig_value_in_range for the multi-angle resultant path: the
// cosine root must be a concrete REAL number in [−1, 1]. A complex root (from a
// negative-discriminant resultant, e.g. sin x + cos x − 5) evaluates to a
// non-Float and is rejected, where trig_value_in_range would accept it as
// "symbolic".
[[nodiscard]] bool is_unit_real(const Expr& r) {
    Expr v = evalf(r, 40);
    if (v->type_id() != TypeId::Float) return false;
    try {
        return std::fabs(std::stod(v->str())) <= 1.0 + 1e-9;
    } catch (...) {
        return false;
    }
}

// Append the representative roots of f(B*var) = c over one principal period
// (var = theta/B), deduplicating into `roots`. Mirrors SymPy's `solve`, which
// inverts the inner function and divides by B rather than enumerating every x
// in [0, 2pi). For sin/cos an out-of-range numeric c contributes nothing.
// Append the representative roots of f(B·var + P) = c, where the inner argument
// is affine in var: B·var + P with B, P var-free and B ≠ 0. Each principal angle
// θ of f gives var = (θ − P)/B. P defaults to 0 (a pure B·var argument).
void append_trig_roots(FunctionId id, const Expr& c, const Expr& B,
                       const Expr& var, std::vector<Expr>& roots,
                       const Expr& P = S::Zero()) {
    std::vector<Expr> thetas;
    switch (id) {
        case FunctionId::Sin:
            if (!trig_value_in_range(c)) return;
            thetas = {asin(c), simplify(S::Pi() - asin(c))};
            break;
        case FunctionId::Cos:
            if (!trig_value_in_range(c)) return;
            thetas = {acos(c), simplify(integer(2) * S::Pi() - acos(c))};
            break;
        case FunctionId::Tan:
            thetas = {atan(c)};
            break;
        default:
            return;
    }
    for (auto& th : thetas) {
        Expr r = simplify((th - P) * pow(B, integer(-1)));
        if (has(r, var)) continue;
        if (std::none_of(roots.begin(), roots.end(),
                         [&](const Expr& u) { return u == r; })) {
            roots.push_back(std::move(r));
        }
    }
}

// Decompose a trig argument that is affine in var: arg = B·var + P with B, P
// var-free and B ≠ 0. Returns {B, P} or nullopt if arg is not affine (degree > 1
// in var, var-dependent coefficients, or no var). Handles a bare B·var (P = 0)
// and an additive phase such as 2·x + π/3.
[[nodiscard]] std::optional<std::pair<Expr, Expr>>
affine_arg(const Expr& arg, const Expr& var) {
    if (!has(arg, var)) return std::nullopt;
    Expr B, P;
    try {
        Poly p(expand(arg), var);
        if (p.degree() != 1) return std::nullopt;
        const auto& cf = p.coeffs();
        P = cf[0];
        B = cf[1];
    } catch (...) {
        return std::nullopt;
    }
    if (has(B, var) || has(P, var) || B == S::Zero()) return std::nullopt;
    return std::make_pair(B, P);
}

// Solve A*f(B*var) + C = 0 for f in {sin, cos, tan} over one principal period,
// returning representative roots. This mirrors SymPy's `solve`, which inverts
// the inner function and divides by B (rather than enumerating every x in
// [0, 2pi)). Returns nullopt unless expr is a single trig term, linear in var
// with a var-free coefficient and no additive phase. Periodic/infinite sets
// remain the domain of solveset.
[[nodiscard]] std::optional<std::vector<Expr>>
solve_trig(const Expr& expr, const Expr& var) {
    // Split off the var-free additive constant C; the rest must be one term.
    Expr dep, cst = S::Zero();
    if (expr->type_id() == TypeId::Add) {
        std::vector<Expr> dv, fv;
        for (const auto& t : expr->args())
            (has(t, var) ? dv : fv).push_back(t);
        if (dv.size() != 1) return std::nullopt;
        dep = dv[0];
        if (!fv.empty()) cst = (fv.size() == 1) ? fv[0] : add(fv);
    } else {
        dep = expr;
    }
    // dep = coeff * core, coeff var-free, core the single var-dependent factor.
    Expr coeff = S::One();
    Expr core;
    if (dep->type_id() == TypeId::Mul) {
        std::vector<Expr> cf, vf;
        for (const auto& f : dep->args()) (has(f, var) ? vf : cf).push_back(f);
        // Zero-product: a product of var-dependent factors vanishes iff some
        // factor does. Solve each factor (recursively, via full solve so
        // polynomial factors like x are handled too) and union the roots.
        if (vf.size() >= 2 && cst == S::Zero()) {
            std::vector<Expr> roots;
            for (const auto& fac : vf) {
                for (auto& r : solve(fac, var)) {
                    if (has(r, var)) continue;
                    if (std::none_of(roots.begin(), roots.end(),
                                     [&](const Expr& u) { return u == r; })) {
                        roots.push_back(std::move(r));
                    }
                }
            }
            return roots;
        }
        if (vf.size() != 1) return std::nullopt;
        core = vf[0];
        coeff = cf.empty() ? Expr{S::One()} : mul(cf);
    } else {
        core = dep;
    }
    // core is f(arg), or f(arg)^n with n a positive integer. The power form is
    // only solved homogeneously: f(arg)^n = 0 reduces to f(arg) = 0.
    Expr fexpr;
    if (core->type_id() == TypeId::Function) {
        fexpr = core;
    } else if (core->type_id() == TypeId::Pow
               && core->args()[0]->type_id() == TypeId::Function
               && core->args()[1]->type_id() == TypeId::Integer
               && !static_cast<const Integer&>(*core->args()[1]).is_negative()
               && core->args()[1] != S::Zero()) {
        if (cst != S::Zero()) return std::nullopt;   // only f(arg)^n = 0
        fexpr = core->args()[0];
    } else {
        return std::nullopt;
    }
    const FunctionId id = static_cast<const Function&>(*fexpr).function_id();
    if (id != FunctionId::Sin && id != FunctionId::Cos
        && id != FunctionId::Tan) {
        return std::nullopt;
    }
    // arg must be affine in var: B*var + P, B var-free nonzero (P is the phase).
    const Expr& arg = fexpr->args()[0];
    auto bp = affine_arg(arg, var);
    if (!bp) return std::nullopt;
    const Expr& B = bp->first;
    const Expr& P = bp->second;
    // f(arg) = c, c = -C/A.
    Expr c = simplify(mul({S::NegativeOne(), cst, pow(coeff, integer(-1))}));
    std::vector<Expr> roots;   // var = (theta - P) / B, deduplicated.
    append_trig_roots(id, c, B, var, roots, P);
    return roots;
}

// Solve a polynomial in a single trig function — e.g. sin(x)^2 - 1,
// 2*cos(x)^2 - cos(x) - 1, sin(x)^2 - sin(x). Substitute u = f(B*var), solve
// the resulting polynomial in u, then invert each in-range root to
// representative angles over one principal period (matching SymPy's `solve`).
// Returns nullopt unless expr is a polynomial in exactly one trig atom
// f(B*var), f in {sin, cos, tan} and B var-free. This subsumes the single-term
// and homogeneous-power cases of solve_trig and additionally handles higher
// degrees; the structural zero-product path stays in solve_trig.
[[nodiscard]] std::optional<std::vector<Expr>>
solve_trig_poly(const Expr& expr, const Expr& var) {
    // Collect the distinct var-dependent trig atoms without descending into
    // their arguments (so sin(x) inside sin(x)^2 counts once).
    std::vector<Expr> atoms;
    std::function<void(const Expr&)> collect = [&](const Expr& e) {
        if (e->type_id() == TypeId::Function && has(e, var)) {
            const FunctionId fid =
                static_cast<const Function&>(*e).function_id();
            if (fid == FunctionId::Sin || fid == FunctionId::Cos
                || fid == FunctionId::Tan) {
                if (std::none_of(atoms.begin(), atoms.end(),
                                 [&](const Expr& a) { return a == e; })) {
                    atoms.push_back(e);
                }
                return;
            }
        }
        for (const auto& a : e->args()) collect(a);
    };
    collect(expr);
    if (atoms.size() != 1) return std::nullopt;
    const Expr g = atoms.front();
    const FunctionId id = static_cast<const Function&>(*g).function_id();
    // arg must be affine in var: B*var + P, B var-free nonzero (P is the phase).
    const Expr& arg = g->args()[0];
    auto bp = affine_arg(arg, var);
    if (!bp) return std::nullopt;
    const Expr& B = bp->first;
    const Expr& P = bp->second;
    // Substitute u for the trig atom; the remainder must be a polynomial in u
    // alone (a leftover var means a mixed poly·trig term such as x*sin(x)).
    Expr u = symbol("u_trig_subst");
    Expr eu = subs(expr, g, u);
    if (has(eu, var)) return std::nullopt;
    std::vector<Expr> cvals;
    try {
        Poly p(expand(eu), u);
        cvals = p.roots();
    } catch (...) {
        return std::nullopt;
    }
    if (cvals.empty()) return std::nullopt;
    std::vector<Expr> roots;
    for (auto& c : cvals) {
        if (has(c, var) || has(c, u)) continue;
        append_trig_roots(id, c, B, var, roots, P);
    }
    return roots;
}

// Sign of a numeric expression via floating-point evaluation: -1, 0, or +1.
// Returns 0 when the value is not a concrete real number (symbolic).
[[nodiscard]] int numeric_sign(const Expr& e) {
    Expr v = evalf(e, 40);
    if (v->type_id() != TypeId::Float) return 0;
    const auto& f = static_cast<const Number&>(*v);
    if (f.is_positive()) return 1;
    if (f.is_negative()) return -1;
    return 0;
}

// Add the representative angle θ ∈ (−π, π] with cos θ = c0, sin θ = s0 to the
// root list, mapped back through the affine argument var = (θ − P)/B. The angle
// is taken as ±acos(c0), the sign chosen from sin's sign — this keeps the clean
// acos special-angle forms (π/6, …) rather than a half-angle atan expression.
void append_angle_from_cos_sin(const Expr& c0, const Expr& s0, const Expr& B,
                               const Expr& P, const Expr& var,
                               std::vector<Expr>& roots) {
    Expr base = simplify(acos(c0));
    Expr theta = (numeric_sign(s0) < 0) ? simplify(mul({S::NegativeOne(), base}))
                                        : base;
    Expr r = simplify((theta - P) * pow(B, integer(-1)));
    if (has(r, var)) return;
    if (std::none_of(roots.begin(), roots.end(),
                     [&](const Expr& u) { return u == r; })) {
        roots.push_back(std::move(r));
    }
}

// Solve a trig equation whose terms carry different integer multiples of a
// common angle — e.g. sin(x) − cos(2x), cos(2x) + cos(x), sin(2x) − sin(x).
// expand_trig rewrites every multiple angle in terms of sin(θ), cos(θ) for the
// base angle θ, after which the equation is a polynomial in s = sin θ and
// c = cos θ. Reducing s² → 1 − c² brings it to P(c) + s·Q(c) = 0:
//   * Q ≡ 0  → P(c) = 0, a polynomial in cos θ (both ± angles per root);
//   * P ≡ 0  → s·Q(c) = 0, i.e. sin θ = 0 together with Q(c) = 0;
//   * else   → s = −P/Q with s² = 1 − c² gives P² − (1−c²)Q² = 0, a polynomial
//              in cos θ; each root fixes sin θ's sign, hence a unique angle.
// Returns representative roots over one period (matching SymPy's solve), or
// nullopt when the equation is not a single-base-angle trig polynomial.
[[nodiscard]] std::optional<std::vector<Expr>>
solve_trig_reduce(const Expr& expr, const Expr& var) {
    Expr e = expand(expand_trig(expr));
    // Collect the var-dependent trig atoms; they must all share one argument θ
    // and be sin/cos (after expansion tan would leave a denominator).
    std::vector<Expr> atoms;
    Expr theta;
    bool ok = true;
    std::function<void(const Expr&)> collect = [&](const Expr& g) {
        if (!ok) return;
        if (g->type_id() == TypeId::Function && has(g, var)) {
            const FunctionId fid =
                static_cast<const Function&>(*g).function_id();
            if (fid == FunctionId::Sin || fid == FunctionId::Cos) {
                if (!theta) theta = g->args()[0];
                else if (g->args()[0] != theta) { ok = false; return; }
                if (std::none_of(atoms.begin(), atoms.end(),
                                 [&](const Expr& a) { return a == g; })) {
                    atoms.push_back(g);
                }
                return;
            }
            ok = false;  // tan/other transcendental of var → not handled here
            return;
        }
        for (const auto& a : g->args()) collect(a);
    };
    collect(e);
    if (!ok || !theta || atoms.size() < 2) return std::nullopt;
    auto bp = affine_arg(theta, var);
    if (!bp) return std::nullopt;
    const Expr& B = bp->first;
    const Expr& Pphase = bp->second;

    // Substitute s, c for sin θ, cos θ. The result must be a pure polynomial.
    Expr s = symbol("s_trig_reduce");
    Expr c = symbol("c_trig_reduce");
    Expr es = subs(subs(e, sin(theta), s), cos(theta), c);
    if (has(es, var)) return std::nullopt;

    // Split into P(c) + s·Q(c) by reducing s² → 1 − c².
    Expr Pc = S::Zero();
    Expr Qc = S::Zero();
    Expr one_minus_c2 = integer(1) - pow(c, integer(2));
    try {
        Poly ps(expand(es), s);
        const auto& co = ps.coeffs();
        for (std::size_t k = 0; k < co.size(); ++k) {
            if (has(co[k], s) || has(co[k], var)) return std::nullopt;
            Expr half = pow(one_minus_c2, integer(static_cast<long>(k / 2)));
            if (k % 2 == 0) {
                Pc = Pc + co[k] * half;
            } else {
                Qc = Qc + co[k] * half;
            }
        }
    } catch (...) {
        return std::nullopt;
    }
    Pc = expand(Pc);
    Qc = expand(Qc);
    const bool p_zero = (simplify(Pc) == S::Zero());
    const bool q_zero = (simplify(Qc) == S::Zero());
    if (p_zero && q_zero) return std::nullopt;  // identically satisfied

    std::vector<Expr> roots;
    auto cos_roots = [&](const Expr& poly_in_c) -> bool {
        try {
            Poly pc(expand(poly_in_c), c);
            for (auto& r : pc.roots()) {
                if (has(r, c) || has(r, var)) continue;
                append_trig_roots(FunctionId::Cos, r, B, var, roots, Pphase);
            }
        } catch (...) {
            return false;
        }
        return true;
    };

    if (q_zero) {
        // P(cos θ) = 0 — both ± angles per in-range cosine root.
        if (!cos_roots(Pc)) return std::nullopt;
        return roots;
    }
    if (p_zero) {
        // sin θ = 0 ∪ Q(cos θ) = 0.
        append_trig_roots(FunctionId::Sin, S::Zero(), B, var, roots, Pphase);
        if (!cos_roots(Qc)) return std::nullopt;
        return roots;
    }
    // Mixed: P² − (1 − c²)·Q² = 0 in cos θ; each root fixes sin θ's sign. Build
    // the resultant by convolving coefficient vectors — robust and independent of
    // how the symbolic squares would canonicalize.
    try {
        auto conv = [](const std::vector<Expr>& a,
                       const std::vector<Expr>& b) {
            std::vector<Expr> r(a.empty() || b.empty() ? 0 : a.size() + b.size() - 1,
                                S::Zero());
            for (std::size_t i = 0; i < a.size(); ++i)
                for (std::size_t j = 0; j < b.size(); ++j)
                    r[i + j] = simplify(r[i + j] + a[i] * b[j]);
            return r;
        };
        std::vector<Expr> pa = Poly(Pc, c).coeffs();
        std::vector<Expr> qa = Poly(Qc, c).coeffs();
        std::vector<Expr> p2 = conv(pa, pa);
        std::vector<Expr> q2omc2 =
            conv(conv(qa, qa), {S::One(), S::Zero(), S::NegativeOne()});
        std::vector<Expr> rc(std::max(p2.size(), q2omc2.size()), S::Zero());
        for (std::size_t i = 0; i < p2.size(); ++i) rc[i] = p2[i];
        for (std::size_t i = 0; i < q2omc2.size(); ++i)
            rc[i] = simplify(rc[i] - q2omc2[i]);
        Poly pr(std::move(rc), c);
        for (auto& r : pr.roots()) {
            if (has(r, c) || has(r, var)) continue;
            if (!is_unit_real(r)) continue;  // reject complex / out-of-range
            Expr qval = simplify(subs(Qc, c, r));
            if (qval == S::Zero()) {
                // s·Q vanishes for any s at this cosine; it is a root only when
                // P(r) = 0 too, and then both ± angles satisfy the equation.
                if (simplify(subs(Pc, c, r)) == S::Zero()) {
                    append_trig_roots(FunctionId::Cos, r, B, var, roots, Pphase);
                }
                continue;
            }
            Expr s0 = simplify(mul({S::NegativeOne(), subs(Pc, c, r),
                                    pow(qval, integer(-1))}));
            append_angle_from_cos_sin(r, s0, B, Pphase, var, roots);
        }
    } catch (...) {
        return std::nullopt;
    }
    return roots;
}

// Solve a polynomial in a single exp or log atom — e.g. exp(x)^2-3*exp(x)+2,
// log(x)^2-1. Substitute u = g(B*var), solve the polynomial in u, then invert
// each root: exp(B*var)=c → var=log(c)/B (c=0 has no solution), log(B*var)=c →
// var=exp(c)/B. Matches SymPy's solve. Requires exactly one exp/log atom with
// arg = B*var (B var-free, no additive constant); a composite like exp(2x)
// alongside exp(x) is two distinct atoms and is left alone.
[[nodiscard]] std::optional<std::vector<Expr>>
solve_exp_log_poly(const Expr& expr, const Expr& var) {
    std::vector<Expr> atoms;
    std::function<void(const Expr&)> collect = [&](const Expr& e) {
        if (e->type_id() == TypeId::Function && has(e, var)) {
            const FunctionId fid =
                static_cast<const Function&>(*e).function_id();
            if (fid == FunctionId::Exp || fid == FunctionId::Log) {
                if (std::none_of(atoms.begin(), atoms.end(),
                                 [&](const Expr& a) { return a == e; })) {
                    atoms.push_back(e);
                }
                return;
            }
        }
        for (const auto& a : e->args()) collect(a);
    };
    collect(expr);
    if (atoms.size() != 1) return std::nullopt;
    const Expr g = atoms.front();
    const FunctionId id = static_cast<const Function&>(*g).function_id();
    const Expr& arg = g->args()[0];
    Expr B = simplify(arg * pow(var, integer(-1)));
    if (has(B, var) || simplify(B * var - arg) != S::Zero()) return std::nullopt;
    // exp(B·var)=c has B complex representatives when B≠1 (the period 2πi/B);
    // enumerating those is out of scope, so only the principal B=1 case is taken.
    // log is single-valued, so any B is fine there.
    if (id == FunctionId::Exp && !(B == S::One())) return std::nullopt;
    Expr u = symbol("u_explog_subst");
    Expr eu = subs(expr, g, u);
    if (has(eu, var)) return std::nullopt;
    std::vector<Expr> cvals;
    try {
        Poly p(expand(eu), u);
        cvals = p.roots();
    } catch (...) {
        return std::nullopt;
    }
    if (cvals.empty()) return std::nullopt;
    std::vector<Expr> roots;
    auto push = [&](Expr r) {
        if (has(r, var) || has(r, u)) return;
        if (std::none_of(roots.begin(), roots.end(),
                         [&](const Expr& w) { return w == r; })) {
            roots.push_back(std::move(r));
        }
    };
    for (auto& c : cvals) {
        if (has(c, var) || has(c, u)) continue;
        if (id == FunctionId::Exp) {
            if (c == S::Zero()) continue;  // exp(B*var) = 0 has no solution
            push(simplify(log(c) * pow(B, integer(-1))));
        } else {  // Log: log(B*var) = c → B*var = exp(c)
            push(simplify(exp(c) * pow(B, integer(-1))));
        }
    }
    return roots;
}

// Solve a sum of logarithms: Σ cᵢ·log(gᵢ(x)) + K = 0 with the cᵢ and K var-free.
// Combine via log(∏ gᵢ^cᵢ) = −K ⇒ ∏ gᵢ^cᵢ = exp(−K), solve that, and keep only
// roots in the log domain (every gᵢ(root) > 0). Mirrors SymPy's logcombine path:
// log(x)+log(x−1)=0 → x(x−1)=1 → x=(1+√5)/2 (the negative root is dropped).
[[nodiscard]] std::optional<std::vector<Expr>>
solve_log_sum(const Expr& expr, const Expr& var) {
    if (expr->type_id() != TypeId::Add) return std::nullopt;
    std::vector<std::pair<Expr, Expr>> logs;  // (g, coeff)
    Expr konst = S::Zero();
    auto as_log = [&](const Expr& t) -> std::optional<std::pair<Expr, Expr>> {
        // log(g) or (var-free coeff)·log(g).
        auto is_log = [&](const Expr& f) {
            return f->type_id() == TypeId::Function
                   && static_cast<const Function&>(*f).function_id()
                          == FunctionId::Log
                   && f->args().size() == 1;
        };
        if (is_log(t)) return std::make_pair(t->args()[0], Expr{S::One()});
        if (t->type_id() == TypeId::Mul) {
            Expr inner;
            std::vector<Expr> coeff;
            for (const auto& f : t->args()) {
                if (!inner && is_log(f)) {
                    inner = f->args()[0];
                    continue;
                }
                if (has(f, var)) return std::nullopt;  // non-log var factor
                coeff.push_back(f);
            }
            if (inner) return std::make_pair(inner, mul(std::move(coeff)));
        }
        return std::nullopt;
    };
    for (const auto& term : expr->args()) {
        if (!has(term, var)) {
            konst = konst + term;
            continue;
        }
        auto lp = as_log(term);
        if (!lp) return std::nullopt;  // a var-dependent non-log term
        logs.push_back(std::move(*lp));
    }
    if (logs.size() < 2) return std::nullopt;  // single log: handled elsewhere

    std::vector<Expr> factors;
    for (const auto& [g, c] : logs) factors.push_back(pow(g, c));
    Expr prod = mul(std::move(factors));
    Expr target = simplify(exp(mul(S::NegativeOne(), konst)));
    std::vector<Expr> raw = solve(simplify(prod - target), var);

    // A log argument is in-domain when it is positive. is_positive can't judge an
    // irrational like (1+√5)/2, so fall back to a numeric sign from evalf.
    auto positive_arg = [&](const Expr& g, const Expr& r) -> bool {
        Expr gv = simplify(subs(g, var, r));
        if (is_positive(gv) == true) return true;
        Expr fv = evalf(gv, 30);
        if (fv->type_id() != TypeId::Float) return false;
        try {
            return std::stod(fv->str()) > 1e-12;
        } catch (...) {
            return false;
        }
    };
    std::vector<Expr> out;
    for (auto& r : raw) {
        if (has(r, var)) continue;
        bool valid = true;
        for (const auto& [g, c] : logs) {
            if (!positive_arg(g, r)) {
                valid = false;
                break;
            }
        }
        if (valid
            && std::none_of(out.begin(), out.end(),
                            [&](const Expr& u) { return u == r; })) {
            out.push_back(std::move(r));
        }
    }
    return out;
}

// True if e is a concrete positive real number (numeric sign via evalf).
[[nodiscard]] bool numeric_positive(const Expr& e) {
    if (is_positive(e) == true) return true;
    Expr v = evalf(e, 30);
    if (v->type_id() != TypeId::Float) return false;
    try {
        return std::stod(v->str()) > 1e-12;
    } catch (...) {
        return false;
    }
}

// Solve a sum of constant-base exponentials Σ cᵢ·∏ⱼ aⱼ^(pⱼ·x+qⱼ) = 0. Each term
// reduces to coeffᵢ·exp(rateᵢ·x) with rateᵢ = Σ pⱼ·log(aⱼ) var-free. After
// combining equal rates: (A) when every rate is an integer multiple of a common
// r₀, substitute u = exp(r₀·x) → a polynomial in u (2^(2x)−5·2^x+4 → u²−5u+4);
// (B) with exactly two rates, d₁·exp(r₁x)+d₂·exp(r₂x)=0 ⇒
// x = log(−d₂/d₁)/(r₁−r₂) when −d₂/d₁ > 0 (2^x−3^x → 0).
[[nodiscard]] std::optional<std::vector<Expr>>
solve_const_base_exp_sum(const Expr& expr, const Expr& var) {
    // Only claim equations with a genuine constant-base power a^(…) (a ≠ e).
    // Pure exp(…) equations are left to solve_exp_sum, which keeps the complex
    // (period 2πi) roots SymPy enumerates for base e.
    bool has_const_pow = false;
    std::function<void(const Expr&)> scan = [&](const Expr& e) {
        if (e->type_id() == TypeId::Pow && has(e->args()[1], var)
            && !has(e->args()[0], var) && !(e->args()[0] == S::E())) {
            has_const_pow = true;
        }
        for (const auto& a : e->args()) scan(a);
    };
    scan(expr);
    if (!has_const_pow) return std::nullopt;

    // Normalize commensurate integer bases to a common base, so a quadratic in 2^x
    // written with 4^x = (2²)^x is recognized: if every constant integer base is a
    // power of the smallest one b, rewrite a^e = b^(k·e). (4^x − 2^x − 2 →
    // 2^(2x) − 2^x − 2, which the u = 2^x substitution below then closes.) This
    // sidesteps needing log(4)/log(2) = 2, which simplify() does not compute.
    {
        std::vector<long> ibases;
        std::function<void(const Expr&)> collect = [&](const Expr& e) {
            if (e->type_id() == TypeId::Pow && has(e->args()[1], var)
                && e->args()[0]->type_id() == TypeId::Integer) {
                const auto& z = static_cast<const Integer&>(*e->args()[0]);
                if (z.fits_long() && z.to_long() >= 2) ibases.push_back(z.to_long());
            }
            for (const auto& a : e->args()) collect(a);
        };
        collect(expr);
        auto power_of = [](long b, long base) -> long {  // k with base^k == b, else 0
            if (base < 2 || b < base) return b == 1 ? 0 : (b == base ? 1 : 0);
            long k = 0;
            long t = 1;
            while (t < b) { t *= base; ++k; }
            return t == b ? k : 0;
        };
        if (ibases.size() >= 2) {
            const long bmin = *std::min_element(ibases.begin(), ibases.end());
            const long bmax = *std::max_element(ibases.begin(), ibases.end());
            bool ok = bmax > bmin;
            for (long b : ibases) {
                if (power_of(b, bmin) == 0) { ok = false; break; }
            }
            if (ok) {
                ExprMap<Expr> repl;
                std::function<void(const Expr&)> rw = [&](const Expr& e) {
                    if (e->type_id() == TypeId::Pow && has(e->args()[1], var)
                        && e->args()[0]->type_id() == TypeId::Integer) {
                        const long b =
                            static_cast<const Integer&>(*e->args()[0]).to_long();
                        const long k = power_of(b, bmin);
                        if (b > bmin && k > 1) {
                            repl.emplace(e, pow(integer(bmin),
                                               mul(integer(k), e->args()[1])));
                        }
                    }
                    for (const auto& a : e->args()) rw(a);
                };
                rw(expr);
                if (!repl.empty()) {
                    Expr norm = xreplace(expr, repl);
                    if (!(norm == expr)) {
                        return solve_const_base_exp_sum(norm, var);
                    }
                }
            }
        }
    }

    auto term_rate =
        [&](const Expr& term) -> std::optional<std::pair<Expr, Expr>> {
        Expr coeff = S::One();
        Expr rate = S::Zero();
        std::vector<Expr> factors;
        if (term->type_id() == TypeId::Mul) {
            for (const auto& f : term->args()) factors.push_back(f);
        } else {
            factors.push_back(term);
        }
        for (const auto& f : factors) {
            if (!has(f, var)) {
                coeff = mul(coeff, f);
                continue;
            }
            Expr base;
            Expr exponent;
            if (f->type_id() == TypeId::Pow) {
                base = f->args()[0];
                exponent = f->args()[1];
            } else if (f->type_id() == TypeId::Function
                       && static_cast<const Function&>(*f).function_id()
                              == FunctionId::Exp) {
                base = S::E();
                exponent = f->args()[0];
            } else {
                return std::nullopt;
            }
            if (has(base, var) || is_positive(base) != std::optional<bool>{true})
                return std::nullopt;
            Expr p;
            Expr q;
            try {
                Poly pe(expand(exponent), var);
                if (pe.degree() > 1) return std::nullopt;
                const auto& cf = pe.coeffs();
                p = cf.size() > 1 ? cf[1] : Expr{S::Zero()};
                q = !cf.empty() ? cf[0] : Expr{S::Zero()};
            } catch (...) {
                return std::nullopt;
            }
            if (has(p, var) || has(q, var)) return std::nullopt;
            coeff = mul(coeff, pow(base, q));
            rate = add(rate, mul(p, log(base)));
        }
        return std::make_pair(simplify(coeff), simplify(rate));
    };

    std::vector<Expr> terms;
    if (expr->type_id() == TypeId::Add) {
        for (const auto& t : expr->args()) terms.push_back(t);
    } else {
        terms.push_back(expr);
    }
    std::vector<std::pair<Expr, Expr>> cr;  // (coeff, rate), combined by rate
    for (const auto& t : terms) {
        auto pr = term_rate(t);
        if (!pr) return std::nullopt;
        bool merged = false;
        for (auto& [c, r] : cr) {
            if (r == pr->second) {
                c = simplify(c + pr->first);
                merged = true;
                break;
            }
        }
        if (!merged) cr.emplace_back(pr->first, pr->second);
    }
    std::erase_if(cr, [](const auto& e) { return e.first == S::Zero(); });
    if (cr.size() < 2) return std::nullopt;

    Expr r0;
    for (const auto& [c, r] : cr) {
        if (!(r == S::Zero())) {
            r0 = r;
            break;
        }
    }
    // (A) commensurate rates → polynomial in u = exp(r₀·x).
    if (r0) {
        std::vector<std::pair<long, Expr>> pows;  // (exponent, coeff)
        bool commensurate = true;
        for (const auto& [c, r] : cr) {
            Expr ratio = simplify(r * pow(r0, integer(-1)));
            if (ratio->type_id() != TypeId::Integer
                || !static_cast<const Integer&>(*ratio).fits_long()) {
                commensurate = false;
                break;
            }
            pows.emplace_back(static_cast<const Integer&>(*ratio).to_long(), c);
        }
        if (commensurate) {
            auto u = symbol("__cbexp_u");
            std::vector<Expr> uterms;
            for (const auto& [nexp, c] : pows)
                uterms.push_back(mul(c, pow(u, integer(nexp))));
            std::vector<Expr> out;
            for (auto& ur : solve(add(std::move(uterms)), u)) {
                if (has(ur, u) || has(ur, var)) continue;
                if (numeric_positive(ur)) {
                    out.push_back(simplify(log(ur) * pow(r0, integer(-1))));
                }
            }
            return out;
        }
    }
    // (B) two incommensurate rates → ratio method.
    if (cr.size() == 2) {
        Expr ratio = simplify(
            mul({S::NegativeOne(), cr[1].first, pow(cr[0].first, integer(-1))}));
        Expr denom = simplify(cr[0].second - cr[1].second);
        if (numeric_positive(ratio) && !(denom == S::Zero())) {
            return std::vector<Expr>{
                simplify(log(ratio) * pow(denom, integer(-1)))};
        }
        return std::vector<Expr>{};
    }
    return std::nullopt;
}

// Solve a (Laurent) polynomial in exp(var) whose terms are exp(m·var) for
// integer m — e.g. exp(x)+exp(-x)-2, exp(2x)-3·exp(x)+2. Substitute u = exp(var)
// (so exp(m·var) → uᵐ), solve the resulting equation in u (the rational/poly
// path clears any negative powers), and invert each root via x = log(u).
// Handles the multi-atom and scaled-exponent cases that solve_exp_log_poly
// (single atom, B=1) leaves alone; complex roots are kept where SymPy keeps them
// (exp(2x)=1 → {0, iπ}).
[[nodiscard]] std::optional<std::vector<Expr>>
solve_exp_sum(const Expr& expr, const Expr& var) {
    std::vector<Expr> nodes;  // exp(m·var) atoms
    bool ok = true;
    std::function<void(const Expr&)> collect = [&](const Expr& e) {
        if (!ok || !has(e, var)) return;
        if (e->type_id() == TypeId::Function
            && static_cast<const Function&>(*e).function_id()
                   == FunctionId::Exp) {
            const Expr& arg = e->args()[0];
            Expr m = simplify(arg * pow(var, integer(-1)));
            if (m->type_id() != TypeId::Integer || m == S::Zero()) {
                ok = false;  // exp of a non-integer-multiple argument
                return;
            }
            if (std::none_of(nodes.begin(), nodes.end(),
                             [&](const Expr& n) { return n == e; })) {
                nodes.push_back(e);
            }
            return;
        }
        for (const auto& a : e->args()) collect(a);
    };
    collect(expr);
    // solve_exp_log_poly already covers a single exp(var) atom with unit rate;
    // this handles the rest (multiple atoms, or a scaled exp(m·var)).
    if (!ok || nodes.empty()) return std::nullopt;
    Expr u = symbol("u_expsum_subst");
    ExprMap<Expr> repl;
    for (const auto& node : nodes) {
        Expr m = simplify(node->args()[0] * pow(var, integer(-1)));
        repl.emplace(node, pow(u, m));
    }
    Expr eu = xreplace(expr, repl);
    if (has(eu, var)) return std::nullopt;  // a bare var survived (e.g. x·eˣ)
    std::vector<Expr> uroots = solve(eu, u);
    std::vector<Expr> roots;
    for (const auto& ur : uroots) {
        if (has(ur, u) || has(ur, var) || ur == S::Zero()) continue;
        Expr r = simplify(log(ur));
        if (std::none_of(roots.begin(), roots.end(),
                         [&](const Expr& w) { return w == r; })) {
            roots.push_back(std::move(r));
        }
    }
    return roots;
}

// Numeric value of a constant Expr, or nullopt if it isn't a real number.
[[nodiscard]] std::optional<double> numeric_value(const Expr& c) {
    Expr v = evalf(c, 30);
    if (v->type_id() != TypeId::Float) return std::nullopt;
    try {
        return std::stod(v->str());
    } catch (...) {
        return std::nullopt;
    }
}

// Solve a polynomial in a single inverse trig/hyperbolic atom — e.g.
// asin(x)-1, atan(2x)-1, asinh(x)-2. Substitute u = g(B*var), solve the
// polynomial, then invert: g⁻¹(B*var)=c → B*var = g(c) → var = g(c)/B, where g
// is the forward function. A root c outside the inverse function's range yields
// no solution (asin(x)=2 → []), matching SymPy. Inverse functions are
// single-valued, so any B is fine.
[[nodiscard]] std::optional<std::vector<Expr>>
solve_inverse_func_poly(const Expr& expr, const Expr& var) {
    static constexpr double kPi = 3.14159265358979323846;
    static constexpr double kHalfPi = kPi / 2.0;
    static constexpr double kEps = 1e-9;
    std::vector<Expr> atoms;
    std::function<void(const Expr&)> collect = [&](const Expr& e) {
        if (e->type_id() == TypeId::Function && has(e, var)) {
            const FunctionId fid =
                static_cast<const Function&>(*e).function_id();
            if (fid == FunctionId::Asin || fid == FunctionId::Acos
                || fid == FunctionId::Atan || fid == FunctionId::Asinh
                || fid == FunctionId::Acosh || fid == FunctionId::Atanh) {
                if (std::none_of(atoms.begin(), atoms.end(),
                                 [&](const Expr& a) { return a == e; })) {
                    atoms.push_back(e);
                }
                return;
            }
        }
        for (const auto& a : e->args()) collect(a);
    };
    collect(expr);
    if (atoms.size() != 1) return std::nullopt;
    const Expr g = atoms.front();
    const FunctionId id = static_cast<const Function&>(*g).function_id();
    const Expr& arg = g->args()[0];
    Expr B = simplify(arg * pow(var, integer(-1)));
    if (has(B, var) || simplify(B * var - arg) != S::Zero()) return std::nullopt;
    Expr u = symbol("u_invfn_subst");
    Expr eu = subs(expr, g, u);
    if (has(eu, var)) return std::nullopt;
    std::vector<Expr> cvals;
    try {
        Poly p(expand(eu), u);
        cvals = p.roots();
    } catch (...) {
        return std::nullopt;
    }
    if (cvals.empty()) return std::nullopt;
    std::vector<Expr> roots;
    auto push = [&](Expr r) {
        if (has(r, var) || has(r, u)) return;
        if (std::none_of(roots.begin(), roots.end(),
                         [&](const Expr& w) { return w == r; })) {
            roots.push_back(std::move(r));
        }
    };
    // c is in range when known-numeric and within bounds. A symbolic c yields the
    // formal principal-branch inverse (matching SymPy: solve(atan(z)−a) = tan(a));
    // only a concrete out-of-range value is rejected. An unbounded range accepts
    // anything.
    auto in_range = [&](const Expr& c, double lo, double hi, bool bounded) {
        if (!bounded) return true;
        auto v = numeric_value(c);
        if (!v) return true;  // symbolic c: return the formal inverse
        return *v >= lo - kEps && *v <= hi + kEps;
    };
    for (auto& c : cvals) {
        if (has(c, var) || has(c, u)) continue;
        const Expr inv_b = pow(B, integer(-1));
        switch (id) {
            case FunctionId::Asin:
                if (in_range(c, -kHalfPi, kHalfPi, true))
                    push(simplify(sin(c) * inv_b));
                break;
            case FunctionId::Acos:
                if (in_range(c, 0.0, kPi, true))
                    push(simplify(cos(c) * inv_b));
                break;
            case FunctionId::Atan:
                if (in_range(c, -kHalfPi, kHalfPi, true))
                    push(simplify(tan(c) * inv_b));
                break;
            case FunctionId::Asinh:
                push(simplify(sinh(c) * inv_b));
                break;
            case FunctionId::Acosh: {
                // acosh range is [0, ∞): need c ≥ 0.
                auto v = numeric_value(c);
                if (v && *v >= -kEps) push(simplify(cosh(c) * inv_b));
                break;
            }
            case FunctionId::Atanh:
                push(simplify(tanh(c) * inv_b));
                break;
            default:
                return std::nullopt;
        }
    }
    return roots;
}

// If `e` is exp(var) or log(var), report which; used for the additive Lambert
// forms var + exp(var) + c and var + log(var) + c.
[[nodiscard]] std::optional<FunctionId> as_exp_or_log_of_var(const Expr& e,
                                                            const Expr& var) {
    if (e->type_id() != TypeId::Function) return std::nullopt;
    const auto& f = static_cast<const Function&>(*e);
    const FunctionId id = f.function_id();
    if ((id == FunctionId::Exp || id == FunctionId::Log)
        && f.args().size() == 1 && f.args()[0] == var) {
        return id;
    }
    return std::nullopt;
}

// Solve the elementary Lambert-W forms via the identity W(z)·e^(W(z)) = z:
//   product-exp:  a·var·exp(b·var) + c = 0   → W(−c·b/a)/b
//   product-log:  a·var·log(var)   + c = 0   → exp(W(−c/a))
//   additive-exp: var + exp(var)   + c = 0   → −c − W(e^(−c))
//   additive-log: var + log(var)   + c = 0   → W(e^(−c))
// (a, b, c var-free; the additive forms require unit coefficients and argument
// var.) Matches SymPy's LambertW results.
[[nodiscard]] std::optional<std::vector<Expr>>
solve_lambert(const Expr& expr, const Expr& var) {
    // Additive forms: var + {exp(var)|log(var)} + c, with a single bare `var`
    // term and a single exp/log(var) term (each coefficient 1).
    if (expr->type_id() == TypeId::Add) {
        Expr c = S::Zero();
        Expr a = S::Zero();  // coefficient of the bare-var term a·var
        bool have_var = false, have_trans = false, trans_is_exp = false;
        bool ok = true;
        for (const auto& t : expr->args()) {
            if (!has(t, var)) {
                c = add({c, t});
                continue;
            }
            // A bare-var term a·var with a var-free (a is recovered as t/var).
            Expr q = simplify(mul(t, pow(var, integer(-1))));
            if (!has(q, var) && !have_var) {
                a = q;
                have_var = true;
                continue;
            }
            // The single exp(var) or log(var) term (unit coefficient, argument var).
            if (auto id = as_exp_or_log_of_var(t, var); id && !have_trans) {
                have_trans = true;
                trans_is_exp = (*id == FunctionId::Exp);
                continue;
            }
            ok = false;
            break;
        }
        if (ok && have_var && have_trans && !(a == S::Zero())) {
            Expr inv_a = pow(a, integer(-1));
            Expr cc = simplify(mul(c, inv_a));  // c/a
            if (trans_is_exp) {
                // a·var + e^var + c = 0 → var = −W(e^(−c/a)/a) − c/a.
                Expr z = simplify(mul(exp(mul(S::NegativeOne(), cc)), inv_a));
                return std::vector<Expr>{simplify(
                    mul(S::NegativeOne(), lambertw(z)) + mul(S::NegativeOne(), cc))};
            }
            // a·var + log(var) + c = 0 → var = W(a·e^(−c))/a.
            Expr z = simplify(mul(a, exp(mul(S::NegativeOne(), c))));
            return std::vector<Expr>{simplify(mul(lambertw(z), inv_a))};
        }
    }

    // Product forms: a·var·exp(b·var) + c  or  a·var·log(var) + c.
    Expr cst = S::Zero();
    Expr dep;
    if (expr->type_id() == TypeId::Add) {
        std::vector<Expr> dv, fv;
        for (const auto& t : expr->args())
            (has(t, var) ? dv : fv).push_back(t);
        if (dv.size() != 1) return std::nullopt;
        dep = dv[0];
        if (!fv.empty()) cst = (fv.size() == 1) ? fv[0] : add(fv);
    } else {
        dep = expr;
    }
    if (dep->type_id() != TypeId::Mul) return std::nullopt;
    Expr a = S::One();
    Expr b;
    bool have_var = false, have_exp = false, have_log = false;
    for (const auto& f : dep->args()) {
        if (!has(f, var)) {
            a = mul({a, f});
            continue;
        }
        if (f == var) {
            if (have_var) return std::nullopt;
            have_var = true;
            continue;
        }
        if (f->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*f);
            if (fn.function_id() == FunctionId::Exp && !have_exp && !have_log) {
                const Expr& arg = f->args()[0];
                Expr bb = simplify(arg * pow(var, integer(-1)));
                if (has(bb, var) || simplify(bb * var - arg) != S::Zero()) {
                    return std::nullopt;  // exp argument is not linear b·var
                }
                b = bb;
                have_exp = true;
                continue;
            }
            if (fn.function_id() == FunctionId::Log && !have_exp && !have_log
                && f->args()[0] == var) {
                have_log = true;
                continue;
            }
        }
        return std::nullopt;  // another var-dependent factor
    }
    if (!have_var) return std::nullopt;
    if (have_exp) {
        if (b == S::Zero()) return std::nullopt;
        // a·var·e^(b·var) = −c → var = W(−c·b/a)/b.
        Expr z = simplify(mul({S::NegativeOne(), cst, b, pow(a, integer(-1))}));
        return std::vector<Expr>{simplify(lambertw(z) * pow(b, integer(-1)))};
    }
    if (have_log) {
        // a·var·log(var) = −c → var·log(var) = −c/a → var = exp(W(−c/a)).
        Expr z = simplify(mul({S::NegativeOne(), cst, pow(a, integer(-1))}));
        return std::vector<Expr>{simplify(exp(lambertw(z)))};
    }
    return std::nullopt;
}

// If value == base^k for a positive integer k (base ≥ 2, value ≥ 1, both
// integers), return k — so 2^x = 8 solves to the clean 3, not log(8)/log(2).
[[nodiscard]] std::optional<long> integer_power_of(const mpz_class& base,
                                                   mpz_class value) {
    if (base < 2 || value < 1) return std::nullopt;
    long k = 0;
    while (value > 1) {
        if (!mpz_divisible_p(value.get_mpz_t(), base.get_mpz_t())) {
            return std::nullopt;
        }
        value /= base;
        ++k;
    }
    return k;
}

// Solve A·a^(b·var) + C = 0 for a constant base a (a positive number ≠ 1):
// a^(b·var) = −C/A → var = log(−C/A)/(b·log a). When −C/A is an exact integer
// power of an integer base, the clean rational exponent is returned instead of
// the log ratio (2^x = 8 → 3). `a^x` is a Pow with a numeric base, not the exp
// function, so the transcendental branch never sees it.
[[nodiscard]] std::optional<std::vector<Expr>>
solve_const_base_exp(const Expr& expr, const Expr& var) {
    Expr C = S::Zero();
    Expr dep;
    if (expr->type_id() == TypeId::Add) {
        std::vector<Expr> dv, fv;
        for (const auto& t : expr->args())
            (has(t, var) ? dv : fv).push_back(t);
        if (dv.size() != 1) return std::nullopt;
        dep = dv[0];
        if (!fv.empty()) C = (fv.size() == 1) ? fv[0] : add(fv);
    } else {
        dep = expr;
    }
    // dep = A · a^(b·var): a var-free coeff A and one Pow(a, b·var) with a
    // var-free and the exponent depending on var.
    Expr A = S::One();
    Expr base, argexp;
    bool have_pow = false;
    std::vector<Expr> factors;
    if (dep->type_id() == TypeId::Mul) {
        for (const auto& f : dep->args()) factors.push_back(f);
    } else {
        factors.push_back(dep);
    }
    for (const auto& f : factors) {
        if (!has(f, var)) {
            A = mul({A, f});
            continue;
        }
        if (f->type_id() == TypeId::Pow && !has(f->args()[0], var)
            && has(f->args()[1], var)) {
            if (have_pow) return std::nullopt;
            base = f->args()[0];
            argexp = f->args()[1];
            have_pow = true;
            continue;
        }
        return std::nullopt;  // another var-dependent factor
    }
    if (!have_pow) return std::nullopt;
    // Restrict to an integer base ≥ 2 that is NOT a perfect power, with the bare
    // exponent var. For a perfect-power base (4 = 2²) or a scaled exponent (b≠1),
    // SymPy enumerates extra complex representatives (the smaller complex period),
    // which this single-root inversion would not reproduce; those stay unsolved.
    if (base->type_id() != TypeId::Integer || argexp != var) return std::nullopt;
    const mpz_class& a = static_cast<const Integer&>(*base).value();
    if (a < 2 || mpz_perfect_power_p(a.get_mpz_t()) != 0) return std::nullopt;
    Expr rhs = simplify(mul({S::NegativeOne(), C, pow(A, integer(-1))}));
    // a^var = 0 has no solution (a^y > 0 for a real base).
    if (rhs == S::Zero()) return std::nullopt;
    // Exact integer power: a^var = a^k → var = k.
    if (rhs->type_id() == TypeId::Integer) {
        if (auto k = integer_power_of(a, static_cast<const Integer&>(*rhs).value())) {
            return std::vector<Expr>{integer(*k)};
        }
    }
    // General: var = log(rhs)/log(a).
    return std::vector<Expr>{simplify(log(rhs) * pow(log(base), integer(-1)))};
}

// Solve a polynomial in a radical x^(1/d) — e.g. x − √x − 2 (a quadratic in √x).
// Collects every x^e (e rational) and the bare var, takes d = lcm of the
// exponent denominators, substitutes t = x^(1/d) (x^e → t^(e·d)), solves the
// resulting polynomial in t, and back-substitutes x = t^d. Each candidate is
// verified against the original equation so extraneous roots (√x = −1) are
// dropped. Matches SymPy.
[[nodiscard]] std::optional<std::vector<Expr>>
solve_radical_poly(const Expr& expr, const Expr& var) {
    // Gather the distinct var-power subexpressions and the exponent denominators.
    std::vector<Expr> nodes;       // bare var or Pow(var, rational)
    std::vector<mpq_class> exps;   // matching exponent
    bool ok = true;
    std::function<void(const Expr&)> collect = [&](const Expr& e) {
        if (!ok || !has(e, var)) return;
        if (e == var) {
            if (std::none_of(nodes.begin(), nodes.end(),
                             [&](const Expr& n) { return n == e; })) {
                nodes.push_back(e);
                exps.emplace_back(1);
            }
            return;
        }
        if (e->type_id() == TypeId::Pow && e->args()[0] == var) {
            const Expr& ex = e->args()[1];
            mpq_class q;
            if (ex->type_id() == TypeId::Integer) {
                q = mpq_class(static_cast<const Integer&>(*ex).value());
            } else if (ex->type_id() == TypeId::Rational) {
                q = static_cast<const Rational&>(*ex).value();
            } else {
                ok = false;  // var^(symbolic) — not a radical polynomial
                return;
            }
            if (std::none_of(nodes.begin(), nodes.end(),
                             [&](const Expr& n) { return n == e; })) {
                nodes.push_back(e);
                exps.push_back(q);
            }
            return;
        }
        for (const auto& a : e->args()) collect(a);
    };
    collect(expr);
    if (!ok || nodes.empty()) return std::nullopt;
    // d = lcm of the exponent denominators; nothing to do when all are integers.
    mpz_class d = 1;
    for (const auto& q : exps) {
        mpz_class den = q.get_den();
        mpz_class g;
        mpz_gcd(g.get_mpz_t(), d.get_mpz_t(), den.get_mpz_t());
        d = d / g * den;
    }
    if (d <= 1 || !d.fits_slong_p()) return std::nullopt;
    const long dl = d.get_si();
    // Substitute each var-power node with t^(e·d).
    Expr t = symbol("t_radical_subst");
    ExprMap<Expr> repl;
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        mpq_class ed = exps[i] * mpq_class(d);  // e·d is an integer
        repl.emplace(nodes[i], pow(t, integer(ed.get_num().get_si())));
    }
    Expr et = xreplace(expr, repl);
    if (has(et, var)) return std::nullopt;  // var survived (e.g. inside a sin)
    std::vector<Expr> troots;
    try {
        Poly p(expand(et), t);
        troots = p.roots();
    } catch (...) {
        return std::nullopt;
    }
    if (troots.empty()) return std::nullopt;
    // Back-substitute x = t^d and keep only candidates that satisfy the original.
    std::vector<Expr> roots;
    for (const auto& tr : troots) {
        if (has(tr, t) || has(tr, var)) continue;
        Expr cand = simplify(pow(tr, integer(dl)));
        if (has(cand, var) || has(cand, t)) continue;
        if (simplify(subs(expr, var, cand)) != S::Zero()) continue;  // extraneous
        if (std::none_of(roots.begin(), roots.end(),
                         [&](const Expr& u) { return u == cand; })) {
            roots.push_back(std::move(cand));
        }
    }
    return roots;
}

// Solve an equation with a single square root of a var-dependent radicand,
// √(g(x)) appearing linearly — e.g. √(x+1) − x + 1 = 0. Isolate the radical and
// square: writing the equation as A(x)·√(g) + B(x) = 0 (A, B radical-free), the
// squared form A²·g − B² = 0 is a plain polynomial equation; its roots are
// filtered back through the original to drop the extraneous ones introduced by
// squaring (the √(g) = +B/A branch). Complements solve_radical_poly, which only
// handles a polynomial in x^(1/d) of the bare variable.
[[nodiscard]] std::optional<std::vector<Expr>>
solve_radical_isolate(const Expr& expr, const Expr& var) {
    // Collect the distinct √(g) subexpressions (exponent exactly 1/2) of var.
    std::vector<Expr> sqrts;
    std::function<void(const Expr&)> collect = [&](const Expr& e) {
        if (!has(e, var)) return;
        if (e->type_id() == TypeId::Pow
            && e->args()[1]->type_id() == TypeId::Rational
            && static_cast<const Rational&>(*e->args()[1]).value()
                   == mpq_class(1, 2)) {
            if (std::none_of(sqrts.begin(), sqrts.end(),
                             [&](const Expr& s) { return s == e; })) {
                sqrts.push_back(e);
            }
            return;  // do not descend into the radicand
        }
        for (const auto& a : e->args()) collect(a);
    };
    collect(expr);
    if (sqrts.empty() || sqrts.size() > 2) return std::nullopt;

    // Keep only candidates that satisfy the *original* equation, dropping the
    // extraneous roots squaring introduces. Verify numerically — a symbolic check
    // can't denest forms like √((3−√5)/2) = (√5−1)/2 — and fall back to the
    // symbolic test only when the root doesn't evaluate to a concrete Float.
    auto satisfies = [&](const Expr& r) -> bool {
        Expr resid = simplify(subs(expr, var, r));
        if (resid == S::Zero()) return true;
        Expr v = evalf(resid, 40);
        if (v == S::Zero()) return true;  // collapsed to an exact 0
        if (v->type_id() == TypeId::Float) {
            try {
                return std::fabs(std::stod(v->str())) < 1e-20;
            } catch (...) {
                // inconclusive — reject (squaring only ever adds roots)
            }
        }
        return false;
    };
    auto verified_roots =
        [&](const Expr& poly_or_radical_eq) -> std::vector<Expr> {
        std::vector<Expr> roots;
        for (const auto& r : solve(poly_or_radical_eq, var)) {
            if (has(r, var) || !satisfies(r)) continue;
            if (std::none_of(roots.begin(), roots.end(),
                             [&](const Expr& u) { return u == r; })) {
                roots.push_back(r);
            }
        }
        return roots;
    };

    if (sqrts.size() == 2) {
        // Two radicals √g1, √g2: isolate one and square once, leaving a single
        // radical the recursive solve (this same path, size 1) then clears.
        // expr = A1·√g1 + A2·√g2 + P (A1, A2, P radical-free) ⇒
        //   A1·√g1 = −(A2·√g2 + P) ⇒ A1²·g1 = A2²·g2 + 2·A2·P·√g2 + P².
        Expr w1 = symbol("__rad_w1");
        Expr w2 = symbol("__rad_w2");
        Expr ew = subs(subs(expr, sqrts[0], w1), sqrts[1], w2);
        if (has_radical_of_var(ew, var)) return std::nullopt;
        Expr A1;
        Expr A2;
        Expr P;
        try {
            Poly p1(expand(ew), w1);
            if (p1.degree() != 1) return std::nullopt;
            A1 = p1.coeffs()[1];
            Poly p2(expand(p1.coeffs()[0]), w2);  // P + A2·w2
            if (p2.degree() != 1) return std::nullopt;
            A2 = p2.coeffs()[1];
            P = p2.coeffs()[0];
        } catch (const std::exception&) {
            return std::nullopt;
        }
        // A1, A2, P must be free of both radical placeholders (no √g1·√g2 cross
        // term, no radical-dependent coefficients).
        if (has(A1, w1) || has(A1, w2) || has(A2, w1) || has(A2, w2)
            || has(P, w1) || has(P, w2)) {
            return std::nullopt;
        }
        const Expr& g1 = sqrts[0]->args()[0];
        const Expr& g2 = sqrts[1]->args()[0];
        Expr new_eq = expand(mul(mul(A1, A1), g1) - mul(mul(A2, A2), g2)
                             - mul(P, P)
                             - mul({integer(2), A2, P, sqrts[1]}));
        if (new_eq == S::Zero()) return std::nullopt;
        auto roots = verified_roots(new_eq);
        if (roots.empty()) return std::nullopt;
        return roots;
    }

    const Expr& s = sqrts[0];
    const Expr& g = s->args()[0];
    // Linearize in the radical: √(g) → w. The result must be radical-free and
    // degree 1 in w, i.e. A·w + B with A, B radical-free polynomials in x.
    Expr w = symbol("__rad_w");
    Expr ew = subs(expr, s, w);
    if (has_radical_of_var(ew, var)) return std::nullopt;  // nested/other radical
    Expr A;
    Expr B;
    try {
        Poly p(expand(ew), w);
        if (p.degree() != 1) return std::nullopt;
        B = p.coeffs()[0];
        A = p.coeffs()[1];
    } catch (const std::exception&) {
        return std::nullopt;
    }

    // A·√(g) = −B  ⇒  A²·g − B² = 0.
    Expr eq = expand(mul(mul(A, A), g) - mul(B, B));
    if (eq == S::Zero()) return std::nullopt;  // underdetermined
    auto roots = verified_roots(eq);
    if (roots.empty()) return std::nullopt;
    return roots;
}

// Solve a*sin(B*var) + b*cos(B*var) + c = 0 (the linear-combination / R-method
// case) for representative roots over one principal period, matching SymPy's
// solve. Requires both sin and cos of the SAME argument B*var, with var-free
// coefficients a, b (not both zero). Homogeneous (c = 0) reduces to
// tan(B*var) = -b/a, a single root. Otherwise a*sin+b*cos = R*sin(B*var+phi)
// with R = sqrt(a^2+b^2), phi = atan2(b, a), so sin(B*var+phi) = -c/R yields two
// representatives; an out-of-range -c/R (|.|>1) has no real solution.
[[nodiscard]] std::optional<std::vector<Expr>>
solve_trig_linear(const Expr& expr, const Expr& var) {
    if (expr->type_id() != TypeId::Add) return std::nullopt;
    Expr a = S::Zero();   // coefficient of sin(theta)
    Expr b = S::Zero();   // coefficient of cos(theta)
    Expr c = S::Zero();   // var-free additive constant
    Expr theta;           // shared trig argument
    bool have_sin = false, have_cos = false;
    auto same_theta = [&](const Expr& t) {
        if (!theta) { theta = t; return true; }
        return theta == t;
    };
    for (const auto& term : expr->args()) {
        if (!has(term, var)) {
            c = (c == S::Zero()) ? term : add({c, term});
            continue;
        }
        // term = coeff * f(theta): split a var-free coeff from one trig factor.
        Expr coeff = S::One();
        Expr trig = term;
        if (term->type_id() == TypeId::Mul) {
            std::vector<Expr> cf, vf;
            for (const auto& f : term->args())
                (has(f, var) ? vf : cf).push_back(f);
            if (vf.size() != 1) return std::nullopt;
            trig = vf.front();
            coeff = cf.empty() ? Expr{S::One()} : mul(cf);
        }
        if (trig->type_id() != TypeId::Function) return std::nullopt;
        const FunctionId id = static_cast<const Function&>(*trig).function_id();
        if (id == FunctionId::Sin) {
            if (have_sin || !same_theta(trig->args()[0])) return std::nullopt;
            a = coeff;
            have_sin = true;
        } else if (id == FunctionId::Cos) {
            if (have_cos || !same_theta(trig->args()[0])) return std::nullopt;
            b = coeff;
            have_cos = true;
        } else {
            return std::nullopt;
        }
    }
    if (!have_sin || !have_cos) return std::nullopt;  // single-term: other paths
    // theta = B*var with B var-free and nonzero.
    Expr B = simplify(theta * pow(var, integer(-1)));
    if (has(B, var) || simplify(B * var - theta) != S::Zero()) return std::nullopt;
    std::vector<Expr> roots;
    if (c == S::Zero()) {
        // tan(theta) = -b/a → theta = atan(-b/a), one representative.
        Expr t = simplify(mul({S::NegativeOne(), b, pow(a, integer(-1))}));
        roots.push_back(simplify(atan(t) * pow(B, integer(-1))));
        return roots;
    }
    Expr R = simplify(pow(add({mul({a, a}), mul({b, b})}), rational(1, 2)));
    Expr cR = simplify(mul({S::NegativeOne(), c, pow(R, integer(-1))}));
    if (!trig_value_in_range(cR)) return roots;  // empty: no real solution
    Expr phi = atan2(b, a);
    for (const Expr& th : {Expr{asin(cR)}, Expr{simplify(S::Pi() - asin(cR))}}) {
        Expr r = simplify((th - phi) * pow(B, integer(-1)));
        if (has(r, var)) continue;
        if (std::none_of(roots.begin(), roots.end(),
                         [&](const Expr& u) { return u == r; })) {
            roots.push_back(std::move(r));
        }
    }
    return roots;
}

// Solve a rational equation N(var)/D(var) = 0 by clearing the denominator:
// combine into a single fraction, solve the numerator, then discard any root
// that vanishes the denominator (an extraneous pole). Returns nullopt unless the
// equation actually has a var-dependent denominator. Matches SymPy's solve.
[[nodiscard]] std::optional<std::vector<Expr>>
solve_rational(const Expr& expr, const Expr& var) {
    Expr t = together(expr);
    // Read base^(-n) (n a positive integer) as a denominator factor base^n.
    auto den_base = [](const Expr& f) -> std::optional<Expr> {
        if (f->type_id() == TypeId::Pow
            && f->args()[1]->type_id() == TypeId::Integer) {
            const auto& z = static_cast<const Integer&>(*f->args()[1]);
            if (z.is_negative() && z.fits_long()) {
                return pow(f->args()[0], integer(-z.to_long()));
            }
        }
        return std::nullopt;
    };
    Expr num = S::One();
    Expr den;
    if (auto d = den_base(t)) {
        den = *d;  // bare reciprocal 1/D
    } else if (t->type_id() == TypeId::Mul) {
        std::vector<Expr> nf, df;
        for (const auto& f : t->args()) {
            if (auto df_part = den_base(f)) {
                df.push_back(*df_part);
            } else {
                nf.push_back(f);
            }
        }
        if (df.empty()) return std::nullopt;  // no denominator
        den = mul(std::move(df));
        num = nf.empty() ? Expr{S::One()} : mul(std::move(nf));
    } else {
        return std::nullopt;
    }
    if (!has(den, var)) return std::nullopt;  // constant denominator — not ours
    std::vector<Expr> roots = solve(num, var);
    std::erase_if(roots, [&](const Expr& r) {
        return has(r, var) || simplify(subs(den, var, r)) == S::Zero();
    });
    return roots;
}

// A factor that can never be zero contributes no roots to a zero-product: exp(·)
// is nonzero over all of ℂ, as is any nonzero constant.
[[nodiscard]] bool is_never_zero(const Expr& f) {
    if (f->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*f);
        if (fn.function_id() == FunctionId::Exp) return true;
    }
    if (is_number(f) && !(f == S::Zero())) return true;
    return false;
}

// Zero-product: a product vanishes iff one of its factors does. Solve each factor
// (recursively, via full solve) and union the roots, skipping factors that can
// never be zero (exp(·)) and denominator factors (negative powers). Handles an
// explicit Mul, and an Add with a common factor — x²·eˣ − eˣ = eˣ·(x²−1) — that
// the polynomial path cannot see through (eˣ is not polynomial). Without this such
// equations returned [] and a Mul like eˣ·(x²−4) emitted a spurious zoo from
// solving eˣ = 0.
[[nodiscard]] std::optional<std::vector<Expr>>
solve_zero_product(const Expr& expr, const Expr& var) {
    std::vector<Expr> factors;  // factors whose =0 we solve and union

    // Split a term into base→exponent (positive integer powers folded together).
    auto term_factors = [](const Expr& t) {
        std::vector<std::pair<Expr, long>> fs;
        auto add_one = [&](const Expr& f) {
            Expr base = f;
            long e = 1;
            if (f->type_id() == TypeId::Pow
                && f->args()[1]->type_id() == TypeId::Integer) {
                const auto& z = static_cast<const Integer&>(*f->args()[1]);
                if (z.fits_long() && z.to_long() > 0) {
                    e = z.to_long();
                    base = f->args()[0];
                }
            }
            for (auto& p : fs) {
                if (p.first == base) { p.second += e; return; }
            }
            fs.push_back({base, e});
        };
        if (t->type_id() == TypeId::Mul) {
            for (const auto& f : t->args()) add_one(f);
        } else {
            add_one(t);
        }
        return fs;
    };

    if (expr->type_id() == TypeId::Mul) {
        for (const auto& f : expr->args()) factors.push_back(f);
    } else if (expr->type_id() == TypeId::Add) {
        const auto& terms = expr->args();
        if (terms.size() < 2) return std::nullopt;
        // Common factors = intersection of the per-term factor maps (min power).
        auto common = term_factors(terms[0]);
        for (std::size_t i = 1; i < terms.size() && !common.empty(); ++i) {
            auto tf = term_factors(terms[i]);
            std::vector<std::pair<Expr, long>> next;
            for (const auto& [b, e] : common) {
                for (const auto& [b2, e2] : tf) {
                    if (b2 == b) { next.push_back({b, std::min(e, e2)}); break; }
                }
            }
            common = std::move(next);
        }
        bool has_var_common = false;
        for (const auto& [b, e] : common) {
            if (has(b, var)) has_var_common = true;
        }
        if (!has_var_common) return std::nullopt;
        // Cofactor = Σ (term ÷ ∏ commonᵢ): subtract the common exponents per term.
        std::vector<Expr> cofactor_terms;
        for (const auto& t : terms) {
            auto tf = term_factors(t);
            std::vector<Expr> remaining;
            for (const auto& [b, e] : tf) {
                long ce = 0;
                for (const auto& [cb, cee] : common) {
                    if (cb == b) { ce = cee; break; }
                }
                const long left = e - ce;
                if (left > 0) {
                    remaining.push_back(left == 1 ? b : pow(b, integer(left)));
                }
            }
            cofactor_terms.push_back(remaining.empty() ? Expr{S::One()}
                                                       : mul(remaining));
        }
        for (const auto& [b, e] : common) factors.push_back(b);
        factors.push_back(add(cofactor_terms));
    } else {
        return std::nullopt;
    }

    if (factors.size() < 2) return std::nullopt;
    std::vector<Expr> roots;
    std::vector<Expr> denoms;  // bases of negative-power (denominator) factors
    bool progressed = false;
    for (const auto& f : factors) {
        if (!has(f, var)) continue;        // constant factor: no roots
        if (is_never_zero(f)) { progressed = true; continue; }
        // A denominator factor (negative power) is never a root, and its zeros are
        // poles that must be excluded from the numerator roots below.
        if (f->type_id() == TypeId::Pow
            && f->args()[1]->type_id() == TypeId::Integer
            && static_cast<const Integer&>(*f->args()[1]).is_negative()) {
            denoms.push_back(f->args()[0]);
            progressed = true;
            continue;
        }
        if (f == expr) return std::nullopt;  // no reduction — avoid infinite loop
        progressed = true;
        for (const auto& r : solve(f, var)) {
            if (std::none_of(roots.begin(), roots.end(),
                             [&](const Expr& u) { return u == r; })) {
                roots.push_back(r);
            }
        }
    }
    if (!progressed) return std::nullopt;
    // Drop poles: a numerator root that also vanishes a denominator factor is not
    // a solution ((x²−1)/(x−1) = 0 has root −1 only, not the cancelled pole 1).
    if (!denoms.empty()) {
        std::erase_if(roots, [&](const Expr& r) {
            return std::any_of(denoms.begin(), denoms.end(), [&](const Expr& d) {
                return simplify(subs(d, var, r)) == S::Zero();
            });
        });
    }
    return roots;
}

}  // namespace

std::vector<Expr> solve(const Expr& expr, const Expr& var) {
    // An equation Eq(lhs, rhs) reduces to solving lhs − rhs = 0 (matching
    // SymPy's solve(Eq(...))). Other relations (inequalities) describe a region,
    // not a discrete root list, so they don't fit this vector API.
    if (expr->type_id() == TypeId::Relational) {
        const auto& r = static_cast<const Relational&>(*expr);
        if (r.kind() == RelKind::Eq) return solve(r.lhs() - r.rhs(), var);
        return {};
    }
    // A product with a transcendental/radical factor (x·cos x, eˣ·(x²−4)) is
    // mis-read by solve_poly as a polynomial in var whose "coefficients" are the
    // functions, yielding only a partial root set. Zero-product over the factors
    // is complete, so prefer it when a function/radical of var is present.
    if (has_function_of_var(expr, var) || has_radical_of_var(expr, var)) {
        if (auto zp = solve_zero_product(expr, var); zp && !zp->empty()) {
            return *zp;
        }
    }
    // Expand first so a factored polynomial reaches the Poly machinery:
    // solve((x-1)*(x-2)) was empty because Poly couldn't build from the Mul.
    auto roots = solve_poly(expand(expr), var);
    // A genuine solution x = c must be free of x. solve_poly treats a
    // var-dependent "coefficient" (e.g. exp(x) in x·exp(x) − 1) as constant and
    // can hand back a rearrangement such as x = exp(x)**(-1); discard any
    // candidate that still contains var rather than return a non-solution. A
    // RootOf is exempt: it legitimately embeds the defining polynomial in var.
    std::erase_if(roots, [&](const Expr& r) {
        return r->type_id() != TypeId::RootOf && has(r, var);
    });
    // Deduplicate: solve_poly emits a root once per factor, so a repeated factor
    // ((x+2)² , x²·(x−1)) yields duplicates. SymPy's solve returns the distinct
    // solution set, so collapse structurally-equal roots (order preserved).
    {
        std::vector<Expr> distinct;
        for (auto& r : roots) {
            bool seen = false;
            for (const auto& u : distinct) {
                if (u == r) { seen = true; break; }
            }
            if (!seen) distinct.push_back(std::move(r));
        }
        roots = std::move(distinct);
    }
    if (!roots.empty()) return roots;
    // Zero-product: a Mul, or an Add with a common factor (x²·eˣ − eˣ), vanishes
    // iff a factor does. Solve each factor, skipping the never-zero ones (eˣ) —
    // this also removes the spurious zoo that solving eˣ = 0 would inject.
    if (auto zp = solve_zero_product(expr, var); zp) return *zp;
    // Rational equation (1/x + … = 0): clear the denominator and solve the
    // numerator, dropping pole roots. The polynomial path above can't build a
    // Poly from a negative-power term, so these reach here empty.
    if (auto rr = solve_rational(expr, var); rr && !rr->empty()) return *rr;
    // Constant-base exponential a^(b·x)=c (a^x is a Pow, not the exp function, so
    // it skips the transcendental branch below).
    if (auto ce = solve_const_base_exp(expr, var); ce && !ce->empty())
        return *ce;
    // Sums of constant-base exponentials (2^x − 3^x, 2^(2x) − 5·2^x + 4, …),
    // which a^x being a Pow (not the exp function) also keeps out of the
    // transcendental branch below.
    if (auto ces = solve_const_base_exp_sum(expr, var); ces && !ces->empty())
        return *ces;
    // The polynomial path can't see through log/exp/sinh/… so transcendental
    // equations like log(x) - 1 = 0 come back empty. solveset() has the
    // _invert chain; route through it and surface any finite solution set as
    // a root vector. (Infinite/periodic solution sets, e.g. sin(x)=0, stay
    // the domain of solveset and yield no finite vector here.)
    if (has_function_of_var(expr, var) || has_radical_of_var(expr, var)) {
        // Single-trig equations have an infinite (periodic) solution set, so
        // solveset returns an ImageSet and the FiniteSet branch below misses
        // them. Surface SymPy-style representative roots over one period first.
        if (auto tp = solve_trig_poly(expr, var); tp && !tp->empty())
            return *tp;
        if (auto el = solve_exp_log_poly(expr, var); el && !el->empty())
            return *el;
        if (auto ls = solve_log_sum(expr, var); ls && !ls->empty()) return *ls;
        if (auto es = solve_exp_sum(expr, var); es && !es->empty()) return *es;
        if (auto inv = solve_inverse_func_poly(expr, var); inv && !inv->empty())
            return *inv;
        if (auto lw = solve_lambert(expr, var); lw && !lw->empty()) return *lw;
        if (auto rp = solve_radical_poly(expr, var); rp && !rp->empty())
            return *rp;
        if (auto ri = solve_radical_isolate(expr, var); ri && !ri->empty())
            return *ri;
        if (auto tl = solve_trig_linear(expr, var); tl && !tl->empty())
            return *tl;
        if (auto tr = solve_trig(expr, var); tr && !tr->empty()) return *tr;
        if (auto tm = solve_trig_reduce(expr, var); tm && !tm->empty())
            return *tm;
        SetPtr s = solveset(expr, var);
        // Flatten a finite solution set — a FiniteSet, or a Union of finite sets
        // (e.g. solveset(|x−1|−2) = {3} ∪ {−1}). Anything containing a
        // non-finite component (ImageSet, Interval, …) is not a discrete root
        // list and is left empty.
        std::function<std::optional<std::vector<Expr>>(const SetPtr&)> flatten =
            [&](const SetPtr& set) -> std::optional<std::vector<Expr>> {
            if (!set) return std::nullopt;
            switch (set->kind()) {
                case SetKind::Empty:
                    return std::vector<Expr>{};
                case SetKind::FiniteSet:
                    return static_cast<const FiniteSet&>(*set).elements();
                case SetKind::Union: {
                    const auto& u = static_cast<const Union&>(*set);
                    auto a = flatten(u.lhs());
                    if (!a) return std::nullopt;
                    auto b = flatten(u.rhs());
                    if (!b) return std::nullopt;
                    a->insert(a->end(), b->begin(), b->end());
                    return a;
                }
                default:
                    return std::nullopt;
            }
        };
        if (auto elems = flatten(s)) {
            std::erase_if(*elems, [&](const Expr& r) { return has(r, var); });
            // Deduplicate (a union may repeat a shared root).
            std::vector<Expr> out;
            for (auto& r : *elems) {
                if (std::none_of(out.begin(), out.end(),
                                 [&](const Expr& u) { return u == r; })) {
                    out.push_back(std::move(r));
                }
            }
            return out;
        }
    }
    return roots;
}

std::vector<Expr> solve(const Expr& lhs, const Expr& rhs, const Expr& var) {
    return solve(lhs - rhs, var);
}

Matrix linsolve(const Matrix& A, const Matrix& b) {
    if (!A.is_square()) {
        throw std::invalid_argument("linsolve: A must be square");
    }
    if (b.rows() != A.rows()) {
        throw std::invalid_argument("linsolve: b dimension mismatch");
    }
    return A.inverse() * b;
}

// ----- solveset --------------------------------------------------------------

SetPtr solveset(const Expr& expr, const Expr& var) {
    return solveset(expr, var, reals());
}

namespace {

// Detect simple trig zero patterns and emit ImageSet representing the
// full periodic solution. Returns nullopt when expr isn't a recognized
// pattern.
//   sin(var)    → {n·π : n ∈ ℤ}
//   cos(var)    → {n·π + π/2 : n ∈ ℤ}
//   tan(var)    → {n·π : n ∈ ℤ}
//   sin(a·var)  → {n·π/a : n ∈ ℤ}
//   cos(a·var)  → {(2n+1)·π/(2a) : n ∈ ℤ}
[[nodiscard]] std::optional<SetPtr> trig_solveset(const Expr& expr,
                                                    const Expr& var) {
    if (expr->type_id() != TypeId::Function) return std::nullopt;
    const auto& fn = static_cast<const Function&>(*expr);
    if (fn.args().size() != 1) return std::nullopt;
    const Expr& inner = fn.args()[0];
    // Extract scaling: inner = a · var with a free of var.
    Expr a;
    if (inner == var) {
        a = S::One();
    } else if (inner->type_id() == TypeId::Mul) {
        std::vector<Expr> rest;
        bool found = false;
        for (const auto& f : inner->args()) {
            if (f == var) {
                if (found) return std::nullopt;
                found = true;
            } else if (has(f, var)) {
                return std::nullopt;
            } else {
                rest.push_back(f);
            }
        }
        if (!found) return std::nullopt;
        a = mul(rest);
    } else {
        return std::nullopt;
    }
    auto n = symbol("__n");
    Expr image;
    switch (fn.function_id()) {
        case FunctionId::Sin:
        case FunctionId::Tan:
            image = mul(n, S::Pi()) / a;
            break;
        case FunctionId::Cos:
            image = (mul(integer(2), n) + integer(1)) * S::Pi()
                    / (integer(2) * a);
            break;
        default:
            return std::nullopt;
    }
    return image_set(n, image, integers());
}

// Split `expr` (treated as expr = 0) into (var-dependent terms, var-free
// terms) over an Add node. For non-Add input, var-dep = expr and free = 0.
[[nodiscard]] std::pair<Expr, Expr>
split_var_free(const Expr& expr, const Expr& var) {
    if (expr->type_id() == TypeId::Add) {
        std::vector<Expr> dep, free_terms;
        for (const auto& t : expr->args()) {
            if (has(t, var)) dep.push_back(t);
            else free_terms.push_back(t);
        }
        Expr d = dep.empty() ? Expr{S::Zero()}
                              : (dep.size() == 1 ? dep[0] : add(dep));
        Expr f = free_terms.empty() ? Expr{S::Zero()}
                                       : (free_terms.size() == 1
                                              ? free_terms[0]
                                              : add(free_terms));
        return {d, f};
    }
    if (has(expr, var)) return {expr, Expr{S::Zero()}};
    return {Expr{S::Zero()}, expr};
}

// _invert(expr, var) — peel a known function from expr = 0 to reduce to
// either a polynomial / direct equation in `var`, or to an ImageSet for
// the periodic-trig case. Returns nullopt when no head is recognized.
//
// Handles:
//   log(g) = c     → g = exp(c)
//   exp(g) = c     → g = log(c)
//   sin(g) = c, cos(g) = c, tan(g) = c → ImageSet over ℤ when g == var
//   sinh(g)=c, tanh(g)=c → g = asinh(c) / atanh(c)
//   cosh(g) = c    → g = ±acosh(c)
//   |g| = c        → g = ±c
//
// For non-trivial g, recurses with solveset(g - target, var). For the
// periodic cases we currently only emit ImageSet when g == var; nested
// linear g would require parametric ImageSet substitution which is
// deferred-deep.
[[nodiscard]] std::optional<SetPtr> invert_solveset(const Expr& expr,
                                                       const Expr& var,
                                                       const SetPtr& domain);

// Forward declaration so invert_solveset can recurse via solveset.
SetPtr solveset_impl(const Expr& expr, const Expr& var, const SetPtr& domain);

[[nodiscard]] std::optional<SetPtr> invert_solveset(const Expr& expr,
                                                       const Expr& var,
                                                       const SetPtr& domain) {
    auto [var_dep, free] = split_var_free(expr, var);
    if (var_dep == var || var_dep == S::Zero()) return std::nullopt;

    // Radical: g^p = c with p a non-integer rational → g = c^(1/p). Integer
    // powers are deliberately left to the polynomial solver, which finds every
    // root (e.g. x² = 4 → ±2) rather than just the principal one.
    if (var_dep->type_id() == TypeId::Pow) {
        const Expr& base = var_dep->args()[0];
        const Expr& p = var_dep->args()[1];
        Expr rhs = mul(S::NegativeOne(), free);  // g^p = rhs
        // g^n = 0 (n a positive integer) ⟺ g = 0. Without this, sin(x)² = 0
        // came back EmptySet (the polynomial path can't see through sin), even
        // though its solution set is exactly that of sin(x) = 0.
        if (p->type_id() == TypeId::Integer
            && is_positive(p) == std::optional<bool>{true}
            && rhs == S::Zero()) {
            return solveset_impl(base, var, domain);
        }
        if (p->type_id() != TypeId::Rational) return std::nullopt;
        // Principal-branch convention (matches SymPy): a negative real rhs has
        // no solution for a fractional power (√x = −2, x^(1/3) = −2 → ∅).
        if (is_negative(rhs) == std::optional<bool>{true}) return empty_set();
        Expr inv = pow(rhs, pow(p, integer(-1)));  // rhs^(1/p)
        if (base == var) return finite_set({inv});
        return solveset_impl(base - inv, var, domain);
    }

    if (var_dep->type_id() != TypeId::Function) return std::nullopt;
    if (var_dep->args().size() != 1) return std::nullopt;
    const auto& fn = static_cast<const Function&>(*var_dep);
    const Expr& g = var_dep->args()[0];
    Expr c = -free;  // var_dep = c
    auto fid = fn.function_id();

    auto recurse_or_finite = [&](const Expr& target) -> SetPtr {
        if (g == var) return finite_set({target});
        return solveset_impl(g - target, var, domain);
    };

    auto n = symbol("__n");

    // If g is linear in var (g == a·var, a free of var, a ≠ 0) return a, so a
    // periodic-trig solution over a·var divides through by a to give the
    // solution over var (cos(2x) = 1 → x ∈ {nπ}). Returns nullopt otherwise.
    auto linear_coeff = [&](const Expr& gg) -> std::optional<Expr> {
        if (gg == var) return Expr{S::One()};
        if (gg->type_id() != TypeId::Mul) return std::nullopt;
        std::vector<Expr> rest;
        bool found = false;
        for (const auto& f : gg->args()) {
            if (f == var) {
                if (found) return std::nullopt;
                found = true;
            } else if (has(f, var)) {
                return std::nullopt;
            } else {
                rest.push_back(f);
            }
        }
        if (!found) return std::nullopt;
        return rest.empty() ? Expr{S::One()} : mul(std::move(rest));
    };
    // Periodic image over var: (base + period·n) / a, simplified.
    auto periodic = [&](const Expr& base, const Expr& period, const Expr& a) {
        return image_set(n, simplify((base + period * n) / a), integers());
    };
    switch (fid) {
        case FunctionId::Log:
            return recurse_or_finite(exp(c));
        case FunctionId::Exp:
            return recurse_or_finite(log(c));
        case FunctionId::Sinh:
            return recurse_or_finite(asinh(c));
        case FunctionId::Tanh:
            return recurse_or_finite(atanh(c));
        case FunctionId::Cosh: {
            Expr p = acosh(c);
            if (g == var) return finite_set({p, -p});
            return set_union(solveset_impl(g - p, var, domain),
                              solveset_impl(g + p, var, domain));
        }
        case FunctionId::Abs: {
            // |g| = c has real solutions only for c ≥ 0; a concrete negative c
            // (e.g. |x+1| = −2 from |x+1|+2 = 0) has none.
            if (is_negative(c) == true) return empty_set();
            if (g == var) return finite_set({c, -c});
            return set_union(solveset_impl(g - c, var, domain),
                              solveset_impl(g + c, var, domain));
        }
        case FunctionId::Sin: {
            // sin(a·var) = c → var ∈ {((-1)^n·asin(c) + n·π) / a : n ∈ ℤ}.
            if (auto a = linear_coeff(g)) {
                return periodic(pow(integer(-1), n) * asin(c), S::Pi(), *a);
            }
            return std::nullopt;
        }
        case FunctionId::Cos: {
            // cos(a·var) = c → var ∈ {(±acos(c) + 2nπ) / a}. Two ImageSets,
            // but when acos(c) ∈ {0, π} (c = ±1) the ± branches coincide, so a
            // single ImageSet suffices (cos(2x)=1 → {nπ}, not {nπ} ∪ {nπ}).
            if (auto a = linear_coeff(g)) {
                Expr ac = simplify(acos(c));
                Expr two_pi = integer(2) * S::Pi();
                if (ac == S::Zero() || ac == S::Pi()) {
                    return periodic(ac, two_pi, *a);
                }
                return set_union(periodic(ac, two_pi, *a),
                                 periodic(mul(S::NegativeOne(), ac), two_pi, *a));
            }
            return std::nullopt;
        }
        case FunctionId::Tan: {
            // tan(a·var) = c → var ∈ {(atan(c) + nπ) / a}.
            if (auto a = linear_coeff(g)) {
                return periodic(atan(c), S::Pi(), *a);
            }
            return std::nullopt;
        }
        default:
            return std::nullopt;
    }
}

// Internal solveset that can be reentered from the inverter.
SetPtr solveset_impl(const Expr& expr, const Expr& var, const SetPtr& domain) {
    if (auto trig = trig_solveset(expr, var); trig) {
        if (domain->kind() == SetKind::Reals) return *trig;
        return set_intersection(*trig, domain);
    }
    if (auto inv = invert_solveset(expr, var, domain); inv) {
        if (domain->kind() == SetKind::Reals) return *inv;
        return set_intersection(*inv, domain);
    }
    // Polynomial-only fallback. Must NOT call the public solve(), which would
    // re-enter solveset() for transcendental input and recurse without bound.
    auto roots = solve_poly(expr, var);
    if (roots.empty()) return empty_set();
    std::vector<Expr> kept;
    kept.reserve(roots.size());
    for (auto& r : roots) {
        // Reject rearrangements that still contain var (solve_poly can return
        // x = exp(x)**(-1) for x·exp(x) − 1 by treating exp(x) as a constant).
        if (has(r, var)) continue;
        auto m = domain->contains(r);
        if (m == std::optional<bool>{false}) continue;
        kept.push_back(std::move(r));
    }
    return finite_set(std::move(kept));
}

}  // namespace

SetPtr solveset(const Expr& expr, const Expr& var, const SetPtr& domain) {
    return solveset_impl(expr, var, domain);
}

// ----- nsolve --------------------------------------------------------------

namespace {

// Convert an Expr to a double via recursive numeric evaluation. Walks
// Add, Mul, Pow, and Function nodes, evalfs their leaves, and combines
// the results in double. Returns nullopt for symbolic-only inputs.
[[nodiscard]] std::optional<double> recursive_evalf_to_double(
    const Expr& e, int dps);

[[nodiscard]] std::optional<double> expr_to_double(const Expr& e, int dps) {
    Expr v = evalf(e, dps);
    if (v->type_id() == TypeId::Float) {
        try { return std::stod(v->str()); }
        catch (...) { return std::nullopt; }
    }
    return recursive_evalf_to_double(e, dps);
}

[[nodiscard]] std::optional<double> recursive_evalf_to_double(
    const Expr& e, int dps) {
    if (!e) return std::nullopt;
    switch (e->type_id()) {
        case TypeId::Integer:
        case TypeId::Rational:
        case TypeId::Float:
        case TypeId::NumberSymbol:
        case TypeId::RootOf: {
            Expr v = evalf(e, dps);
            if (v->type_id() != TypeId::Float) return std::nullopt;
            try { return std::stod(v->str()); }
            catch (...) { return std::nullopt; }
        }
        case TypeId::Add: {
            double sum = 0.0;
            for (const auto& a : e->args()) {
                auto v = recursive_evalf_to_double(a, dps);
                if (!v) return std::nullopt;
                sum += *v;
            }
            return sum;
        }
        case TypeId::Mul: {
            double prod = 1.0;
            for (const auto& a : e->args()) {
                auto v = recursive_evalf_to_double(a, dps);
                if (!v) return std::nullopt;
                prod *= *v;
            }
            return prod;
        }
        case TypeId::Pow: {
            auto b = recursive_evalf_to_double(e->args()[0], dps);
            auto x = recursive_evalf_to_double(e->args()[1], dps);
            if (!b || !x) return std::nullopt;
            return std::pow(*b, *x);
        }
        case TypeId::Function: {
            // Use evalf on the wrapped arg, then dispatch by name.
            // Simpler: substitute Float values into args and let evalf
            // handle the function. But that requires reconstruction;
            // use mpfr directly when possible.
            const auto& fn = static_cast<const Function&>(*e);
            if (fn.args().size() != 1) return std::nullopt;
            auto arg = recursive_evalf_to_double(fn.args()[0], dps);
            if (!arg) return std::nullopt;
            switch (fn.function_id()) {
                case FunctionId::Sin: return std::sin(*arg);
                case FunctionId::Cos: return std::cos(*arg);
                case FunctionId::Tan: return std::tan(*arg);
                case FunctionId::Exp: return std::exp(*arg);
                case FunctionId::Log: return std::log(*arg);
                case FunctionId::Abs: return std::fabs(*arg);
                default: return std::nullopt;
            }
        }
        default:
            return std::nullopt;
    }
}

}  // namespace

Expr nsolve(const Expr& expr, const Expr& var, const Expr& x0, int dps) {
    if (dps < 6) dps = 6;
    const int work_dps = dps + 10;
    const mpfr_prec_t work_prec = dps_to_prec(work_dps);
    Expr fprime = diff(expr, var);

    // Initial guess as MPFR.
    mpfr_t x_mpfr, step_mpfr, fv_mpfr, fpv_mpfr, tol_mpfr, abs_step;
    mpfr_init2(x_mpfr, work_prec);
    mpfr_init2(step_mpfr, work_prec);
    mpfr_init2(fv_mpfr, work_prec);
    mpfr_init2(fpv_mpfr, work_prec);
    mpfr_init2(tol_mpfr, work_prec);
    mpfr_init2(abs_step, work_prec);
    mpfr_set_d(tol_mpfr, 10.0, MPFR_RNDN);
    mpfr_pow_si(tol_mpfr, tol_mpfr, -dps, MPFR_RNDN);

    Expr x0_eval = evalf(x0, work_dps);
    if (x0_eval->type_id() != TypeId::Float) {
        // Try numeric coercion via str round-trip.
        try {
            mpfr_set_str(x_mpfr, x0->str().c_str(), 10, MPFR_RNDN);
        } catch (...) {
            mpfr_clear(x_mpfr); mpfr_clear(step_mpfr); mpfr_clear(fv_mpfr);
            mpfr_clear(fpv_mpfr); mpfr_clear(tol_mpfr); mpfr_clear(abs_step);
            throw std::runtime_error("nsolve: initial guess must be numeric");
        }
    } else {
        mpfr_set_str(x_mpfr, x0_eval->str().c_str(), 10, MPFR_RNDN);
    }

    const int max_iters = 200;
    bool converged = false;
    for (int i = 0; i < max_iters; ++i) {
        // Build a Float Expr from x_mpfr and substitute.
        char buf[256];
        mpfr_snprintf(buf, sizeof(buf), "%.*Re", work_dps, x_mpfr);
        Expr xe = make<Float>(std::string(buf), work_dps);
        Expr fv = evalf(subs(expr, var, xe), work_dps);
        Expr fpv = evalf(subs(fprime, var, xe), work_dps);
        if (fv->type_id() != TypeId::Float
            || fpv->type_id() != TypeId::Float) {
            mpfr_clear(x_mpfr); mpfr_clear(step_mpfr); mpfr_clear(fv_mpfr);
            mpfr_clear(fpv_mpfr); mpfr_clear(tol_mpfr); mpfr_clear(abs_step);
            throw std::runtime_error("nsolve: non-numeric step");
        }
        mpfr_set_str(fv_mpfr, fv->str().c_str(), 10, MPFR_RNDN);
        mpfr_set_str(fpv_mpfr, fpv->str().c_str(), 10, MPFR_RNDN);
        // step = fv / fpv
        if (mpfr_zero_p(fpv_mpfr)) {
            mpfr_clear(x_mpfr); mpfr_clear(step_mpfr); mpfr_clear(fv_mpfr);
            mpfr_clear(fpv_mpfr); mpfr_clear(tol_mpfr); mpfr_clear(abs_step);
            throw std::runtime_error("nsolve: derivative vanished");
        }
        mpfr_div(step_mpfr, fv_mpfr, fpv_mpfr, MPFR_RNDN);
        mpfr_sub(x_mpfr, x_mpfr, step_mpfr, MPFR_RNDN);
        mpfr_abs(abs_step, step_mpfr, MPFR_RNDN);
        if (mpfr_cmp(abs_step, tol_mpfr) < 0) { converged = true; break; }
    }
    if (!converged) {
        mpfr_clear(x_mpfr); mpfr_clear(step_mpfr); mpfr_clear(fv_mpfr);
        mpfr_clear(fpv_mpfr); mpfr_clear(tol_mpfr); mpfr_clear(abs_step);
        throw std::runtime_error("nsolve: failed to converge");
    }

    Expr result = make<Float>(static_cast<mpfr_srcptr>(x_mpfr), dps);
    mpfr_clear(x_mpfr); mpfr_clear(step_mpfr); mpfr_clear(fv_mpfr);
    mpfr_clear(fpv_mpfr); mpfr_clear(tol_mpfr); mpfr_clear(abs_step);
    return result;
}

// ----- solve_univariate_inequality -----------------------------------------

namespace {

[[nodiscard]] std::optional<double> sample_at(const Expr& f, const Expr& var,
                                                double x_val, int dps) {
    Expr xe = make<Float>(std::to_string(x_val), dps);
    Expr v = evalf(subs(f, var, xe), dps);
    return expr_to_double(v, dps);
}

}  // namespace

SetPtr solve_univariate_inequality(const Expr& lhs, const Expr& rel_op_str,
                                    const Expr& rhs, const Expr& var) {
    Expr f = lhs - rhs;
    auto roots = solve(f, var);
    // Keep each real root as its EXACT expression, paired with a numeric value
    // for ordering / sampling — so the interval endpoints are exact (−2, not
    // −2.0000…) and the unbounded ends are the real ±∞, not a 1e30 proxy.
    std::vector<std::pair<Expr, double>> rv;
    for (const auto& r : roots) {
        if (auto d = expr_to_double(r, 30)) rv.push_back({r, *d});
    }
    std::sort(rv.begin(), rv.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });
    rv.erase(std::unique(rv.begin(), rv.end(),
                         [](const auto& a, const auto& b) {
                             return a.second == b.second;
                         }),
             rv.end());

    std::string op = rel_op_str->str();
    bool want_pos = (op == ">" || op == ">=");
    bool want_neg = (op == "<" || op == "<=");
    bool include_eq = (op == ">=" || op == "<=");
    bool not_equal = (op == "!=");

    if (not_equal) {
        std::vector<Expr> root_exprs;
        for (const auto& [e, d] : rv) root_exprs.push_back(e);
        return root_exprs.empty()
                   ? reals()
                   : set_complement(reals(), finite_set(std::move(root_exprs)));
    }

    // Sample the sign of f on each open region between consecutive roots
    // (and the two unbounded ends), keeping the regions whose sign matches.
    const std::size_t n = rv.size();
    SetPtr result = empty_set();
    for (std::size_t i = 0; i <= n; ++i) {
        // Sample point strictly inside region i.
        double sample;
        if (n == 0) sample = 0.0;
        else if (i == 0) sample = rv.front().second - 1.0;
        else if (i == n) sample = rv.back().second + 1.0;
        else sample = (rv[i - 1].second + rv[i].second) / 2.0;

        auto v = sample_at(f, var, sample, 30);
        if (!v) continue;
        if (!((want_pos && *v > 0) || (want_neg && *v < 0))) continue;

        Expr lo = (i == 0) ? S::NegativeInfinity() : rv[i - 1].first;
        Expr hi = (i == n) ? S::Infinity() : rv[i].first;
        bool lo_open = (i == 0) || !include_eq;
        bool hi_open = (i == n) || !include_eq;
        SetPtr piece = (i == 0 && i == n)
                           ? reals()
                           : interval(lo, hi, lo_open, hi_open);
        result = set_union(result, piece);
    }
    return result;
}

// ----- rsolve --------------------------------------------------------------

Expr rsolve(const std::vector<Expr>& coeffs, const Expr& n) {
    if (coeffs.size() < 2) {
        throw std::invalid_argument("rsolve: need at least order 1 recurrence");
    }
    auto x = symbol("__rsolve_x");
    Poly p(std::vector<Expr>(coeffs.begin(), coeffs.end()), x);
    auto roots = p.roots();
    if (roots.empty()) {
        throw std::runtime_error(
            "rsolve: characteristic polynomial has no closed-form roots");
    }
    // Build Σ C_i * r_i^n. Track multiplicities for repeated roots:
    //   mult-k root r contributes (C_0 + C_1 n + ... + C_{k-1} n^{k-1}) r^n.
    std::vector<std::pair<Expr, std::size_t>> distinct;
    for (const auto& r : roots) {
        bool merged = false;
        for (auto& [val, count] : distinct) {
            if (val == r) {
                ++count;
                merged = true;
                break;
            }
        }
        if (!merged) distinct.emplace_back(r, 1);
    }
    Expr result = S::Zero();
    int next_c = 0;
    for (const auto& [root, mult] : distinct) {
        Expr poly_coef = S::Zero();
        for (std::size_t k = 0; k < mult; ++k) {
            auto C = symbol("__C" + std::to_string(next_c++));
            poly_coef = poly_coef + C * pow(n, integer(static_cast<long>(k)));
        }
        result = result + poly_coef * pow(root, n);
    }
    return result;
}

// ----- nonlinsolve --------------------------------------------------------

std::vector<std::vector<Expr>> nonlinsolve(
    const std::vector<Expr>& eqs, const std::vector<Expr>& vars) {
    if (eqs.size() != 2 || vars.size() != 2) {
        throw std::invalid_argument(
            "nonlinsolve: minimal implementation handles 2 eqs in 2 vars");
    }
    const Expr& x = vars[0];
    const Expr& y = vars[1];
    // Eliminate y via resultant in y; the resultant is a polynomial in x.
    Expr res = resultant(eqs[0], eqs[1], y);
    Expr res_simp = simplify(res);
    auto x_vals = solve(res_simp, x);
    std::vector<std::vector<Expr>> out;
    for (auto& xv : x_vals) {
        // Substitute into eqs[0] and solve for y.
        Expr eq_y = simplify(subs(eqs[0], x, xv));
        auto y_vals = solve(eq_y, y);
        for (auto& yv : y_vals) {
            // Verify the candidate (xv, yv) actually satisfies eqs[1] —
            // resultant elimination introduces extraneous solutions when
            // either equation has higher degree. Symbolic simplify can't
            // always reduce nested radicals, so fall back to a high-
            // precision numeric check.
            Expr check = simplify(subs(subs(eqs[1], x, xv), y, yv));
            if (!(check == S::Zero())) {
                auto v = expr_to_double(check, 30);
                if (!v || std::fabs(*v) > 1e-9) continue;
            }
            // De-duplicate.
            bool dup = false;
            for (const auto& prev : out) {
                if (prev[0] == xv && prev[1] == yv) { dup = true; break; }
            }
            if (!dup) out.push_back({xv, yv});
        }
    }
    return out;
}

// ----- linear_diophantine --------------------------------------------------

namespace {

// Extended Euclidean: gcd(a, b) = a*x + b*y. Operates on mpz.
[[nodiscard]] mpz_class ext_gcd_mpz(const mpz_class& a, const mpz_class& b,
                                     mpz_class& x, mpz_class& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    mpz_class x1, y1;
    mpz_class q = a / b;
    mpz_class g = ext_gcd_mpz(b, a - q * b, x1, y1);
    x = y1;
    y = x1 - q * y1;
    return g;
}

[[nodiscard]] std::optional<mpz_class> as_mpz(const Expr& e) {
    if (e->type_id() == TypeId::Integer) {
        return static_cast<const Integer&>(*e).value();
    }
    return std::nullopt;
}

}  // namespace

// ----- reduce_inequalities --------------------------------------------------

namespace {

[[nodiscard]] SetPtr reduce_single_rel(const Expr& rel, const Expr& var) {
    if (rel->type_id() != TypeId::Relational) {
        // Pass through; treat as Reals (assume vacuously true).
        return reals();
    }
    const auto& r = static_cast<const Relational&>(*rel);
    std::string op_name;
    switch (r.kind()) {
        case RelKind::Lt: op_name = "<"; break;
        case RelKind::Le: op_name = "<="; break;
        case RelKind::Gt: op_name = ">"; break;
        case RelKind::Ge: op_name = ">="; break;
        case RelKind::Ne: op_name = "!="; break;
        case RelKind::Eq:
            // Equality reduces to FiniteSet of solutions.
            return solveset(r.lhs() - r.rhs(), var);
    }
    auto op_sym = symbol(op_name);
    return solve_univariate_inequality(r.lhs(), op_sym, r.rhs(), var);
}

}  // namespace

SetPtr reduce_inequalities(const Expr& rel, const Expr& var) {
    return reduce_single_rel(rel, var);
}

SetPtr reduce_inequalities(const std::vector<Expr>& rels, const Expr& var,
                             bool conjunction) {
    if (rels.empty()) return conjunction ? reals() : empty_set();
    SetPtr acc = reduce_single_rel(rels[0], var);
    for (std::size_t i = 1; i < rels.size(); ++i) {
        SetPtr next = reduce_single_rel(rels[i], var);
        acc = conjunction ? set_intersection(acc, next)
                          : set_union(acc, next);
    }
    return acc;
}

// ----- Pythagorean Diophantine ---------------------------------------------

std::vector<std::array<Expr, 3>> pythagorean_triples(long max_z) {
    std::vector<std::array<Expr, 3>> out;
    if (max_z < 5) return out;  // smallest is (3, 4, 5)
    for (long m = 2; m * m <= max_z; ++m) {
        for (long n = 1; n < m; ++n) {
            // gcd(m, n) == 1 and opposite parity → primitive.
            mpz_class mz(m), nz(n);
            mpz_class g;
            mpz_gcd(g.get_mpz_t(), mz.get_mpz_t(), nz.get_mpz_t());
            if (g != 1) continue;
            if ((m % 2) == (n % 2)) continue;
            long x_val = m * m - n * n;
            long y_val = 2 * m * n;
            long z_val = m * m + n * n;
            if (z_val > max_z) break;
            out.push_back({integer(x_val), integer(y_val), integer(z_val)});
        }
    }
    return out;
}

// ----- Gröbner basis (Buchberger) and multivariate nonlinsolve --------------
//
// Minimal Buchberger over ℚ with lex order. Polynomials are represented
// as vectors of (coefficient, exponent_tuple) terms where the tuple is
// indexed by variable position. We build dense vectors keyed on
// monomials.

namespace {

// Multivariate monomial: vector of long exponents, one per variable.
struct Monomial {
    std::vector<long> e;
    // Lex ordering: compare exponents left-to-right. Returns -1/0/+1.
    [[nodiscard]] int cmp(const Monomial& o) const {
        for (std::size_t i = 0; i < e.size(); ++i) {
            if (e[i] != o.e[i]) return e[i] < o.e[i] ? -1 : 1;
        }
        return 0;
    }
    [[nodiscard]] bool operator==(const Monomial& o) const { return cmp(o) == 0; }
    [[nodiscard]] bool operator<(const Monomial& o) const { return cmp(o) < 0; }
};

// Multivariate polynomial as map monomial → mpq coefficient.
struct MPoly {
    std::vector<std::pair<Monomial, mpq_class>> terms;  // sorted by monomial desc

    [[nodiscard]] bool is_zero() const { return terms.empty(); }
    [[nodiscard]] const Monomial& leading_mono() const { return terms.front().first; }
    [[nodiscard]] const mpq_class& leading_coef() const { return terms.front().second; }

    void canonicalize() {
        // Sort descending by monomial.
        std::sort(terms.begin(), terms.end(),
                   [](const auto& a, const auto& b) { return b.first < a.first; });
        // Combine equal monomials.
        std::vector<std::pair<Monomial, mpq_class>> out;
        for (auto& t : terms) {
            if (!out.empty() && out.back().first == t.first) {
                out.back().second += t.second;
                if (out.back().second == 0) out.pop_back();
            } else if (!(t.second == 0)) {
                out.push_back(std::move(t));
            }
        }
        terms = std::move(out);
    }
};

// Convert an Expr polynomial in `vars` to an MPoly via expand + walk.
[[nodiscard]] MPoly expr_to_mpoly(const Expr& expr,
                                    const std::vector<Expr>& vars) {
    MPoly p;
    Expr e = expand(expr);
    auto add_term = [&](const mpq_class& coef,
                        std::vector<long> exps) {
        p.terms.push_back({Monomial{std::move(exps)}, coef});
    };
    auto handle_term = [&](const Expr& term) -> bool {
        // term must factor as (numeric coef) * (Π vars[i]^exp_i).
        std::vector<long> exps(vars.size(), 0);
        mpq_class coef = 1;
        auto handle_factor = [&](const Expr& f) -> bool {
            if (f->type_id() == TypeId::Integer) {
                coef *= mpq_class(static_cast<const Integer&>(*f).value());
                return true;
            }
            if (f->type_id() == TypeId::Rational) {
                coef *= static_cast<const Rational&>(*f).value();
                return true;
            }
            // Match var or var^k.
            for (std::size_t i = 0; i < vars.size(); ++i) {
                if (f == vars[i]) { exps[i] += 1; return true; }
                if (f->type_id() == TypeId::Pow && f->args()[0] == vars[i]
                    && f->args()[1]->type_id() == TypeId::Integer) {
                    long k = static_cast<const Integer&>(*f->args()[1]).to_long();
                    if (k < 0) return false;
                    exps[i] += k;
                    return true;
                }
            }
            return false;
        };
        if (term->type_id() == TypeId::Mul) {
            for (const auto& f : term->args()) {
                if (!handle_factor(f)) return false;
            }
        } else {
            if (!handle_factor(term)) return false;
        }
        if (coef == 0) return true;
        add_term(coef, std::move(exps));
        return true;
    };
    if (e->type_id() == TypeId::Add) {
        for (const auto& t : e->args()) {
            if (!handle_term(t)) {
                p.terms.clear();  // malformed; signal failure below
                return p;
            }
        }
    } else {
        if (!handle_term(e)) p.terms.clear();
    }
    p.canonicalize();
    return p;
}

[[nodiscard]] Expr mpoly_to_expr(const MPoly& p,
                                   const std::vector<Expr>& vars) {
    if (p.terms.empty()) return S::Zero();
    std::vector<Expr> term_exprs;
    for (const auto& [mono, coef] : p.terms) {
        std::vector<Expr> factors;
        // coefficient
        mpq_class cc = coef;
        cc.canonicalize();
        if (cc.get_den() == 1) {
            factors.push_back(make<Integer>(cc.get_num()));
        } else {
            factors.push_back(make<Rational>(cc));
        }
        for (std::size_t i = 0; i < vars.size(); ++i) {
            if (mono.e[i] == 0) continue;
            if (mono.e[i] == 1) factors.push_back(vars[i]);
            else factors.push_back(pow(vars[i], integer(mono.e[i])));
        }
        term_exprs.push_back(mul(factors));
    }
    return add(term_exprs);
}

// MPoly arithmetic helpers.
[[nodiscard]] MPoly mpoly_sub(const MPoly& a, const MPoly& b) {
    MPoly out = a;
    for (const auto& t : b.terms) {
        out.terms.push_back({t.first, -t.second});
    }
    out.canonicalize();
    return out;
}

[[nodiscard]] MPoly mpoly_scale_mono(const MPoly& a, const mpq_class& coef,
                                        const Monomial& mono) {
    MPoly out;
    out.terms.reserve(a.terms.size());
    for (const auto& t : a.terms) {
        Monomial m;
        m.e.resize(mono.e.size());
        for (std::size_t i = 0; i < mono.e.size(); ++i) {
            m.e[i] = t.first.e[i] + mono.e[i];
        }
        out.terms.push_back({std::move(m), t.second * coef});
    }
    return out;  // already sorted (multiplication preserves order)
}

[[nodiscard]] bool mono_divides(const Monomial& a, const Monomial& b) {
    for (std::size_t i = 0; i < a.e.size(); ++i) {
        if (a.e[i] > b.e[i]) return false;
    }
    return true;
}

[[nodiscard]] Monomial mono_sub(const Monomial& a, const Monomial& b) {
    Monomial m;
    m.e.resize(a.e.size());
    for (std::size_t i = 0; i < a.e.size(); ++i) m.e[i] = a.e[i] - b.e[i];
    return m;
}

[[nodiscard]] Monomial mono_lcm(const Monomial& a, const Monomial& b) {
    Monomial m;
    m.e.resize(a.e.size());
    for (std::size_t i = 0; i < a.e.size(); ++i) {
        m.e[i] = std::max(a.e[i], b.e[i]);
    }
    return m;
}

// Multivariate polynomial reduction: reduce f modulo a list G of
// polynomials. Returns the remainder.
[[nodiscard]] MPoly mpoly_reduce(MPoly f, const std::vector<MPoly>& G) {
    MPoly r;
    while (!f.is_zero()) {
        bool reduced = false;
        for (const auto& g : G) {
            if (g.is_zero()) continue;
            if (mono_divides(g.leading_mono(), f.leading_mono())) {
                Monomial diff = mono_sub(f.leading_mono(), g.leading_mono());
                mpq_class coef = f.leading_coef() / g.leading_coef();
                MPoly scaled = mpoly_scale_mono(g, coef, diff);
                f = mpoly_sub(f, scaled);
                reduced = true;
                break;
            }
        }
        if (!reduced) {
            r.terms.push_back(f.terms.front());
            f.terms.erase(f.terms.begin());
        }
    }
    return r;
}

// S-polynomial of f, g.
[[nodiscard]] MPoly s_poly(const MPoly& f, const MPoly& g) {
    Monomial l = mono_lcm(f.leading_mono(), g.leading_mono());
    Monomial diff_f = mono_sub(l, f.leading_mono());
    Monomial diff_g = mono_sub(l, g.leading_mono());
    MPoly fterm = mpoly_scale_mono(f, mpq_class(1) / f.leading_coef(), diff_f);
    MPoly gterm = mpoly_scale_mono(g, mpq_class(1) / g.leading_coef(), diff_g);
    return mpoly_sub(fterm, gterm);
}

}  // namespace

std::vector<Expr> groebner(const std::vector<Expr>& eqs,
                             const std::vector<Expr>& vars) {
    std::vector<MPoly> G;
    G.reserve(eqs.size());
    for (const auto& e : eqs) {
        auto p = expr_to_mpoly(e, vars);
        if (!p.is_zero()) G.push_back(std::move(p));
    }
    // Buchberger's algorithm.
    std::size_t i = 0;
    while (i < G.size()) {
        for (std::size_t j = 0; j < i; ++j) {
            MPoly s = s_poly(G[j], G[i]);
            MPoly r = mpoly_reduce(s, G);
            if (!r.is_zero()) G.push_back(std::move(r));
        }
        ++i;
        if (G.size() > 200) break;  // safety bound
    }
    // Reduce basis: drop polynomials whose leading mono is divisible
    // by another's, normalize each to monic.
    std::vector<MPoly> reduced;
    for (std::size_t k = 0; k < G.size(); ++k) {
        bool drop = false;
        for (std::size_t l = 0; l < G.size(); ++l) {
            if (l == k || G[l].is_zero()) continue;
            if (mono_divides(G[l].leading_mono(), G[k].leading_mono())
                && !(G[l].leading_mono() == G[k].leading_mono())) {
                drop = true;
                break;
            }
        }
        if (!drop && !G[k].is_zero()) {
            MPoly p = G[k];
            mpq_class lc = p.leading_coef();
            for (auto& t : p.terms) t.second /= lc;
            reduced.push_back(std::move(p));
        }
    }
    // Final pass: reduce each polynomial against the others.
    for (std::size_t k = 0; k < reduced.size(); ++k) {
        std::vector<MPoly> others;
        for (std::size_t l = 0; l < reduced.size(); ++l) {
            if (l != k) others.push_back(reduced[l]);
        }
        MPoly red = mpoly_reduce(reduced[k], others);
        // Re-normalize.
        if (!red.is_zero()) {
            mpq_class lc = red.leading_coef();
            for (auto& t : red.terms) t.second /= lc;
        }
        reduced[k] = std::move(red);
    }
    std::vector<Expr> out;
    out.reserve(reduced.size());
    for (auto& p : reduced) {
        if (!p.is_zero()) out.push_back(mpoly_to_expr(p, vars));
    }
    return out;
}

std::vector<std::vector<Expr>> nonlinsolve_groebner(
    const std::vector<Expr>& eqs, const std::vector<Expr>& vars) {
    auto basis = groebner(eqs, vars);
    if (basis.empty()) return {};
    // Lex order Gröbner basis with vars = [x, y, z, ...] is triangular
    // with the last variable's univariate polynomial *last* in the
    // elimination ideal. Solve in reverse order: z first, then y given
    // z, then x given y, z. Build solutions as a vector indexed by the
    // original `vars` order.
    const std::size_t n = vars.size();
    std::function<std::vector<std::vector<Expr>>(std::vector<Expr>, long)>
        rec = [&](std::vector<Expr> partial,
                  long var_idx) -> std::vector<std::vector<Expr>> {
        if (var_idx < 0) {
            return {partial};
        }
        // Substitute already-known values (for indices > var_idx) into
        // each basis element; pick the first that becomes univariate
        // in vars[var_idx].
        for (const auto& g : basis) {
            Expr g_sub = g;
            for (std::size_t k = static_cast<std::size_t>(var_idx + 1);
                 k < n; ++k) {
                g_sub = subs(g_sub, vars[k], partial[k]);
            }
            g_sub = expand(simplify(g_sub));
            // Must contain vars[var_idx] and no earlier variable
            // vars[k] for k < var_idx.
            if (!has(g_sub, vars[static_cast<std::size_t>(var_idx)])) continue;
            bool ok = true;
            for (long k = 0; k < var_idx; ++k) {
                if (has(g_sub, vars[static_cast<std::size_t>(k)])) {
                    ok = false; break;
                }
            }
            if (!ok) continue;
            auto roots = solve(g_sub,
                                vars[static_cast<std::size_t>(var_idx)]);
            std::vector<std::vector<Expr>> out;
            for (auto& r : roots) {
                auto next = partial;
                next[static_cast<std::size_t>(var_idx)] = r;
                auto sub_sols = rec(next, var_idx - 1);
                for (auto& s : sub_sols) out.push_back(std::move(s));
            }
            return out;
        }
        return {};  // no triangular polynomial found
    };
    std::vector<Expr> partial(n, S::Zero());
    return rec(partial, static_cast<long>(n) - 1);
}

// ----- linear_diophantine --------------------------------------------------

std::optional<std::pair<Expr, Expr>>
linear_diophantine(const Expr& a, const Expr& b, const Expr& c) {
    auto a_z = as_mpz(a);
    auto b_z = as_mpz(b);
    auto c_z = as_mpz(c);
    if (!a_z || !b_z || !c_z) return std::nullopt;
    mpz_class x, y;
    mpz_class g = ext_gcd_mpz(*a_z, *b_z, x, y);
    // Reduce sign: we want g > 0.
    if (sgn(g) < 0) { g = -g; x = -x; y = -y; }
    if ((*c_z) % g != 0) return std::nullopt;
    mpz_class scale = (*c_z) / g;
    x *= scale;
    y *= scale;
    return std::pair<Expr, Expr>{make<Integer>(x), make<Integer>(y)};
}

}  // namespace sympp
