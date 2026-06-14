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
#include <sympp/simplify/simplify.hpp>

namespace sympp {

namespace {

// Polynomial-only root finder. This is the original solve() body; it is
// kept separate so solveset()'s internal fallback can use it WITHOUT
// re-entering the transcendental solveset path below (which would recurse).
[[nodiscard]] std::vector<Expr> solve_poly(const Expr& expr, const Expr& var) {
    Poly p(expr, var);
    return p.roots();
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

// Append the representative roots of f(B*var) = c over one principal period
// (var = theta/B), deduplicating into `roots`. Mirrors SymPy's `solve`, which
// inverts the inner function and divides by B rather than enumerating every x
// in [0, 2pi). For sin/cos an out-of-range numeric c contributes nothing.
void append_trig_roots(FunctionId id, const Expr& c, const Expr& B,
                       const Expr& var, std::vector<Expr>& roots) {
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
        Expr r = simplify(th * pow(B, integer(-1)));
        if (has(r, var)) continue;
        if (std::none_of(roots.begin(), roots.end(),
                         [&](const Expr& u) { return u == r; })) {
            roots.push_back(std::move(r));
        }
    }
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
    // arg must be B*var with B var-free and nonzero (no additive phase).
    const Expr& arg = fexpr->args()[0];
    Expr B = simplify(arg * pow(var, integer(-1)));
    if (has(B, var) || simplify(B * var - arg) != S::Zero()) return std::nullopt;
    // f(arg) = c, c = -C/A.
    Expr c = simplify(mul({S::NegativeOne(), cst, pow(coeff, integer(-1))}));
    std::vector<Expr> roots;   // var = theta / B, deduplicated.
    append_trig_roots(id, c, B, var, roots);
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
    // arg must be B*var with B var-free and nonzero (no additive phase).
    const Expr& arg = g->args()[0];
    Expr B = simplify(arg * pow(var, integer(-1)));
    if (has(B, var) || simplify(B * var - arg) != S::Zero()) return std::nullopt;
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
        append_trig_roots(id, c, B, var, roots);
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
    // c is in range when known-numeric and within bounds; an unbounded range
    // accepts anything (including a symbolic c).
    auto in_range = [&](const Expr& c, double lo, double hi, bool bounded) {
        if (!bounded) return true;
        auto v = numeric_value(c);
        if (!v) return false;  // symbolic c with a bounded range — conservative
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

// Solve the canonical Lambert-W form a·var·exp(b·var) + c = 0 → var =
// W(−c·b/a)/b, using the defining identity W(z)·e^(W(z)) = z. Detects a single
// term a·var·exp(b·var) (a, b var-free, b ≠ 0, var to the first power) plus a
// var-free constant. Matches SymPy's LambertW result.
[[nodiscard]] std::optional<std::vector<Expr>>
solve_lambert(const Expr& expr, const Expr& var) {
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
    // dep = a · var · exp(b·var): a var-free coeff, one bare var, one exp(b·var).
    Expr a = S::One();
    Expr b;
    bool have_var = false, have_exp = false;
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
        if (f->type_id() == TypeId::Function
            && static_cast<const Function&>(*f).function_id()
                   == FunctionId::Exp) {
            if (have_exp) return std::nullopt;
            const Expr& arg = f->args()[0];
            Expr bb = simplify(arg * pow(var, integer(-1)));
            if (has(bb, var) || simplify(bb * var - arg) != S::Zero()) {
                return std::nullopt;  // exp argument is not linear b·var
            }
            b = bb;
            have_exp = true;
            continue;
        }
        return std::nullopt;  // another var-dependent factor
    }
    if (!have_var || !have_exp || b == S::Zero()) return std::nullopt;
    // a·var·e^(b·var) = −c → (b·var)·e^(b·var) = −c·b/a → var = W(−c·b/a)/b.
    Expr z = simplify(mul({S::NegativeOne(), cst, b, pow(a, integer(-1))}));
    return std::vector<Expr>{simplify(lambertw(z) * pow(b, integer(-1)))};
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

}  // namespace

std::vector<Expr> solve(const Expr& expr, const Expr& var) {
    // Expand first so a factored polynomial reaches the Poly machinery:
    // solve((x-1)*(x-2)) was empty because Poly couldn't build from the Mul.
    auto roots = solve_poly(expand(expr), var);
    // A genuine solution x = c must be free of x. solve_poly treats a
    // var-dependent "coefficient" (e.g. exp(x) in x·exp(x) − 1) as constant and
    // can hand back a rearrangement such as x = exp(x)**(-1); discard any
    // candidate that still contains var rather than return a non-solution.
    std::erase_if(roots, [&](const Expr& r) { return has(r, var); });
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
    // Rational equation (1/x + … = 0): clear the denominator and solve the
    // numerator, dropping pole roots. The polynomial path above can't build a
    // Poly from a negative-power term, so these reach here empty.
    if (auto rr = solve_rational(expr, var); rr && !rr->empty()) return *rr;
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
        if (auto inv = solve_inverse_func_poly(expr, var); inv && !inv->empty())
            return *inv;
        if (auto lw = solve_lambert(expr, var); lw && !lw->empty()) return *lw;
        if (auto tl = solve_trig_linear(expr, var); tl && !tl->empty())
            return *tl;
        if (auto tr = solve_trig(expr, var); tr && !tr->empty()) return *tr;
        SetPtr s = solveset(expr, var);
        if (s && s->kind() == SetKind::FiniteSet) {
            auto elems = static_cast<const FiniteSet&>(*s).elements();
            std::erase_if(elems, [&](const Expr& r) { return has(r, var); });
            return elems;
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
