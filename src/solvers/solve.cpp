#include <sympp/solvers/solve.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include <gmpxx.h>
#include <mpfr.h>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/add.hpp>
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
#include <sympp/polys/poly.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp {

std::vector<Expr> solve(const Expr& expr, const Expr& var) {
    Poly p(expr, var);
    return p.roots();
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

SetPtr solveset(const Expr& expr, const Expr& var, const SetPtr& domain) {
    auto roots = solve(expr, var);
    if (roots.empty()) return empty_set();
    // Filter against the domain when membership is decidable.
    std::vector<Expr> kept;
    kept.reserve(roots.size());
    for (auto& r : roots) {
        auto m = domain->contains(r);
        // Keep when membership is true OR unknown (be inclusive on
        // symbolic results).
        if (m == std::optional<bool>{false}) continue;
        kept.push_back(std::move(r));
    }
    return finite_set(std::move(kept));
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
    // Filter to numeric roots; sort.
    std::vector<double> numeric_roots;
    for (const auto& r : roots) {
        auto d = expr_to_double(r, 30);
        if (d) numeric_roots.push_back(*d);
    }
    std::sort(numeric_roots.begin(), numeric_roots.end());
    numeric_roots.erase(std::unique(numeric_roots.begin(), numeric_roots.end()),
                         numeric_roots.end());

    // Decode the relation: rel_op_str is a Symbol with one of these names.
    std::string op = rel_op_str->str();
    bool want_pos = (op == ">" || op == ">=");
    bool want_neg = (op == "<" || op == "<=");
    bool include_eq = (op == ">=" || op == "<=");
    bool not_equal = (op == "!=");

    // Build the disjoint sample intervals.
    const double huge = 1e30;
    std::vector<double> bps;
    bps.push_back(-huge);
    for (auto r : numeric_roots) bps.push_back(r);
    bps.push_back(huge);

    SetPtr result = empty_set();
    for (std::size_t i = 0; i + 1 < bps.size(); ++i) {
        double mid = (bps[i] + bps[i + 1]) / 2.0;
        auto v = sample_at(f, var, mid, 30);
        if (!v) continue;
        bool keep = false;
        if (not_equal) keep = true;
        else if (want_pos && *v > 0) keep = true;
        else if (want_neg && *v < 0) keep = true;
        if (keep) {
            // Translate bps[i], bps[i+1] back to Exprs. -huge / +huge
            // act as ±∞ proxies; we still emit numeric Floats since we
            // don't have an Infinity singleton.
            Expr lo = (i == 0) ? make<Float>(std::to_string(-huge), 30)
                                : make<Float>(std::to_string(bps[i]), 30);
            Expr hi = (i + 1 == bps.size() - 1)
                          ? make<Float>(std::to_string(huge), 30)
                          : make<Float>(std::to_string(bps[i + 1]), 30);
            // Endpoint inclusion: at internal boundaries the boundary
            // is a root → strict inequality excludes, ≤/≥ includes.
            // For ±∞ proxies we always exclude.
            bool lo_open = (i == 0) || !include_eq;
            bool hi_open = (i + 1 == bps.size() - 1) || !include_eq;
            result = set_union(result,
                                interval(lo, hi, lo_open, hi_open));
        }
    }
    // For !=, the "kept" intervals collectively form Reals \ {roots}.
    if (not_equal && result->kind() != SetKind::Empty) {
        // Collapse to Reals \ FiniteSet(roots).
        std::vector<Expr> root_exprs;
        for (auto r : numeric_roots) {
            root_exprs.push_back(make<Float>(std::to_string(r), 30));
        }
        return set_complement(reals(), finite_set(std::move(root_exprs)));
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
