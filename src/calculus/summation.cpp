#include <sympp/calculus/summation.hpp>

#include <optional>
#include <utility>
#include <vector>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/add.hpp>
#include <sympp/core/expand.hpp>
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
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/core/undefined_function.hpp>
#include <sympp/functions/combinatorial.hpp>
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/functions/miscellaneous.hpp>
#include <sympp/functions/special.hpp>
#include <sympp/functions/trigonometric.hpp>
#include <sympp/polys/poly.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp {

namespace {

// Unevaluated sum marker — an UndefinedFunction `Sum(expr, var, lo, hi)`,
// mirroring the `Integral(_, _)` marker. Returned when no closed form is
// found, so a sum is never silently dropped to its bare summand.
[[nodiscard]] Expr sum_marker(const Expr& expr, const Expr& var,
                              const Expr& lo, const Expr& hi) {
    return function_symbol("Sum")(std::vector<Expr>{expr, var, lo, hi});
}

// Telescoping sum of a rational summand c / D(var) with a var-free numerator c,
// where D factors into distinct linear factors whose rational roots pairwise
// differ by integers (degree ≥ 2). Partial fractions + telescoping give a closed
// form; see the body for the derivation. Covers 1/(k(k+1)), 1/(4k²−1),
// 1/(k(k+1)(k+2)), … Returns nullopt unless the summand has this shape. The closed
// form is exact (equivalent to SymPy's, though it may be presented differently).
[[nodiscard]] std::optional<Expr> telescope_rational(
    const Expr& expr, const Expr& var, const Expr& lo, const Expr& hi) {
    // Split expr = c · D^{-1}: c var-free, D one or more reciprocal factors.
    Expr c = S::One();
    Expr denom;
    auto take_recip = [&](const Expr& f) -> bool {
        if (f->type_id() == TypeId::Pow && f->args()[1] == S::NegativeOne()) {
            denom = denom ? mul({denom, f->args()[0]}) : f->args()[0];
            return true;
        }
        return false;
    };
    if (expr->type_id() == TypeId::Pow) {
        if (!take_recip(expr)) return std::nullopt;
    } else if (expr->type_id() == TypeId::Mul) {
        std::vector<Expr> cf;
        for (const auto& f : expr->args()) {
            if (take_recip(f)) continue;
            if (has(f, var)) return std::nullopt;  // var in the numerator
            cf.push_back(f);
        }
        if (!cf.empty()) c = mul(std::move(cf));
    } else {
        return std::nullopt;
    }
    if (!denom || !has(denom, var)) return std::nullopt;
    // D must factor into distinct linear factors with rational roots that pairwise
    // differ by integers (degree ≥ 2). Partial fractions give c/D = Σ Aᵢ/(k−rᵢ);
    // with the largest root as a reference, each 1/(k−rᵢ) = u(k+mᵢ) (mᵢ = ref−rᵢ ≥ 0,
    // u(k)=1/(k−ref)), so the summand is Σ Aᵢ(u(k+mᵢ) − u(k)) (the −u(k) terms cancel
    // since Σ Aᵢ = 0 for a constant numerator over degree ≥ 2). Each piece
    // telescopes to Σ Aᵢ[ Σ_{j=1}^{mᵢ} u(hi+j) − Σ_{j=0}^{mᵢ−1} u(lo+j) ].
    Poly p(expand(denom), var);
    const std::size_t deg = p.degree();
    if (deg < 2) return std::nullopt;
    std::vector<Expr> rts = p.roots();
    if (rts.size() != deg) return std::nullopt;  // need every root in closed form
    auto is_rat = [](const Expr& e) {
        return e->type_id() == TypeId::Integer
               || e->type_id() == TypeId::Rational;
    };
    for (const auto& r : rts) {
        if (!is_rat(r)) return std::nullopt;
    }
    // Pairwise: distinct, and differences must be (long-representable) integers.
    for (std::size_t i = 0; i < rts.size(); ++i) {
        for (std::size_t j = i + 1; j < rts.size(); ++j) {
            Expr diff = simplify(rts[i] - rts[j]);
            if (diff == S::Zero()) return std::nullopt;        // repeated root
            if (diff->type_id() != TypeId::Integer
                || !static_cast<const Integer&>(*diff).fits_long()) {
                return std::nullopt;  // non-integer gap → digamma, not rational
            }
        }
    }
    const Expr lead = p.leading_coeff();
    // Reference = the maximum root, so every mᵢ = ref − rᵢ ≥ 0.
    Expr ref = rts[0];
    for (const auto& r : rts) {
        if (is_positive(simplify(r - ref)) == std::optional<bool>{true}) ref = r;
    }
    // Guard against a pole inside the range: an integer root ≥ lo makes a summand
    // undefined. A non-integer root is never hit by an integer index.
    auto safe_root = [&](const Expr& r) {
        if (r->type_id() != TypeId::Integer) return true;
        return is_positive(simplify(lo - r)) == std::optional<bool>{true};
    };
    for (const auto& r : rts) {
        if (!safe_root(r)) return std::nullopt;
    }
    auto u = [&](const Expr& k) { return pow(k - ref, integer(-1)); };
    const bool infinite = hi->type_id() == TypeId::Infinity;
    Expr total = S::Zero();
    for (std::size_t i = 0; i < rts.size(); ++i) {
        // Aᵢ = c / (lead · ∏_{j≠i}(rᵢ − rⱼ)).
        Expr prod = lead;
        for (std::size_t j = 0; j < rts.size(); ++j) {
            if (j != i) prod = mul(prod, simplify(rts[i] - rts[j]));
        }
        Expr a_i = mul(c, pow(prod, integer(-1)));
        Expr m_expr = simplify(ref - rts[i]);  // ≥ 0 integer
        const long m = static_cast<const Integer&>(*m_expr).to_long();
        Expr head = S::Zero();
        if (!infinite) {
            for (long j = 1; j <= m; ++j) head = add({head, u(hi + integer(j))});
        }
        Expr tail = S::Zero();
        for (long j = 0; j < m; ++j) tail = add({tail, u(lo + integer(j))});
        total = add({total,
                     mul({a_i, add({head, mul({S::NegativeOne(), tail})})})});
    }
    return simplify(total);
}

// Σ_{k=lo}^{hi} P(k)·r^k for a polynomial P (degree ≥ 1) and a concrete numeric
// ratio r ≠ 1. Finds the antidifference S(k) = Q(k)·r^k, where Q is a polynomial
// of the same degree satisfying r·Q(k+1) − Q(k) = P(k); the coefficients solve
// top-down since q_j = [P_j − r·Σ_{i>j} C(i,j) q_i]/(r−1). Then the sum is
// S(hi+1) − S(lo). Returns nullopt unless the summand has this shape.
[[nodiscard]] std::optional<Expr> sum_poly_geometric(
    const Expr& expr, const Expr& var, const Expr& lo, const Expr& hi) {
    if (expr->type_id() != TypeId::Mul) return std::nullopt;
    // Split off one geometric factor base^(c·var+d); the rest is the polynomial.
    Expr geo;
    std::vector<Expr> poly_factors;
    for (const auto& f : expr->args()) {
        if (!geo && f->type_id() == TypeId::Pow && !has(f->args()[0], var)
            && has(f->args()[1], var)) {
            geo = f;
        } else {
            poly_factors.push_back(f);
        }
    }
    if (!geo || poly_factors.empty()) return std::nullopt;
    Expr c = diff(geo->args()[1], var);
    if (has(c, var)) return std::nullopt;  // exponent not linear in var
    Expr d = simplify(geo->args()[1] - c * var);
    if (has(d, var)) return std::nullopt;
    Expr ratio = pow(geo->args()[0], c);
    // The ratio (and base) may be symbolic — e.g. Σ k·xᵏ = x(1−(n+1)xⁿ+n·xⁿ⁺¹)/(x−1)²,
    // the generating-function identity. A finite sum yields a clean closed form; an
    // infinite sum with a symbolic ratio fails the |r| < 1 convergence test below and
    // is left unevaluated (rather than emitting x**oo terms). r = 1 has no closed
    // form here (the poly·1ᵏ case is a power sum handled elsewhere).
    if (ratio == S::One()) return std::nullopt;
    Expr prefactor = pow(geo->args()[0], d);
    // Coefficients of P (P_i = coeff of var^i); must be a clean polynomial.
    std::vector<Expr> pc;
    try {
        Poly pp(expand(mul(poly_factors)), var);
        for (const auto& cc : pp.coeffs()) {
            if (has(cc, var)) return std::nullopt;
            pc.push_back(cc);
        }
    } catch (...) {
        return std::nullopt;
    }
    const long p = static_cast<long>(pc.size()) - 1;
    if (p < 1) return std::nullopt;  // degree 0 is plain geometric
    const Expr inv_rm1 = pow(simplify(ratio - integer(1)), integer(-1));
    std::vector<Expr> q(static_cast<std::size_t>(p) + 1);
    for (long j = p; j >= 0; --j) {
        Expr s = pc[static_cast<std::size_t>(j)];
        for (long i = j + 1; i <= p; ++i) {
            s = s - ratio * binomial(integer(i), integer(j))
                        * q[static_cast<std::size_t>(i)];
        }
        q[static_cast<std::size_t>(j)] = simplify(s * inv_rm1);
    }
    auto eval_q = [&](const Expr& m) -> Expr {
        std::vector<Expr> terms;
        terms.reserve(static_cast<std::size_t>(p) + 1);
        for (long i = 0; i <= p; ++i) {
            terms.push_back(q[static_cast<std::size_t>(i)]
                            * pow(m, integer(i)));
        }
        return add(std::move(terms));
    };
    Expr s_lo = eval_q(lo) * pow(ratio, lo);
    Expr s_hi;
    if (hi->type_id() == TypeId::Infinity) {
        // Σ_{k=lo}^∞ P(k)·r^k converges iff |r| < 1 (⇔ r² < 1); then the upper
        // boundary term Q(k)·r^k → 0 (geometric decay dominates the polynomial Q),
        // so the sum is −S(lo). A divergent or undecidable ratio is left alone.
        if (!(is_negative(simplify(ratio * ratio - S::One()))
              == std::optional<bool>{true})) {
            return std::nullopt;
        }
        s_hi = S::Zero();
    } else {
        s_hi = eval_q(hi + integer(1)) * pow(ratio, hi + integer(1));
    }
    return simplify(prefactor * (s_hi - s_lo));
}

// Σ_{k=lo}^∞ c·r^k / k! = c·e^r, minus the omitted head Σ_{k=0}^{lo-1} when
// lo > 0. The exponential series converges for every r, so no convergence test
// is needed. Recognizes c · P(var) · (∏ baseᵢ^(aᵢ·var + bᵢ)) · factorial(var)^(-1):
// the bases give the rate r = ∏ baseᵢ^{aᵢ} (and a constant from the bᵢ), and the
// polynomial numerator P(var) is folded in through the falling-factorial basis,
// since Σ_{k≥0} k^{(m)} rᵏ/k! = rᵐ e^r. Writing P = Σ_m c_m·k^{(m)} then gives
// Σ P(k) rᵏ/k! = (Σ_m c_m rᵐ) e^r = Q(r) e^r. So Σ k/k! = e, Σ k²/k! = 2e, etc.
//
// Falling-factorial transform: returns Q(r) for a polynomial whose power-basis
// coefficients are `coeffs` (ascending). Extracts the monic falling factorials
// k^{(m)} = k(k−1)…(k−m+1) top-down (a triangular solve, no Stirling numbers).
[[nodiscard]] Expr exp_series_poly_transform(std::vector<Expr> coeffs,
                                             const Expr& var, const Expr& rate) {
    std::vector<Expr> q_terms;
    for (long m = static_cast<long>(coeffs.size()) - 1; m >= 0; --m) {
        Expr cm = simplify(coeffs[static_cast<std::size_t>(m)]);
        if (cm == S::Zero()) continue;
        q_terms.push_back(mul(cm, pow(rate, integer(m))));
        if (m >= 1) {
            // Subtract cm · ∏_{i=0}^{m−1}(var − i) from the remaining coeffs.
            std::vector<Expr> ff_factors;
            for (long i = 0; i < m; ++i) {
                ff_factors.push_back(add(var, integer(-i)));
            }
            Poly ff(expand(mul(std::move(ff_factors))), var);
            const auto& fc = ff.coeffs();
            for (std::size_t j = 0; j < fc.size(); ++j) {
                coeffs[j] = simplify(
                    add(coeffs[j], mul(S::NegativeOne(), mul(cm, fc[j]))));
            }
        }
    }
    return q_terms.empty() ? Expr{S::Zero()} : add(std::move(q_terms));
}

[[nodiscard]] std::optional<Expr> sum_exponential_series(
    const Expr& expr, const Expr& var, const Expr& lo, const Expr& hi) {
    if (hi->type_id() != TypeId::Infinity) return std::nullopt;
    if (lo->type_id() != TypeId::Integer) return std::nullopt;
    const auto& loz = static_cast<const Integer&>(*lo);
    if (loz.is_negative() || !loz.fits_long()) return std::nullopt;

    std::vector<Expr> factors;
    if (expr->type_id() == TypeId::Mul) {
        for (const auto& f : expr->args()) factors.push_back(f);
    } else {
        factors.push_back(expr);
    }

    // A shifted factorial (k+m)! re-indexes to k! via j = k+m: the sum
    // Σ_{k=lo}^∞ P(k)/(k+m)! equals Σ_{j=lo+m}^∞ P(j−m)/j!, which is the m=0 case
    // below. Reduce to it by substituting var → var−m and shifting the lower bound.
    // (factorial(2k) and the like — a non-unit var coefficient — are left alone.)
    for (const auto& f : factors) {
        if (f->type_id() == TypeId::Pow && f->args()[1] == S::NegativeOne()
            && f->args()[0]->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*f->args()[0]);
            if (fn.function_id() == FunctionId::Factorial
                && fn.args().size() == 1) {
                Expr mm = simplify(fn.args()[0] + mul(S::NegativeOne(), var));
                if (!has(mm, var) && mm->type_id() == TypeId::Integer
                    && static_cast<const Integer&>(*mm).fits_long()
                    && static_cast<const Integer&>(*mm).to_long() != 0) {
                    const long m = static_cast<const Integer&>(*mm).to_long();
                    if (loz.to_long() + m < 0) return std::nullopt;  // factorial pole
                    Expr shifted =
                        subs(expr, var, var + mul(S::NegativeOne(), integer(m)));
                    return sum_exponential_series(shifted, var,
                                                  lo + integer(m), hi);
                }
                break;  // the bare factorial(var) case proceeds below
            }
        }
    }
    bool has_factorial = false;
    Expr cst = S::One();                // var-free constant prefactor c
    Expr rate = S::One();               // r
    std::vector<Expr> poly_factors;     // var-dependent numerator P(var)
    for (const auto& f : factors) {
        // 1/var! = factorial(var)^(-1).
        if (f->type_id() == TypeId::Pow && f->args()[1] == S::NegativeOne()
            && f->args()[0]->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*f->args()[0]);
            if (fn.function_id() == FunctionId::Factorial
                && fn.args().size() == 1 && fn.args()[0] == var) {
                if (has_factorial) return std::nullopt;  // (k!)^(-2): not this
                has_factorial = true;
                continue;
            }
        }
        if (!has(f, var)) {
            cst = mul(cst, f);
            continue;
        }
        // base^(a·var + b) with base var-free: base^a → rate, base^b → const.
        if (f->type_id() == TypeId::Pow && !has(f->args()[0], var)) {
            const Expr& base = f->args()[0];
            try {
                Poly p(expand(f->args()[1]), var);
                if (p.degree() > 1) return std::nullopt;
                const auto& cf = p.coeffs();
                Expr a = cf.size() > 1 ? cf[1] : Expr{S::Zero()};
                Expr b = !cf.empty() ? cf[0] : Expr{S::Zero()};
                if (has(a, var) || has(b, var)) return std::nullopt;
                rate = mul(rate, pow(base, a));
                cst = mul(cst, pow(base, b));
                continue;
            } catch (const std::exception&) {
                return std::nullopt;
            }
        }
        // Any other var-dependent factor must be polynomial in var (the
        // numerator P): k, (k+1), k², … are folded via the falling-factorial
        // transform; a non-polynomial factor (sin(k), 1/k, …) bails.
        poly_factors.push_back(f);
    }
    if (!has_factorial) return std::nullopt;
    rate = simplify(rate);
    cst = simplify(cst);

    // Build P(var) and its closed form Q(r): Σ P(k) r^k/k! = Q(r) e^r.
    Expr poly = poly_factors.empty() ? Expr{S::One()}
                                     : mul(std::move(poly_factors));
    std::vector<Expr> pcoeffs;
    try {
        Poly p(expand(poly), var);
        pcoeffs = p.coeffs();
    } catch (const std::exception&) {
        return std::nullopt;  // numerator is not a polynomial in var
    }
    // Poly() treats a non-polynomial factor (cos(k), √k, …) as an opaque degree-0
    // coefficient, which would leak the summation variable into the result
    // (Σ cos(k)/k! → e·cos(k)). Reject any var-dependent coefficient.
    for (const auto& cf : pcoeffs) {
        if (has(cf, var)) return std::nullopt;
    }
    Expr q_of_r = exp_series_poly_transform(pcoeffs, var, rate);
    Expr e_r = exp(rate);
    Expr tail = mul(q_of_r, e_r);  // Σ_{k≥0} P(k) r^k/k!

    const long L = loz.to_long();
    if (L == 0) return simplify(mul(cst, tail));
    // Subtract the omitted head terms Σ_{k=0}^{L-1} P(k) r^k/k!.
    std::vector<Expr> head;
    for (long kk = 0; kk < L; ++kk) {
        head.push_back(mul({subs(poly, var, integer(kk)),
                            pow(rate, integer(kk)),
                            pow(factorial(integer(kk)), integer(-1))}));
    }
    return simplify(
        mul(cst, add(tail, mul(S::NegativeOne(), add(std::move(head))))));
}

}  // namespace

// The binomial theorem Σ_{k=0}^n C(n,k)·rᵏ = (1+r)ⁿ. Matches a summand
// binomial(n,var)·base^(a·var+b) (the geometric factor optional) over k = 0…n,
// where n is exactly the binomial's first argument, and returns
// const·base^b·(1 + base^a)ⁿ. Covers Σ C(n,k) = 2ⁿ, Σ(−1)ᵏC(n,k) = 0,
// Σ 2ᵏC(n,k) = 3ⁿ, Σ xᵏC(n,k) = (1+x)ⁿ.
[[nodiscard]] std::optional<Expr> sum_binomial_theorem(const Expr& expr,
                                                       const Expr& var,
                                                       const Expr& lo,
                                                       const Expr& hi) {
    if (!(lo == S::Zero())) return std::nullopt;
    std::vector<Expr> factors;
    if (expr->type_id() == TypeId::Mul) {
        for (const auto& f : expr->args()) factors.push_back(f);
    } else {
        factors.push_back(expr);
    }
    bool have_binom = false;
    Expr ratio = S::One();       // ∏ base^a from the geometric factors
    Expr prefactor = S::One();   // ∏ base^b
    Expr coeff = S::One();
    for (const auto& f : factors) {
        if (!have_binom && f->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*f);
            if (fn.function_id() == FunctionId::Binomial
                && fn.args().size() == 2 && fn.args()[1] == var
                && fn.args()[0] == hi) {
                have_binom = true;
                continue;
            }
        }
        // base^(a·var+b): a constant base raised to a linear exponent in var.
        if (f->type_id() == TypeId::Pow && !has(f->args()[0], var)
            && has(f->args()[1], var)) {
            try {
                Poly pe(expand(f->args()[1]), var);
                if (pe.degree() != 1) return std::nullopt;
                const Expr& a = pe.coeffs()[1];
                const Expr& b = pe.coeffs()[0];
                if (has(a, var) || has(b, var)) return std::nullopt;
                ratio = mul(ratio, pow(f->args()[0], a));
                prefactor = mul(prefactor, pow(f->args()[0], b));
                continue;
            } catch (const std::exception&) {
                return std::nullopt;
            }
        }
        if (!has(f, var)) {
            coeff = mul(coeff, f);
            continue;
        }
        return std::nullopt;  // a var factor that isn't the binomial or geometric
    }
    if (!have_binom) return std::nullopt;
    Expr base_sum = simplify(integer(1) + ratio);
    // (1 + r)ⁿ, with (1−1)ⁿ = 0 for the alternating sum (n ≥ 1).
    Expr power = (base_sum == S::Zero()) ? Expr{S::Zero()} : pow(base_sum, hi);
    return simplify(mul(mul(coeff, prefactor), power));
}

// Σ_{k=lo}^{hi} c·P(k)/(k+m)! for a polynomial P of degree ≥ 1 and integer m ≥ 0 —
// Gosper's algorithm specialized to a factorial denominator. The antidifference, if
// it exists, is g(k) = Q(k)/(k+m−1)! with P(k)/(k+m)! = g(k) − g(k+1); multiplying
// by (k+m)! gives the polynomial identity Q(k)·(k+m) − Q(k+1) = P(k). Q has degree
// deg(P)−1 and is solved top-down; the constant-term equation is a consistency check
// that fails for non-telescoping terms (1/(k+1)! → e−2, handled elsewhere). The sum
// is then g(lo) − g(hi+1) (g(∞)=0 since the factorial dominates Q). Closes
// Σ_{k=1}^∞ k/(k+1)! = 1, Σ (k²−1)/(k+1)! = 1, …
[[nodiscard]] std::optional<Expr> sum_factorial_telescope(
    const Expr& expr, const Expr& var, const Expr& lo, const Expr& hi) {
    if (lo->type_id() != TypeId::Integer
        || !static_cast<const Integer&>(*lo).fits_long()
        || static_cast<const Integer&>(*lo).is_negative()) {
        return std::nullopt;
    }
    std::vector<Expr> factors;
    if (expr->type_id() == TypeId::Mul) {
        for (const auto& f : expr->args()) factors.push_back(f);
    } else {
        factors.push_back(expr);
    }
    Expr cst = S::One();
    std::vector<Expr> poly_factors;
    long m = 0;
    bool have_fact = false;
    for (const auto& f : factors) {
        if (f->type_id() == TypeId::Pow && f->args()[1] == S::NegativeOne()
            && f->args()[0]->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*f->args()[0]);
            if (fn.function_id() == FunctionId::Factorial
                && fn.args().size() == 1) {
                if (have_fact) return std::nullopt;  // (…!)^(-2) etc.
                Expr mm = simplify(fn.args()[0] + mul(S::NegativeOne(), var));
                if (has(mm, var) || mm->type_id() != TypeId::Integer
                    || !static_cast<const Integer&>(*mm).fits_long()) {
                    return std::nullopt;  // factorial arg not var + integer
                }
                m = static_cast<const Integer&>(*mm).to_long();
                if (m < 0) return std::nullopt;
                have_fact = true;
                continue;
            }
            return std::nullopt;  // some other reciprocal function of var
        }
        if (!has(f, var)) {
            cst = mul(cst, f);
            continue;
        }
        poly_factors.push_back(f);  // numerator polynomial part
    }
    if (!have_fact || poly_factors.empty()) return std::nullopt;
    // g(lo) needs (lo+m−1)! defined.
    if (static_cast<const Integer&>(*lo).to_long() + m - 1 < 0) return std::nullopt;

    std::vector<Expr> p;  // P coefficients, ascending
    try {
        Poly pp(expand(mul(poly_factors)), var);
        for (const auto& cc : pp.coeffs()) {
            if (has(cc, var)) return std::nullopt;
            p.push_back(cc);
        }
    } catch (const std::exception&) {
        return std::nullopt;
    }
    const long d = static_cast<long>(p.size()) - 1;  // deg P
    if (d < 1) return std::nullopt;  // constant numerator is not this telescoping

    // Solve Q(k)·(k+m) − Q(k+1) = P(k), Q of degree d−1, top-down:
    //   q_{d-1} = p_d ;  q_{t-1} = p_t − (m−1)·q_t + Σ_{j=t+1}^{d-1} C(j,t)·q_j .
    std::vector<Expr> q(static_cast<std::size_t>(d), Expr{S::Zero()});
    q[static_cast<std::size_t>(d - 1)] = p[static_cast<std::size_t>(d)];
    for (long t = d - 1; t >= 1; --t) {
        Expr s = p[static_cast<std::size_t>(t)];
        s = s + mul(integer(-(m - 1)), q[static_cast<std::size_t>(t)]);
        for (long j = t + 1; j <= d - 1; ++j) {
            s = s + mul(binomial(integer(j), integer(t)),
                        q[static_cast<std::size_t>(j)]);
        }
        q[static_cast<std::size_t>(t - 1)] = simplify(s);
    }
    // Consistency: m·q_0 − Σ_{j} q_j == p_0.
    Expr lhs0 = mul(integer(m), q[0]);
    for (long j = 0; j <= d - 1; ++j) {
        lhs0 = lhs0 + mul(S::NegativeOne(), q[static_cast<std::size_t>(j)]);
    }
    if (simplify(lhs0 + mul(S::NegativeOne(), p[0])) != S::Zero()) {
        return std::nullopt;  // not telescoping in this form
    }

    auto g = [&](const Expr& k) {
        std::vector<Expr> terms;
        for (long j = 0; j <= d - 1; ++j) {
            terms.push_back(
                mul(q[static_cast<std::size_t>(j)], pow(k, integer(j))));
        }
        return mul(add(std::move(terms)),
                   pow(factorial(k + integer(m - 1)), integer(-1)));
    };
    Expr result = (hi->type_id() == TypeId::Infinity)
                      ? g(lo)
                      : (g(lo) + mul(S::NegativeOne(), g(hi + integer(1))));
    return simplify(mul(cst, result));
}

// Σ_{k=0}^∞ c·[(−1)^k]·z^(2k+b)/(2k+b)! for b ∈ {0,1} — the even/odd bisection of
// the exponential series. Without the sign it gives cosh/sinh; with (−1)^k, cos/sin:
//   Σ z^(2k)/(2k)!     = cosh z      Σ (−1)^k z^(2k)/(2k)!     = cos z
//   Σ z^(2k+1)/(2k+1)! = sinh z      Σ (−1)^k z^(2k+1)/(2k+1)! = sin z
// z is the numerator base (z = 1 when the numerator is constant); its exponent must
// equal the factorial argument 2k+b. A lower bound lo > 0 subtracts the head terms.
[[nodiscard]] std::optional<Expr> sum_cosh_sinh_series(
    const Expr& expr, const Expr& var, const Expr& lo, const Expr& hi) {
    if (hi->type_id() != TypeId::Infinity) return std::nullopt;
    if (lo->type_id() != TypeId::Integer
        || static_cast<const Integer&>(*lo).is_negative()
        || !static_cast<const Integer&>(*lo).fits_long()) {
        return std::nullopt;
    }
    std::vector<Expr> factors;
    if (expr->type_id() == TypeId::Mul) {
        for (const auto& f : expr->args()) factors.push_back(f);
    } else {
        factors.push_back(expr);
    }
    Expr cst = S::One();
    Expr base;       // z
    Expr num_exp;    // exponent of the numerator base (must equal 2k+b)
    Expr fact_arg;   // 2k + b
    bool have_fact = false;
    bool have_sign = false;
    long b = 0;
    for (const auto& f : factors) {
        if (f->type_id() == TypeId::Pow && f->args()[1] == S::NegativeOne()
            && f->args()[0]->type_id() == TypeId::Function) {
            const auto& fn = static_cast<const Function&>(*f->args()[0]);
            if (fn.function_id() == FunctionId::Factorial
                && fn.args().size() == 1) {
                if (have_fact) return std::nullopt;
                Poly pa(expand(fn.args()[0]), var);
                if (pa.degree() != 1 || pa.coeffs()[1] != integer(2)) {
                    return std::nullopt;  // need a (2k+b)! denominator
                }
                const Expr& bb = pa.coeffs()[0];
                if (bb == S::Zero()) {
                    b = 0;
                } else if (bb == S::One()) {
                    b = 1;
                } else {
                    return std::nullopt;
                }
                fact_arg = fn.args()[0];
                have_fact = true;
                continue;
            }
            return std::nullopt;
        }
        // (−1)^(a·k + c): an odd a gives (−1)^k; (−1)^c is a constant sign.
        if (f->type_id() == TypeId::Pow && f->args()[0] == integer(-1)
            && has(f->args()[1], var)) {
            Poly ps(expand(f->args()[1]), var);
            if (ps.degree() != 1
                || ps.coeffs()[1]->type_id() != TypeId::Integer
                || ps.coeffs()[0]->type_id() != TypeId::Integer) {
                return std::nullopt;
            }
            const long a = static_cast<const Integer&>(*ps.coeffs()[1]).to_long();
            const long c = static_cast<const Integer&>(*ps.coeffs()[0]).to_long();
            if (a % 2 != 0) have_sign = !have_sign;  // toggle: (−1)^k present
            if (((c % 2) + 2) % 2 == 1) cst = mul(cst, S::NegativeOne());
            continue;
        }
        // z^(exponent in k), z var-free.
        if (f->type_id() == TypeId::Pow && !has(f->args()[0], var)
            && has(f->args()[1], var)) {
            if (base) return std::nullopt;
            base = f->args()[0];
            num_exp = f->args()[1];
            continue;
        }
        if (!has(f, var)) {
            cst = mul(cst, f);
            continue;
        }
        return std::nullopt;  // an unrecognized var-dependent factor
    }
    if (!have_fact) return std::nullopt;

    Expr z;
    if (base) {
        if (simplify(num_exp + mul(S::NegativeOne(), fact_arg)) != S::Zero()) {
            return std::nullopt;  // numerator exponent ≠ factorial argument
        }
        z = base;
    } else {
        z = S::One();
    }
    Expr full = have_sign ? (b == 0 ? cos(z) : sin(z))
                          : (b == 0 ? cosh(z) : sinh(z));
    full = mul(cst, full);

    const long L = static_cast<const Integer&>(*lo).to_long();
    Expr head = S::Zero();
    for (long kk = 0; kk < L; ++kk) {
        head = add({head, subs(expr, var, integer(kk))});
    }
    return simplify(full + mul(S::NegativeOne(), head));
}

// General convergent rational sum via partial fractions. apart() splits the summand
// into terms c·(a·k+b)^(-j): terms with j ≥ 2 each sum on their own (the p-series
// path gives ζ-values), and the j = 1 terms collectively telescope (their residues
// sum to zero for a convergent series), recombined and handed to telescope_rational.
// Closes Σ1/(k²(k+1)) = π²/6 − 1, mixing a ζ part and a telescoping part. Infinite
// range only — a finite j ≥ 2 part would need harmonic numbers.
[[nodiscard]] std::optional<Expr> sum_rational_via_apart(
    const Expr& expr, const Expr& var, const Expr& lo, const Expr& hi) {
    if (hi->type_id() != TypeId::Infinity) return std::nullopt;
    Expr a;
    try {
        a = apart(expr, var);
    } catch (const std::exception&) {
        return std::nullopt;
    }
    if (a->type_id() != TypeId::Add || a->args().size() < 2) return std::nullopt;

    // Classify a term as coeff·(linear)^(-j): returns 1 and sets j for a pole term,
    // 0 for a var-free constant, -1 for anything else (a positive var power, a
    // higher-degree base, or two poles in one term).
    auto classify = [&](const Expr& term, long& j) -> int {
        j = 0;
        std::vector<Expr> facs =
            term->type_id() == TypeId::Mul
                ? std::vector<Expr>(term->args().begin(), term->args().end())
                : std::vector<Expr>{term};
        bool have_pole = false;
        for (const auto& f : facs) {
            if (!has(f, var)) continue;
            if (f->type_id() == TypeId::Pow
                && f->args()[1]->type_id() == TypeId::Integer
                && static_cast<const Integer&>(*f->args()[1]).is_negative()
                && has(f->args()[0], var)) {
                if (have_pole) return -1;
                try {
                    Poly pb(expand(f->args()[0]), var);
                    if (pb.degree() != 1) return -1;
                } catch (const std::exception&) {
                    return -1;
                }
                j = -static_cast<const Integer&>(*f->args()[1]).to_long();
                have_pole = true;
            } else {
                return -1;
            }
        }
        if (have_pole) return 1;
        return has(term, var) ? -1 : 0;
    };

    std::vector<Expr> j1;
    Expr hi_sum = S::Zero();
    bool summed_any = false;
    for (const auto& term : a->args()) {
        long j = 0;
        const int cls = classify(term, j);
        if (cls < 0) return std::nullopt;
        if (cls == 0) {
            if (!(simplify(term) == S::Zero())) return std::nullopt;  // diverges
            continue;
        }
        if (j == 1) {
            j1.push_back(term);
            continue;
        }
        Expr s = summation(term, var, lo, hi);
        if (s->str().find("Sum(") != std::string::npos) return std::nullopt;
        hi_sum = add({hi_sum, s});
        summed_any = true;
    }
    if (!j1.empty()) {
        // simplify so the recombined numerator collapses to the var-free constant
        // telescope_rational expects (together alone may leave it as e.g. k−(k+1)).
        auto t = telescope_rational(simplify(together(add(j1))), var, lo, hi);
        if (!t) return std::nullopt;  // residues don't sum to zero → divergent
        hi_sum = add({hi_sum, *t});
        summed_any = true;
    }
    if (!summed_any) return std::nullopt;
    return simplify(hi_sum);
}

// Σ_{k=1}^∞ c/(a·k²+b) for a var-free quadratic denominator with no linear term
// and B = b/a > 0 (so the summand has no real pole). The classic Mittag-Leffler /
// cotangent result Σ_{k=1}^∞ 1/(k²+B) = (π·√B·coth(π·√B) − 1)/(2B). Apart can't
// reach it (the poles are an irreducible-quadratic conjugate pair, not rational),
// so it would otherwise stay unevaluated. Closes Σ1/(n²+1), Σ1/(n²+4), …
[[nodiscard]] std::optional<Expr> sum_inverse_quadratic(
    const Expr& expr, const Expr& var, const Expr& lo, const Expr& hi) {
    if (lo != S::One() || hi->type_id() != TypeId::Infinity) return std::nullopt;

    // Peel a var-free coefficient c off a (a·k²+b)^(-1) reciprocal factor.
    Expr coeff = S::One();
    const Expr* recip = nullptr;
    std::vector<Expr> factors =
        expr->type_id() == TypeId::Mul
            ? std::vector<Expr>(expr->args().begin(), expr->args().end())
            : std::vector<Expr>{expr};
    for (const auto& f : factors) {
        if (!has(f, var)) {
            coeff = mul(coeff, f);
        } else if (!recip && f->type_id() == TypeId::Pow
                   && f->args()[1] == S::NegativeOne()) {
            recip = &f;
        } else {
            return std::nullopt;  // a second var-bearing factor — not this form
        }
    }
    if (!recip) return std::nullopt;

    Poly den((*recip)->args()[0], var);
    if (den.degree() != 2) return std::nullopt;
    const auto& dc = den.coeffs();  // [b, lin, a]
    if (has(dc[0], var) || has(dc[2], var)) return std::nullopt;  // opaque coeffs
    if (!(dc[1] == S::Zero())) return std::nullopt;                // no linear term

    // B = b/a must be a positive number for the real-coth closed form.
    Expr B = simplify(dc[0] / dc[2]);
    if (is_positive(B) != std::optional<bool>{true}) return std::nullopt;

    Expr sq = sqrt(B);
    Expr value = (S::Pi() * sq * coth(S::Pi() * sq) - integer(1)) / (integer(2) * B);
    return simplify(coeff * value / dc[2]);
}

Expr summation(const Expr& expr, const Expr& var, const Expr& lo, const Expr& hi) {
    if (!expr) return expr;

    // Single-term range (hi == lo): Σ_{k=a}^{a} f(k) = f(a).
    if (hi == lo) return simplify(subs(expr, var, lo));

    // Explicit telescoping difference: a 2-term Add g(k) − g(k+1) sums to
    // g(lo) − g(hi+1). Checked before the linearity split, which would otherwise
    // separate it into two non-closing sums. Closes Σ(1/k − 1/(k+1)) = 1−1/(n+1),
    // Σ(1/k! − 1/(k+1)!), Σ(1/k² − 1/(k+1)²).
    if (expr->type_id() == TypeId::Add && expr->args().size() == 2) {
        const Expr& t0 = expr->args()[0];
        const Expr& t1 = expr->args()[1];
        Expr g;
        if (simplify(t1 + subs(t0, var, var + integer(1))) == S::Zero()) {
            g = t0;  // expr = g(k) − g(k+1)
        } else if (simplify(t0 + subs(t1, var, var + integer(1)))
                   == S::Zero()) {
            g = t1;
        }
        if (g) {
            return simplify(subs(g, var, lo) - subs(g, var, hi + integer(1)));
        }
    }

    // Linearity: split Add into separate sums.
    if (expr->type_id() == TypeId::Add) {
        std::vector<Expr> partial;
        partial.reserve(expr->args().size());
        for (const auto& term : expr->args()) {
            partial.push_back(summation(term, var, lo, hi));
        }
        return add(std::move(partial));
    }

    // Pull constant factors out of a Mul.
    if (expr->type_id() == TypeId::Mul) {
        std::vector<Expr> const_factors;
        std::vector<Expr> var_factors;
        for (const auto& f : expr->args()) {
            if (has(f, var)) var_factors.push_back(f);
            else const_factors.push_back(f);
        }
        if (!const_factors.empty() && !var_factors.empty()) {
            Expr inner = summation(mul(std::move(var_factors)), var, lo, hi);
            return mul(mul(std::move(const_factors)), inner);
        }
    }

    // Binomial theorem Σ_{k=0}^n C(n,k)·rᵏ = (1+r)ⁿ.
    if (auto bt = sum_binomial_theorem(expr, var, lo, hi)) return *bt;

    // Even/odd-index exponential series Σ z^(2k+b)/(2k+b)! → cosh/sinh/cos/sin.
    if (auto r = sum_cosh_sinh_series(expr, var, lo, hi)) return *r;

    // Factorial telescoping Σ P(k)/(k+m)! (deg P ≥ 1) — must run before the
    // expand-and-recurse below, which would split P(k) and destroy the telescoping
    // structure ((k²−1)/(k+1)! telescopes as a whole, not term by term).
    if (auto r = sum_factorial_telescope(expr, var, lo, hi)) return *r;

    // Repeated-root rational that telescopes after partial fractions: e.g.
    // (2k+1)/(k²(k+1)²) = 1/k² − 1/(k+1)², (3k²+3k+1)/(k³(k+1)³) = 1/k³ − 1/(k+1)³.
    // apart() exposes the g(k) − g(k+1) form that the explicit-difference check at the
    // top misses on the combined fraction (and telescope_rational skips for repeated
    // roots). Must run before the expand-and-recurse below, which would split the
    // numerator into non-telescoping pieces. Only the clean 2-term case is taken; a
    // residual ζ part (1/(k²(k+1)) → −1/k+1/(k+1)+1/k², a 3-term apart) falls through.
    if (expr->type_id() == TypeId::Mul || expr->type_id() == TypeId::Pow) {
        try {
            Expr a = apart(expr, var);
            if (a->type_id() == TypeId::Add && a->args().size() == 2) {
                const Expr& t0 = a->args()[0];
                const Expr& t1 = a->args()[1];
                Expr g;
                if (simplify(t1 + subs(t0, var, var + integer(1))) == S::Zero()) {
                    g = t0;
                } else if (simplify(t0 + subs(t1, var, var + integer(1)))
                           == S::Zero()) {
                    g = t1;
                }
                // Pole guard: an integer root of the summand's denominator that is
                // ≥ lo lies in the range, making the summand undefined; the
                // telescoping formula would otherwise emit a bogus singular value.
                bool safe = bool(g);
                if (safe) {
                    Expr denom = S::One();
                    std::vector<Expr> facs =
                        expr->type_id() == TypeId::Mul
                            ? std::vector<Expr>(expr->args().begin(),
                                                expr->args().end())
                            : std::vector<Expr>{expr};
                    for (const auto& f : facs) {
                        if (f->type_id() == TypeId::Pow
                            && f->args()[1]->type_id() == TypeId::Integer
                            && static_cast<const Integer&>(*f->args()[1])
                                   .is_negative()) {
                            const long e2 = -static_cast<const Integer&>(
                                                 *f->args()[1])
                                                 .to_long();
                            denom = mul(denom, pow(f->args()[0], integer(e2)));
                        }
                    }
                    if (has(denom, var)) {
                        try {
                            Poly pd(expand(denom), var);
                            for (const auto& r : pd.roots()) {
                                if (r->type_id() != TypeId::Integer) continue;
                                if (is_positive(simplify(lo - r))
                                    != std::optional<bool>{true}) {
                                    safe = false;  // root ≥ lo: pole in range
                                    break;
                                }
                            }
                        } catch (const std::exception&) {
                            // denominator roots not resolvable — keep it safe-ish:
                            // verify the endpoints are finite below.
                        }
                    }
                }
                if (safe) {
                    Expr res = simplify(subs(g, var, lo)
                                        - subs(g, var, hi + integer(1)));
                    if (res->type_id() != TypeId::ComplexInfinity
                        && res->type_id() != TypeId::NaN && !has(res, var)) {
                        return res;
                    }
                }
            }
        } catch (const std::exception&) {
            // not a rational function — fall through.
        }
    }

    // A product or power that expands to a sum — e.g. k·(k+1) → k²+k,
    // (k+1)² → k²+2k+1, (k+1)·2ᵏ → k·2ᵏ+2ᵏ — isn't matched by the closed forms
    // below. Expand and recurse so linearity splits it into terms each of which
    // (a monomial kᵖ, or a poly·geometric) is summable.
    if (expr->type_id() == TypeId::Mul || expr->type_id() == TypeId::Pow) {
        Expr ex = expand(expr);
        if (!(ex == expr) && ex->type_id() == TypeId::Add) {
            return summation(ex, var, lo, hi);
        }
    }

    // Constant w.r.t. var: c * (hi - lo + 1).
    if (!has(expr, var)) {
        return simplify(expr * (hi - lo + integer(1)));
    }

    // Identity Σ k from lo to hi = (hi-lo+1)(lo+hi)/2.
    if (expr == var) {
        return simplify((hi - lo + integer(1)) * (lo + hi) / integer(2));
    }

    // Σ kᵖ for a positive integer p via Faulhaber's formula, using Bernoulli
    // numbers (SymPy's B₁ = +1/2 convention):
    //   Σ_{k=1}^n kᵖ = 1/(p+1) · Σ_{j=0}^{p} C(p+1, j) B_j n^(p+1−j).
    if (expr->type_id() == TypeId::Pow && expr->args()[0] == var
        && expr->args()[1]->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*expr->args()[1]);
        if (z.fits_long()) {
            long p = z.to_long();
            if (p >= 2 && p <= 200) {  // p=1 handled above; cap for cost
                auto sum_to_n = [&](const Expr& n) -> Expr {
                    std::vector<Expr> terms;
                    terms.reserve(static_cast<std::size_t>(p) + 1);
                    for (long j = 0; j <= p; ++j) {
                        terms.push_back(mul(
                            {binomial(integer(p + 1), integer(j)),
                             bernoulli(integer(j)),
                             pow(n, integer(p + 1 - j))}));
                    }
                    return mul(pow(integer(p + 1), integer(-1)),
                               add(std::move(terms)));
                };
                Expr s_hi = sum_to_n(hi);
                Expr s_lo = sum_to_n(lo - integer(1));
                return simplify(s_hi - s_lo);
            }
        }
    }

    // Convergent p-series Σ_{k=1}^∞ 1/k^s = ζ(s) for an integer s ≥ 2. zeta()
    // closes the even cases (ζ(2)=π²/6, …) and keeps odd s as a symbolic ζ(s)
    // (matching SymPy's Sum(1/k**3).doit() = zeta(3)). The divergent harmonic
    // s=1 (m=-1) is excluded and falls through unevaluated.
    if (lo == S::One() && hi->type_id() == TypeId::Infinity
        && expr->type_id() == TypeId::Pow && expr->args()[0] == var
        && expr->args()[1]->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*expr->args()[1]);
        if (z.fits_long()) {
            const long m = z.to_long();  // summand is var^m
            if (m <= -2) return zeta(integer(-m));
        }
    }

    // Arithmetic-argument p-series Σ_{k=1}^∞ c/(a·k+b)^s for integer s ≥ 2 and
    // a ∈ {1, 2}. The denominator runs over one residue class, so the value is the
    // matching slice of ζ(s) minus the finitely many terms before the start:
    //   a=1, b ≥ 0      : ζ(s) − Σ_{n=1}^{b} n^(−s)
    //   a=2, b odd  ≥ −1: (1 − 2^(−s))·ζ(s) − Σ_{j=1}^{(b+1)/2} (2j−1)^(−s)   (odd n)
    //   a=2, b even ≥ 0 : 2^(−s)·ζ(s)        − Σ_{j=1}^{b/2}     (2j)^(−s)     (even n)
    // ζ(even) closes to a π^s rational (Σ1/(2k−1)²=π²/8, Σ1/(2k)²=π²/24); odd s
    // stays a symbolic ζ(s), as SymPy does. a ≥ 3 needs Hurwitz ζ and falls through.
    if (lo == S::One() && hi->type_id() == TypeId::Infinity) {
        Expr coeff = S::One();
        Expr powf;
        bool ok = true;
        bool have_pow = false;
        std::vector<Expr> factors;
        if (expr->type_id() == TypeId::Mul) {
            for (const auto& f : expr->args()) factors.push_back(f);
        } else {
            factors.push_back(expr);
        }
        for (const auto& f : factors) {
            if (f->type_id() == TypeId::Pow && has(f->args()[0], var)
                && f->args()[1]->type_id() == TypeId::Integer) {
                if (have_pow) { ok = false; break; }
                powf = f;
                have_pow = true;
            } else if (!has(f, var)) {
                coeff = mul(coeff, f);
            } else {
                ok = false;
                break;
            }
        }
        if (ok && have_pow) {
            const auto& zexp = static_cast<const Integer&>(*powf->args()[1]);
            if (zexp.fits_long() && zexp.to_long() <= -2) {
                const long s = -zexp.to_long();
                try {
                    Poly pb(expand(powf->args()[0]), var);
                    if (pb.degree() == 1
                        && pb.coeffs()[1]->type_id() == TypeId::Integer
                        && pb.coeffs()[0]->type_id() == TypeId::Integer
                        && static_cast<const Integer&>(*pb.coeffs()[1]).fits_long()
                        && static_cast<const Integer&>(*pb.coeffs()[0])
                               .fits_long()) {
                        const long a =
                            static_cast<const Integer&>(*pb.coeffs()[1]).to_long();
                        const long b =
                            static_cast<const Integer&>(*pb.coeffs()[0]).to_long();
                        auto term = [&](long n) {
                            return pow(integer(n), integer(-s));
                        };
                        std::optional<Expr> base_sum;
                        if (a == 1 && b >= 0) {
                            Expr sub = S::Zero();
                            for (long n = 1; n <= b; ++n) sub = add({sub, term(n)});
                            base_sum =
                                simplify(zeta(integer(s))
                                         + mul(S::NegativeOne(), sub));
                        } else if (a == 2 && b >= -1) {
                            const bool odd_n = (((b % 2) + 2) % 2) == 1;
                            Expr full = odd_n
                                ? mul(simplify(S::One()
                                               + mul(S::NegativeOne(),
                                                     pow(integer(2), integer(-s)))),
                                      zeta(integer(s)))
                                : mul(pow(integer(2), integer(-s)),
                                      zeta(integer(s)));
                            Expr sub = S::Zero();
                            if (odd_n) {
                                for (long j = 1; j <= (b + 1) / 2; ++j)
                                    sub = add({sub, term(2 * j - 1)});
                            } else {
                                for (long j = 1; j <= b / 2; ++j)
                                    sub = add({sub, term(2 * j)});
                            }
                            base_sum =
                                simplify(full + mul(S::NegativeOne(), sub));
                        }
                        if (base_sum) return simplify(mul(coeff, *base_sum));
                    }
                } catch (const std::exception&) {
                    // not a clean affine denominator — fall through
                }
            }
        }
    }

    // Alternating p-series Σ_{k=1}^∞ (−1)^k / k^s = (2^(1−s) − 1)·ζ(s) for s ≥ 2,
    // and −log 2 for s = 1 (the Dirichlet eta values η(s)). A (−1)^(a·k+b) factor
    // with odd a is (−1)^k up to the constant sign (−1)^b; a leading constant and
    // the sign multiply through.
    if (lo == S::One() && hi->type_id() == TypeId::Infinity) {
        Expr sign_exp;
        long m = 0;
        bool have_sign = false;
        bool have_pow = false;
        Expr coeff = S::One();
        bool ok = true;
        std::vector<Expr> factors;
        if (expr->type_id() == TypeId::Mul) {
            for (const auto& f : expr->args()) factors.push_back(f);
        } else {
            factors.push_back(expr);
        }
        for (const auto& f : factors) {
            if (f->type_id() == TypeId::Pow
                && f->args()[0] == integer(-1) && has(f->args()[1], var)) {
                if (have_sign) { ok = false; break; }
                sign_exp = f->args()[1];
                have_sign = true;
            } else if (f->type_id() == TypeId::Pow && f->args()[0] == var
                       && f->args()[1]->type_id() == TypeId::Integer) {
                const auto& zp = static_cast<const Integer&>(*f->args()[1]);
                if (have_pow || !zp.fits_long()) { ok = false; break; }
                m = zp.to_long();
                have_pow = true;
            } else if (!has(f, var)) {
                coeff = mul(coeff, f);
            } else {
                ok = false;
                break;
            }
        }
        if (ok && have_sign && have_pow && m <= -1) {
            // sign_exp = a·k + b, a an odd integer (so (−1)^(a·k) = (−1)^k) and b
            // an integer (so (−1)^b is a concrete ±1).
            try {
                Poly pe(expand(sign_exp), var);
                if (pe.degree() == 1
                    && pe.coeffs()[1]->type_id() == TypeId::Integer
                    && pe.coeffs()[0]->type_id() == TypeId::Integer) {
                    const long a = static_cast<const Integer&>(*pe.coeffs()[1])
                                       .to_long();
                    const long b = static_cast<const Integer&>(*pe.coeffs()[0])
                                       .to_long();
                    if (a % 2 != 0) {
                        const long s = -m;
                        Expr sign_b = (b % 2 == 0) ? S::One() : S::NegativeOne();
                        Expr base_sum =
                            (s == 1)
                                ? Expr{mul(S::NegativeOne(), log(integer(2)))}
                                : mul(simplify(pow(integer(2), integer(1 - s))
                                               - S::One()),
                                      zeta(integer(s)));
                        return simplify(mul(mul(coeff, sign_b), base_sum));
                    }
                }
            } catch (const std::exception&) {
                // not a clean affine sign exponent — fall through
            }
        }
    }

    // Dirichlet beta Σ_{k=0}^∞ (−1)^k/(2k+1)^s: β(1) = π/4 (Leibniz), β(2) =
    // Catalan's constant. Higher s have no elementary closed form (SymPy returns a
    // polylog), so only s ∈ {1, 2} are recognized. A (−1)^(a·k+b) factor with odd
    // a is (−1)^k up to the sign (−1)^b; a leading constant multiplies through.
    if (lo == S::Zero() && hi->type_id() == TypeId::Infinity) {
        Expr sign_exp;
        long s = 0;
        bool have_sign = false;
        bool have_base = false;
        Expr coeff = S::One();
        bool ok = true;
        std::vector<Expr> factors;
        if (expr->type_id() == TypeId::Mul) {
            for (const auto& f : expr->args()) factors.push_back(f);
        } else {
            factors.push_back(expr);
        }
        for (const auto& f : factors) {
            if (f->type_id() == TypeId::Pow && f->args()[0] == integer(-1)
                && has(f->args()[1], var)) {
                if (have_sign) { ok = false; break; }
                sign_exp = f->args()[1];
                have_sign = true;
            } else if (f->type_id() == TypeId::Pow
                       && f->args()[1]->type_id() == TypeId::Integer
                       && has(f->args()[0], var)) {
                // base must be the odd-denominator affine 2·var + 1.
                const auto& ze = static_cast<const Integer&>(*f->args()[1]);
                if (have_base || !ze.fits_long() || ze.to_long() >= 0) {
                    ok = false;
                    break;
                }
                try {
                    Poly pb(expand(f->args()[0]), var);
                    if (pb.degree() == 1 && pb.coeffs()[0] == S::One()
                        && pb.coeffs()[1] == integer(2)) {
                        s = -ze.to_long();
                        have_base = true;
                    } else {
                        ok = false;
                        break;
                    }
                } catch (const std::exception&) {
                    ok = false;
                    break;
                }
            } else if (!has(f, var)) {
                coeff = mul(coeff, f);
            } else {
                ok = false;
                break;
            }
        }
        if (ok && have_sign && have_base && (s == 1 || s == 2)) {
            try {
                Poly pe(expand(sign_exp), var);
                if (pe.degree() == 1
                    && pe.coeffs()[1]->type_id() == TypeId::Integer
                    && pe.coeffs()[0]->type_id() == TypeId::Integer) {
                    const long a = static_cast<const Integer&>(*pe.coeffs()[1])
                                       .to_long();
                    const long b = static_cast<const Integer&>(*pe.coeffs()[0])
                                       .to_long();
                    if (a % 2 != 0) {
                        Expr sign_b = (b % 2 == 0) ? S::One() : S::NegativeOne();
                        Expr val = (s == 1)
                                       ? Expr{S::Pi() / integer(4)}
                                       : S::Catalan();
                        return simplify(mul(mul(coeff, sign_b), val));
                    }
                }
            } catch (const std::exception&) {
                // not a clean affine sign exponent — fall through
            }
        }
    }

    // Exponential series Σ_{k=lo}^∞ r^k/k! = e^r (e.g. Σ 1/k! = e,
    // Σ x^k/k! = e^x, Σ 2^k/k! = e²). Convergent for every r.
    if (auto r = sum_exponential_series(expr, var, lo, hi)) return *r;

    // Geometric: base^(c*var + d) with base, c, d all independent of var.
    // Rewrites as base^d * (base^c)^var = A * ratio^var, a geometric series
    // with ratio = base^c and constant prefactor A = base^d. This subsumes
    // the plain base^var case (c=1, d=0) and also handles base^(-var),
    // base^(2*var), etc. — exponents the previous exact-match check missed.
    if (expr->type_id() == TypeId::Pow) {
        const Expr& base = expr->args()[0];
        const Expr& exponent = expr->args()[1];
        if (!has(base, var) && has(exponent, var)) {
            Expr c = diff(exponent, var);
            if (!has(c, var)) {  // exponent is linear in var
                Expr d = simplify(exponent - c * var);
                if (!has(d, var)) {
                    Expr ratio = pow(base, c);
                    Expr prefactor = pow(base, d);
                    // A * (ratio^lo - ratio^(hi+1)) / (1 - ratio)
                    Expr num = pow(ratio, lo) - pow(ratio, hi + integer(1));
                    Expr den = integer(1) - ratio;
                    return simplify(prefactor * num / den);
                }
            }
        }
    }

    // Arithmetic-geometric: Σ k·r^k where r = base^c is a concrete numeric
    // ratio ≠ 1 (k times a geometric term). Restricted to a numeric base so the
    // result is the plain closed form, matching SymPy (a symbolic ratio would
    // need a Piecewise on r = 1). An infinite upper bound is left to the general
    // poly·geometric path below, which has the |r| < 1 convergence handling.
    if (expr->type_id() == TypeId::Mul && expr->args().size() == 2
        && hi->type_id() != TypeId::Infinity) {
        const Expr& f0 = expr->args()[0];
        const Expr& f1 = expr->args()[1];
        // Identify the bare `var` factor and the geometric `base^(c·var+d)`.
        const Expr* geo = nullptr;
        if (f0 == var && f1->type_id() == TypeId::Pow) {
            geo = &f1;
        } else if (f1 == var && f0->type_id() == TypeId::Pow) {
            geo = &f0;
        }
        if (geo != nullptr) {
            const Expr& base = (*geo)->args()[0];
            const Expr& exponent = (*geo)->args()[1];
            if (is_number(base) && has(exponent, var)) {
                Expr c = diff(exponent, var);
                if (!has(c, var)) {  // exponent linear in var
                    Expr d = simplify(exponent - c * var);
                    Expr ratio = pow(base, c);
                    if (!has(d, var) && is_number(ratio)
                        && !(ratio == S::One())) {
                        // Σ_{k=0}^{N} k·r^k = r(1 − (N+1)r^N + N·r^{N+1})/(1−r)².
                        // The base^d prefactor of the geometric term factors out.
                        Expr prefactor = pow(base, d);
                        Expr one_minus_r_sq =
                            pow(integer(1) - ratio, integer(2));
                        auto partial = [&](const Expr& n) -> Expr {
                            return ratio
                                   * (integer(1) - (n + integer(1)) * pow(ratio, n)
                                      + n * pow(ratio, n + integer(1)))
                                   / one_minus_r_sq;
                        };
                        Expr s = partial(hi) - partial(lo - integer(1));
                        return simplify(prefactor * s);
                    }
                }
            }
        }
    }

    // Polynomial × geometric: Σ P(k)·r^k for a polynomial P and concrete ratio
    // r ≠ 1 (covers k²·2^k, (2k+1)·3^k, k³/2^k — the general arithmetic-geometric
    // case the degree-1 block above doesn't reach).
    if (auto pg = sum_poly_geometric(expr, var, lo, hi); pg) return *pg;

    // Telescoping rational summand (e.g. 1/(k(k+1)), 1/(4k²−1)).
    if (auto t = telescope_rational(expr, var, lo, hi); t) return *t;

    // Σ c/(a·k²+b) with b/a > 0 — irreducible-quadratic denominator (cotangent
    // closed form). Checked before apart, which can't split the complex poles.
    if (auto r = sum_inverse_quadratic(expr, var, lo, hi)) return *r;

    // General convergent rational sum: apart into a ζ part (j ≥ 2 poles) and a
    // telescoping part (j = 1 poles). Closes Σ1/(k²(k+1)) = π²/6 − 1.
    if (auto r = sum_rational_via_apart(expr, var, lo, hi)) return *r;

    // Index-shift fallback for an infinite sum with an integer start ≠ 1:
    // Σ_{n=lo}^∞ f(n) = Σ_{m=1}^∞ f(m + lo − 1). Re-expressing from m=1 lets the
    // many lo=1 handlers above (arithmetic p-series, cotangent, …) reach a series
    // merely written from a different start — e.g. the standard odd p-series
    // Σ_{n=0}^∞ 1/(2n+1)² = Σ_{m=1}^∞ 1/(2m−1)² = π²/8. The shifted call has lo=1
    // so it cannot re-enter here; its result is adopted only when it is a true
    // closed form (var-free), never an unevaluated Sum (which still carries var).
    if (hi->type_id() == TypeId::Infinity && lo->type_id() == TypeId::Integer
        && !(lo == S::One())) {
        Expr shifted = summation(subs(expr, var, var + (lo - integer(1))), var,
                                 S::One(), hi);
        if (!has(shifted, var)) return shifted;
    }

    // No closed form found — return the unevaluated Sum marker rather than the
    // bare summand (Σ 1/k² must not collapse to 1/k²).
    return sum_marker(expr, var, lo, hi);
}

}  // namespace sympp
