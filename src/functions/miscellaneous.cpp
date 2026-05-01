#include <sympp/functions/miscellaneous.hpp>

#include <algorithm>
#include <optional>
#include <utility>
#include <vector>

#include <gmpxx.h>
#include <mpfr.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

// Detect a leading minus sign on a Mul: -1 * rest. Returns the stripped tail
// if so, std::nullopt otherwise. Number negatives stay as-is for the abs path
// (we handle them separately).
[[nodiscard]] std::optional<Expr> strip_neg_factor(const Expr& e) {
    if (e->type_id() == TypeId::Mul) {
        const auto& args = e->args();
        if (!args.empty() && args[0]->type_id() == TypeId::Integer) {
            const auto& n = static_cast<const Integer&>(*args[0]);
            if (n.value() == -1) {
                if (args.size() == 2) return args[1];
                std::vector<Expr> rest(args.begin() + 1, args.end());
                return mul(std::move(rest));
            }
        }
    }
    return std::nullopt;
}

// Negate a Number — return its absolute value.
[[nodiscard]] Expr abs_number(const Number& n) {
    switch (n.type_id()) {
        case TypeId::Integer: {
            mpz_class v = static_cast<const Integer&>(n).value();
            if (sgn(v) < 0) v = -v;
            return make<Integer>(std::move(v));
        }
        case TypeId::Rational: {
            mpq_class q = static_cast<const Rational&>(n).value();
            if (sgn(q) < 0) q = -q;
            q.canonicalize();
            if (q.get_den() == 1) return make<Integer>(q.get_num());
            return make<Rational>(std::move(q));
        }
        case TypeId::Float: {
            const auto& f = static_cast<const Float&>(n);
            int dps = f.precision_dps();
            mpfr_t r;
            mpfr_init2(r, dps_to_prec(dps));
            mpfr_abs(r, f.value(), MPFR_RNDN);
            auto out = make<Float>(static_cast<mpfr_srcptr>(r), dps);
            mpfr_clear(r);
            return out;
        }
        default:
            // Should not happen in practice — Number subclasses are exhaustive.
            return Expr{};
    }
}

}  // namespace

// ----- Abs -------------------------------------------------------------------

Abs::Abs(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Abs);
}

Expr Abs::rebuild(std::vector<Expr> new_args) const {
    return abs(new_args[0]);
}

std::optional<bool> Abs::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Real:
            return true;  // |x| is always real (for any complex x, |x| ∈ ℝ_≥0)
        case AssumptionKey::Negative:
            return false;
        case AssumptionKey::Nonnegative:
            return true;
        case AssumptionKey::Positive:
            if (is_nonzero(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Zero:
            if (is_zero(a) == true) return true;
            if (is_nonzero(a) == true) return false;
            return std::nullopt;
        case AssumptionKey::Integer:
            if (is_integer(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Rational:
            if (is_rational(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Finite:
            if (is_finite(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Expr abs(const Expr& arg) {
    // Numeric: take the absolute value directly.
    if (is_number(arg)) {
        return abs_number(static_cast<const Number&>(*arg));
    }

    // NumberSymbol — Pi, E, EulerGamma, Catalan are all positive.
    if (arg->type_id() == TypeId::NumberSymbol) {
        return arg;
    }

    // Abs(x) for nonnegative x → x.
    if (is_nonnegative(arg) == true) return arg;
    // Abs(x) for nonpositive x → -x.
    if (is_nonpositive(arg) == true) return mul(S::NegativeOne(), arg);

    // Abs(-x) = Abs(x).
    if (auto pos = strip_neg_factor(arg); pos.has_value()) {
        return make<Abs>(*pos);
    }

    return make<Abs>(arg);
}

// ----- Sign ------------------------------------------------------------------

Sign::Sign(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Sign);
}

Expr Sign::rebuild(std::vector<Expr> new_args) const {
    return sign(new_args[0]);
}

std::optional<bool> Sign::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Real:
            if (is_real(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Integer:
            if (is_real(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Finite:
            return true;
        default:
            return std::nullopt;
    }
}

Expr sign(const Expr& arg) {
    // Numeric: use the sign predicate.
    if (is_number(arg)) {
        int s = static_cast<const Number&>(*arg).sign();
        if (s > 0) return S::One();
        if (s < 0) return S::NegativeOne();
        return S::Zero();
    }

    // NumberSymbol → +1 (Pi, E, etc.).
    if (arg->type_id() == TypeId::NumberSymbol) {
        return S::One();
    }

    if (is_zero(arg) == true) return S::Zero();
    if (is_positive(arg) == true) return S::One();
    if (is_negative(arg) == true) return S::NegativeOne();

    // sign(-x) = -sign(x)
    if (auto pos = strip_neg_factor(arg); pos.has_value()) {
        return mul(S::NegativeOne(), make<Sign>(*pos));
    }

    return make<Sign>(arg);
}

// ----- sqrt — uses Pow under the hood -----------------------------------------

// ============================================================================
// Complex part extraction
// ============================================================================

Re::Re(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Re);
}
Expr Re::rebuild(std::vector<Expr> new_args) const { return re(new_args[0]); }
std::optional<bool> Re::ask(AssumptionKey k) const noexcept {
    switch (k) {
        case AssumptionKey::Real:
            return true;  // re(z) is always real for any complex z
        case AssumptionKey::Finite:
            if (is_finite(args_[0]) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Im::Im(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Im);
}
Expr Im::rebuild(std::vector<Expr> new_args) const { return im(new_args[0]); }
std::optional<bool> Im::ask(AssumptionKey k) const noexcept {
    switch (k) {
        case AssumptionKey::Real:
            return true;  // im(z) ∈ ℝ
        case AssumptionKey::Finite:
            if (is_finite(args_[0]) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Conjugate::Conjugate(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Conjugate);
}
Expr Conjugate::rebuild(std::vector<Expr> new_args) const {
    return conjugate(new_args[0]);
}
std::optional<bool> Conjugate::ask(AssumptionKey k) const noexcept {
    // Conjugation preserves real-ness, integer-ness, etc. on real arguments;
    // mostly mirrors arg's properties.
    return args_[0]->ask(k);
}

Arg::Arg(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Arg);
}
Expr Arg::rebuild(std::vector<Expr> new_args) const { return arg_(new_args[0]); }
std::optional<bool> Arg::ask(AssumptionKey k) const noexcept {
    switch (k) {
        case AssumptionKey::Real:
            return true;
        default:
            return std::nullopt;
    }
}

Expr re(const Expr& arg) {
    // Numeric: real -> identity; (no Complex type yet, so no a+bi case)
    if (is_real(arg) == true) return arg;
    // Linearity: re(-x) = -re(x)
    if (auto p = strip_neg_factor(arg); p.has_value()) {
        return mul(S::NegativeOne(), make<Re>(*p));
    }
    return make<Re>(arg);
}

Expr im(const Expr& arg) {
    if (is_real(arg) == true) return S::Zero();
    if (auto p = strip_neg_factor(arg); p.has_value()) {
        return mul(S::NegativeOne(), make<Im>(*p));
    }
    return make<Im>(arg);
}

Expr conjugate(const Expr& arg) {
    if (is_real(arg) == true) return arg;
    // conjugate(conjugate(x)) = x
    if (arg->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*arg);
        if (fn.function_id() == FunctionId::Conjugate) {
            return arg->args()[0];
        }
    }
    return make<Conjugate>(arg);
}

Expr arg_(const Expr& arg) {
    if (is_positive(arg) == true) return S::Zero();
    if (is_negative(arg) == true) return S::Pi();
    return make<Arg>(arg);
}

// ============================================================================
// Min, Max
// ============================================================================

MinFn::MinFn(std::vector<Expr> args)
    : Function(std::move(args)) {
    compute_hash(FunctionId::Min);
}
Expr MinFn::rebuild(std::vector<Expr> new_args) const { return min(std::move(new_args)); }
std::optional<bool> MinFn::ask(AssumptionKey k) const noexcept {
    // Min preserves real / integer / rational under all-true closure.
    switch (k) {
        case AssumptionKey::Real:
        case AssumptionKey::Integer:
        case AssumptionKey::Rational:
        case AssumptionKey::Finite:
            for (const auto& a : args_) {
                if (a->ask(k) != true) return std::nullopt;
            }
            return true;
        default:
            return std::nullopt;
    }
}

MaxFn::MaxFn(std::vector<Expr> args)
    : Function(std::move(args)) {
    compute_hash(FunctionId::Max);
}
Expr MaxFn::rebuild(std::vector<Expr> new_args) const { return max(std::move(new_args)); }
std::optional<bool> MaxFn::ask(AssumptionKey k) const noexcept {
    switch (k) {
        case AssumptionKey::Real:
        case AssumptionKey::Integer:
        case AssumptionKey::Rational:
        case AssumptionKey::Finite:
            for (const auto& a : args_) {
                if (a->ask(k) != true) return std::nullopt;
            }
            return true;
        default:
            return std::nullopt;
    }
}

namespace {

// Compare two Number Exprs. Returns -1, 0, +1. Caller verifies both are Numbers.
[[nodiscard]] int compare_numbers(const Number& a, const Number& b) noexcept {
    // Promote to mpq for the exact tower; Float comparisons compute via mpfr.
    if (a.type_id() == TypeId::Float || b.type_id() == TypeId::Float) {
        // At least one is Float — promote the other to Float at the higher dps.
        int dps = kDefaultDps;
        if (a.type_id() == TypeId::Float) {
            dps = std::max(dps, static_cast<const Float&>(a).precision_dps());
        }
        if (b.type_id() == TypeId::Float) {
            dps = std::max(dps, static_cast<const Float&>(b).precision_dps());
        }
        mpfr_t ax, bx;
        mpfr_init2(ax, dps_to_prec(dps));
        mpfr_init2(bx, dps_to_prec(dps));

        if (a.type_id() == TypeId::Float) {
            mpfr_set(ax, static_cast<const Float&>(a).value(), MPFR_RNDN);
        } else if (a.type_id() == TypeId::Integer) {
            mpfr_set_z(ax, static_cast<const Integer&>(a).value().get_mpz_t(), MPFR_RNDN);
        } else {
            mpfr_set_q(ax, static_cast<const Rational&>(a).value().get_mpq_t(), MPFR_RNDN);
        }
        if (b.type_id() == TypeId::Float) {
            mpfr_set(bx, static_cast<const Float&>(b).value(), MPFR_RNDN);
        } else if (b.type_id() == TypeId::Integer) {
            mpfr_set_z(bx, static_cast<const Integer&>(b).value().get_mpz_t(), MPFR_RNDN);
        } else {
            mpfr_set_q(bx, static_cast<const Rational&>(b).value().get_mpq_t(), MPFR_RNDN);
        }
        int c = mpfr_cmp(ax, bx);
        mpfr_clear(ax);
        mpfr_clear(bx);
        return c < 0 ? -1 : (c > 0 ? 1 : 0);
    }
    // Both exact: promote to mpq.
    mpq_class qa, qb;
    if (a.type_id() == TypeId::Integer) qa = static_cast<const Integer&>(a).value();
    else qa = static_cast<const Rational&>(a).value();
    if (b.type_id() == TypeId::Integer) qb = static_cast<const Integer&>(b).value();
    else qb = static_cast<const Rational&>(b).value();
    int c = cmp(qa, qb);
    return c < 0 ? -1 : (c > 0 ? 1 : 0);
}

// Reduce a Min/Max args list by combining all numeric args into one extreme,
// keeping symbolic args as-is.
template <bool IsMin>
[[nodiscard]] Expr min_max_impl(std::vector<Expr> args) {
    std::vector<Expr> non_numeric;
    Expr extreme;  // running min or max numeric
    for (auto& a : args) {
        if (is_number(a)) {
            if (!extreme) {
                extreme = a;
            } else {
                int c = compare_numbers(static_cast<const Number&>(*extreme),
                                       static_cast<const Number&>(*a));
                if constexpr (IsMin) {
                    if (c > 0) extreme = a;  // a < current min
                } else {
                    if (c < 0) extreme = a;  // a > current max
                }
            }
        } else {
            non_numeric.push_back(std::move(a));
        }
    }

    // Drop duplicate symbolic args — Min(x, x) = x.
    std::sort(non_numeric.begin(), non_numeric.end(),
              [](const Expr& l, const Expr& r) noexcept {
                  return l.get() < r.get();
              });
    non_numeric.erase(std::unique(non_numeric.begin(), non_numeric.end()),
                      non_numeric.end());

    std::vector<Expr> out;
    out.reserve(non_numeric.size() + 1);
    if (extreme) out.push_back(std::move(extreme));
    for (auto& a : non_numeric) out.push_back(std::move(a));

    if (out.empty()) {
        // Empty Min/Max — by convention, undefined; SymPy raises. Return
        // a NaN-equivalent? For now build the empty form to surface the
        // misuse; callers rarely hit this with proper inputs.
    }
    if (out.size() == 1) return std::move(out[0]);

    if constexpr (IsMin) return make<MinFn>(std::move(out));
    else return make<MaxFn>(std::move(out));
}

}  // namespace

Expr min(std::vector<Expr> args) { return min_max_impl<true>(std::move(args)); }
Expr min(const Expr& a, const Expr& b) {
    return min_max_impl<true>(std::vector<Expr>{a, b});
}
Expr max(std::vector<Expr> args) { return min_max_impl<false>(std::move(args)); }
Expr max(const Expr& a, const Expr& b) {
    return min_max_impl<false>(std::vector<Expr>{a, b});
}

// ============================================================================
// sqrt — uses Pow under the hood
// ============================================================================

Expr sqrt(const Expr& arg) {
    if (arg == S::Zero()) return S::Zero();
    if (arg == S::One()) return S::One();

    // Perfect-square integers fold to integers when they're nonnegative.
    if (arg->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*arg);
        if (!z.is_negative()) {
            mpz_class s;
            mpz_class r;
            mpz_sqrtrem(s.get_mpz_t(), r.get_mpz_t(), z.value().get_mpz_t());
            if (sgn(r) == 0) {
                return make<Integer>(std::move(s));
            }
        }
    }

    return pow(arg, S::Half());
}

}  // namespace sympp
