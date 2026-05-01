#include <sympp/functions/miscellaneous.hpp>

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
