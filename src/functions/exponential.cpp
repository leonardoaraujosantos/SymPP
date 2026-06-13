#include <sympp/functions/exponential.hpp>

#include <optional>
#include <utility>
#include <vector>

#include <mpfr.h>

#include <sympp/core/add.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/infinity.hpp>
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

[[nodiscard]] Expr unary_evalf(int (*op)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t),
                               const Expr& arg) {
    const auto& f = static_cast<const Float&>(*arg);
    int dps = f.precision_dps();
    mpfr_t r;
    mpfr_init2(r, dps_to_prec(dps));
    op(r, f.value(), MPFR_RNDN);
    auto out = make<Float>(static_cast<mpfr_srcptr>(r), dps);
    mpfr_clear(r);
    return out;
}

// If `arg` is r·I·π for a rational r, return r; else nullopt. Used for the
// Euler identity exp(r·I·π).
[[nodiscard]] std::optional<mpq_class> imaginary_pi_coeff(const Expr& arg) {
    if (arg->type_id() != TypeId::Mul) return std::nullopt;
    bool has_i = false;
    bool has_pi = false;
    Expr coeff = S::One();
    for (const auto& f : arg->args()) {
        if (f == S::I()) {
            if (has_i) return std::nullopt;  // I² etc.
            has_i = true;
        } else if (f == S::Pi()) {
            if (has_pi) return std::nullopt;  // π² etc.
            has_pi = true;
        } else if (is_number(f)) {
            coeff = mul(coeff, f);
        } else {
            return std::nullopt;  // a non-numeric extra factor
        }
    }
    if (!has_i || !has_pi) return std::nullopt;
    if (coeff->type_id() == TypeId::Integer) {
        return mpq_class(static_cast<const Integer&>(*coeff).value());
    }
    if (coeff->type_id() == TypeId::Rational) {
        return static_cast<const Rational&>(*coeff).value();
    }
    return std::nullopt;
}

// If `arg` is b·I for a nonzero rational b (including I itself), return b as a
// numeric Expr; else nullopt. Used for log on the imaginary axis.
[[nodiscard]] std::optional<Expr> imaginary_coeff(const Expr& arg) {
    if (arg == S::I()) return S::One();
    if (arg->type_id() != TypeId::Mul) return std::nullopt;
    bool has_i = false;
    Expr coeff = S::One();
    for (const auto& f : arg->args()) {
        if (f == S::I()) {
            if (has_i) return std::nullopt;  // I² etc.
            has_i = true;
        } else if (f->type_id() == TypeId::Integer
                   || f->type_id() == TypeId::Rational) {
            coeff = mul(coeff, f);
        } else {
            return std::nullopt;  // non-rational extra factor
        }
    }
    if (!has_i) return std::nullopt;
    return coeff;
}

}  // namespace

// ----- Exp -------------------------------------------------------------------

Exp::Exp(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Exp);
}

Expr Exp::rebuild(std::vector<Expr> new_args) const {
    return exp(new_args[0]);
}

std::optional<bool> Exp::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Real:
            if (is_real(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Positive:
            if (is_real(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Nonzero:
            if (is_real(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Finite:
            if (is_real(a) == true && is_finite(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

// ----- Log -------------------------------------------------------------------

Log::Log(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Log);
}

Expr Log::rebuild(std::vector<Expr> new_args) const {
    return log(new_args[0]);
}

std::optional<bool> Log::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Real:
            // log is real iff arg > 0.
            if (is_positive(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

// ----- Factories -------------------------------------------------------------

Expr exp(const Expr& arg) {
    if (arg == S::Zero()) return S::One();
    if (arg == S::One()) return S::E();

    // Infinity: exp(oo) = oo, exp(-oo) = 0, exp(zoo)/exp(nan) = nan.
    if (arg->type_id() == TypeId::Infinity) return S::Infinity();
    if (arg->type_id() == TypeId::NegativeInfinity) return S::Zero();
    if (arg->type_id() == TypeId::ComplexInfinity
        || arg->type_id() == TypeId::NaN) {
        return S::NaN();
    }

    // exp(log(x)) = x, when x is positive (otherwise needs branch handling).
    if (arg->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*arg);
        if (fn.function_id() == FunctionId::Log) {
            const auto& inner = arg->args()[0];
            if (is_positive(inner) == true) return inner;
        }
    }

    // exp(c·log(p)) = p^c for positive p (any real exponent c). Same positivity
    // gate as exp(log p) above — for non-positive p this is branch-cut sensitive.
    if (arg->type_id() == TypeId::Mul) {
        Expr log_inner;
        std::vector<Expr> coeff_factors;
        for (const auto& f : arg->args()) {
            if (!log_inner && f->type_id() == TypeId::Function) {
                const auto& fn = static_cast<const Function&>(*f);
                if (fn.function_id() == FunctionId::Log && fn.args().size() == 1) {
                    log_inner = fn.args()[0];
                    continue;
                }
            }
            coeff_factors.push_back(f);
        }
        if (log_inner && !coeff_factors.empty()
            && is_positive(log_inner) == true) {
            return pow(log_inner, mul(std::move(coeff_factors)));
        }
    }

    // Euler: exp(r·I·π) = i^(2r) when 2r is an integer (q ∈ {1,2}: ±1, ±I).
    // Matches SymPy, which keeps π/3, π/4, … exponents symbolic. pow(I, n)
    // already cycles I^n through {1, I, -1, -I}.
    if (auto r = imaginary_pi_coeff(arg); r.has_value()) {
        mpq_class two_r = 2 * *r;
        if (two_r.get_den() == 1) {
            mpz_class n = two_r.get_num() % 4;
            if (n < 0) n += 4;
            return pow(S::I(), integer(n.get_si()));
        }
    }

    if (arg->type_id() == TypeId::Float) {
        return unary_evalf(mpfr_exp, arg);
    }

    return make<Exp>(arg);
}

Expr log(const Expr& arg) {
    if (arg == S::One()) return S::Zero();
    if (arg == S::E()) return S::One();

    // Infinity: log(oo) = oo, log(0) = zoo, log(-oo) = oo, log(zoo) = oo.
    if (arg->type_id() == TypeId::Infinity
        || arg->type_id() == TypeId::NegativeInfinity
        || arg->type_id() == TypeId::ComplexInfinity) {
        return S::Infinity();
    }
    if (arg->type_id() == TypeId::NaN) return S::NaN();
    if (arg == S::Zero()) return S::ComplexInfinity();

    // log(exp(x)) = x, when x is real.
    if (arg->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*arg);
        if (fn.function_id() == FunctionId::Exp) {
            const auto& inner = arg->args()[0];
            if (is_real(inner) == true) return inner;
        }
    }

    // log(negative real) = log(|x|) + I·π  (covers −1, −2, −1/2, −E, …).
    if (is_real(arg) == true && is_negative(arg) == true) {
        return add(log(mul(S::NegativeOne(), arg)), mul(S::I(), S::Pi()));
    }

    // log(b·I) = log(|b|) + sign(b)·I·π/2 for a nonzero rational b.
    if (auto b = imaginary_coeff(arg); b.has_value()) {
        bool neg = is_negative(*b) == true;
        Expr magnitude = neg ? mul(S::NegativeOne(), *b) : *b;
        Expr quarter_turn = mul(rational(neg ? -1 : 1, 2), mul(S::I(), S::Pi()));
        return add(log(magnitude), quarter_turn);
    }

    if (arg->type_id() == TypeId::Float) {
        // log of negative or zero would produce NaN/-inf via MPFR; let it
        // through — the result Float will encode that and downstream can
        // surface it.
        return unary_evalf(mpfr_log, arg);
    }

    return make<Log>(arg);
}

// ----- Derivatives ----------------------------------------------------------

Expr Exp::diff_arg(std::size_t /*i*/) const {
    // d/dx exp(x) = exp(x)
    return exp(args_[0]);
}
Expr Log::diff_arg(std::size_t /*i*/) const {
    // d/dx log(x) = 1/x
    return pow(args_[0], S::NegativeOne());
}

}  // namespace sympp
