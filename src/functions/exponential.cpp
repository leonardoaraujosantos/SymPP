#include <sympp/functions/exponential.hpp>

#include <optional>
#include <utility>
#include <vector>

#include <mpfr.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
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

    // exp(log(x)) = x, when x is positive (otherwise needs branch handling).
    if (arg->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*arg);
        if (fn.function_id() == FunctionId::Log) {
            const auto& inner = arg->args()[0];
            if (is_positive(inner) == true) return inner;
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

    // log(exp(x)) = x, when x is real.
    if (arg->type_id() == TypeId::Function) {
        const auto& fn = static_cast<const Function&>(*arg);
        if (fn.function_id() == FunctionId::Exp) {
            const auto& inner = arg->args()[0];
            if (is_real(inner) == true) return inner;
        }
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
