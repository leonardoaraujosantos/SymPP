#include <sympp/functions/hyperbolic.hpp>

#include <optional>
#include <utility>
#include <vector>

#include <mpfr.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

// Detect a leading -1 factor on a Mul: -1 * rest. Returns the stripped tail
// if so, std::nullopt otherwise.
[[nodiscard]] std::optional<Expr> strip_neg(const Expr& e) {
    if (is_number(e)) {
        const auto& n = static_cast<const Number&>(*e);
        if (n.is_negative()) {
            return mul(S::NegativeOne(), e);
        }
        return std::nullopt;
    }
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

[[nodiscard]] std::optional<bool> ask_real_finite_for_real_arg(
        const Expr& arg, AssumptionKey k) noexcept {
    switch (k) {
        case AssumptionKey::Real:
        case AssumptionKey::Finite:
            if (is_real(arg) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

}  // namespace

// ----- Sinh ------------------------------------------------------------------

Sinh::Sinh(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Sinh);
}
Expr Sinh::rebuild(std::vector<Expr> new_args) const { return sinh(new_args[0]); }
std::optional<bool> Sinh::ask(AssumptionKey k) const noexcept {
    return ask_real_finite_for_real_arg(args_[0], k);
}

// ----- Cosh ------------------------------------------------------------------

Cosh::Cosh(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Cosh);
}
Expr Cosh::rebuild(std::vector<Expr> new_args) const { return cosh(new_args[0]); }
std::optional<bool> Cosh::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Real:
        case AssumptionKey::Finite:
            if (is_real(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Positive:
            // cosh(real) >= 1 > 0
            if (is_real(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Nonzero:
            if (is_real(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

// ----- Tanh ------------------------------------------------------------------

Tanh::Tanh(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Tanh);
}
Expr Tanh::rebuild(std::vector<Expr> new_args) const { return tanh(new_args[0]); }
std::optional<bool> Tanh::ask(AssumptionKey k) const noexcept {
    return ask_real_finite_for_real_arg(args_[0], k);
}

// ----- Asinh -----------------------------------------------------------------

Asinh::Asinh(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Asinh);
}
Expr Asinh::rebuild(std::vector<Expr> new_args) const { return asinh(new_args[0]); }
std::optional<bool> Asinh::ask(AssumptionKey k) const noexcept {
    return ask_real_finite_for_real_arg(args_[0], k);
}

// ----- Acosh -----------------------------------------------------------------

Acosh::Acosh(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Acosh);
}
Expr Acosh::rebuild(std::vector<Expr> new_args) const { return acosh(new_args[0]); }
std::optional<bool> Acosh::ask(AssumptionKey /*k*/) const noexcept {
    // acosh real iff arg >= 1; without that knowledge we cannot guarantee.
    return std::nullopt;
}

// ----- Atanh -----------------------------------------------------------------

Atanh::Atanh(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Atanh);
}
Expr Atanh::rebuild(std::vector<Expr> new_args) const { return atanh(new_args[0]); }
std::optional<bool> Atanh::ask(AssumptionKey /*k*/) const noexcept {
    return std::nullopt;
}

// ----- Factories -------------------------------------------------------------

Expr sinh(const Expr& arg) {
    if (arg == S::Zero()) return S::Zero();
    if (arg->type_id() == TypeId::Float) {
        return unary_evalf(mpfr_sinh, arg);
    }
    if (auto pos = strip_neg(arg); pos.has_value()) {
        return mul(S::NegativeOne(), make<Sinh>(*pos));
    }
    return make<Sinh>(arg);
}

Expr cosh(const Expr& arg) {
    if (arg == S::Zero()) return S::One();
    if (arg->type_id() == TypeId::Float) {
        return unary_evalf(mpfr_cosh, arg);
    }
    if (auto pos = strip_neg(arg); pos.has_value()) {
        return make<Cosh>(*pos);  // even
    }
    return make<Cosh>(arg);
}

Expr tanh(const Expr& arg) {
    if (arg == S::Zero()) return S::Zero();
    if (arg->type_id() == TypeId::Float) {
        return unary_evalf(mpfr_tanh, arg);
    }
    if (auto pos = strip_neg(arg); pos.has_value()) {
        return mul(S::NegativeOne(), make<Tanh>(*pos));
    }
    return make<Tanh>(arg);
}

Expr asinh(const Expr& arg) {
    if (arg == S::Zero()) return S::Zero();
    if (arg->type_id() == TypeId::Float) {
        return unary_evalf(mpfr_asinh, arg);
    }
    if (auto pos = strip_neg(arg); pos.has_value()) {
        return mul(S::NegativeOne(), make<Asinh>(*pos));
    }
    return make<Asinh>(arg);
}

Expr acosh(const Expr& arg) {
    if (arg == S::One()) return S::Zero();
    if (arg->type_id() == TypeId::Float) {
        return unary_evalf(mpfr_acosh, arg);
    }
    return make<Acosh>(arg);
}

Expr atanh(const Expr& arg) {
    if (arg == S::Zero()) return S::Zero();
    if (arg->type_id() == TypeId::Float) {
        return unary_evalf(mpfr_atanh, arg);
    }
    if (auto pos = strip_neg(arg); pos.has_value()) {
        return mul(S::NegativeOne(), make<Atanh>(*pos));
    }
    return make<Atanh>(arg);
}

}  // namespace sympp
