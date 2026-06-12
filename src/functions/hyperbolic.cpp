#include <sympp/functions/hyperbolic.hpp>

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

// If `arg` is the single-argument application of function `id`, return its
// inner argument; else nullopt. Used for the f(f⁻¹(x)) = x simplifications.
[[nodiscard]] std::optional<Expr> arg_of(const Expr& arg, FunctionId id) {
    if (arg->type_id() != TypeId::Function) return std::nullopt;
    const auto& fn = static_cast<const Function&>(*arg);
    if (fn.function_id() == id && fn.args().size() == 1) return fn.args()[0];
    return std::nullopt;
}

// √e and 1/√e, for the cross-function inverse-composition algebraic forms.
[[nodiscard]] Expr sqrtp(const Expr& e) { return pow(e, rational(1, 2)); }
[[nodiscard]] Expr invsqrt(const Expr& e) { return pow(e, rational(-1, 2)); }
// x² + 1 and 1 − x².
[[nodiscard]] Expr x2_plus_1(const Expr& x) {
    return add(pow(x, integer(2)), S::One());
}
[[nodiscard]] Expr one_minus_x2(const Expr& x) {
    return add(S::One(), mul(S::NegativeOne(), pow(x, integer(2))));
}
// √(x−1)·√(x+1) — the form SymPy emits for the acosh compositions.
[[nodiscard]] Expr sqrt_x2_minus_1(const Expr& x) {
    return mul(sqrtp(add(x, S::NegativeOne())), sqrtp(add(x, S::One())));
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

// True if `e` still contains a bare application of inverse function `id` (the
// underlying factory left it unevaluated, possibly under a leading −1). Drives
// the acoth/asech/acsch factories — fold to the computed value when the inner
// inverse simplified, else keep the reciprocal node.
[[nodiscard]] bool is_bare_inverse(const Expr& e, FunctionId id) {
    auto is_it = [id](const Expr& f) {
        return f->type_id() == TypeId::Function
               && static_cast<const Function&>(*f).function_id() == id;
    };
    if (is_it(e)) return true;
    if (e->type_id() == TypeId::Mul) {
        for (const auto& f : e->args()) {
            if (is_it(f)) return true;
        }
    }
    return false;
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

// ----- Coth / Sech / Csch ----------------------------------------------------
// The reciprocal trio. coth/csch have poles at 0 (not finite everywhere), so we
// report only Real-for-real-arg (as Tan does), accepting the implicit non-pole
// caveat rather than claiming Finite.

Coth::Coth(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Coth);
}
Expr Coth::rebuild(std::vector<Expr> new_args) const { return coth(new_args[0]); }
std::optional<bool> Coth::ask(AssumptionKey k) const noexcept {
    if (k == AssumptionKey::Real && is_real(args_[0]) == true) return true;
    return std::nullopt;
}

Sech::Sech(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Sech);
}
Expr Sech::rebuild(std::vector<Expr> new_args) const { return sech(new_args[0]); }
std::optional<bool> Sech::ask(AssumptionKey k) const noexcept {
    if (k == AssumptionKey::Real && is_real(args_[0]) == true) return true;
    return std::nullopt;
}

Csch::Csch(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Csch);
}
Expr Csch::rebuild(std::vector<Expr> new_args) const { return csch(new_args[0]); }
std::optional<bool> Csch::ask(AssumptionKey k) const noexcept {
    if (k == AssumptionKey::Real && is_real(args_[0]) == true) return true;
    return std::nullopt;
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

// ----- Inverse reciprocal trio: acoth / asech / acsch ------------------------
// Domain-restricted (acoth: |x|>1, asech: 0<x≤1, acsch: x≠0), so ask stays
// agnostic.

Acoth::Acoth(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Acoth);
}
Expr Acoth::rebuild(std::vector<Expr> new_args) const { return acoth(new_args[0]); }
std::optional<bool> Acoth::ask(AssumptionKey /*k*/) const noexcept {
    return std::nullopt;
}

Asech::Asech(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Asech);
}
Expr Asech::rebuild(std::vector<Expr> new_args) const { return asech(new_args[0]); }
std::optional<bool> Asech::ask(AssumptionKey /*k*/) const noexcept {
    return std::nullopt;
}

Acsch::Acsch(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Acsch);
}
Expr Acsch::rebuild(std::vector<Expr> new_args) const { return acsch(new_args[0]); }
std::optional<bool> Acsch::ask(AssumptionKey /*k*/) const noexcept {
    return std::nullopt;
}

// ----- Factories -------------------------------------------------------------

Expr sinh(const Expr& arg) {
    if (arg == S::Zero()) return S::Zero();
    // sinh(±oo) = ±oo.
    if (arg->type_id() == TypeId::Infinity) return S::Infinity();
    if (arg->type_id() == TypeId::NegativeInfinity) return S::NegativeInfinity();
    if (arg->type_id() == TypeId::Float) {
        return unary_evalf(mpfr_sinh, arg);
    }
    // sinh(asinh(x)) = x.
    if (auto v = arg_of(arg, FunctionId::Asinh); v.has_value()) return *v;
    // sinh(acosh(x)) = √(x−1)·√(x+1).
    if (auto v = arg_of(arg, FunctionId::Acosh); v.has_value()) {
        return sqrt_x2_minus_1(*v);
    }
    // sinh(atanh(x)) = x / √(1 − x²).
    if (auto v = arg_of(arg, FunctionId::Atanh); v.has_value()) {
        return mul(*v, invsqrt(one_minus_x2(*v)));
    }
    if (auto pos = strip_neg(arg); pos.has_value()) {
        return mul(S::NegativeOne(), make<Sinh>(*pos));
    }
    return make<Sinh>(arg);
}

Expr cosh(const Expr& arg) {
    if (arg == S::Zero()) return S::One();
    // cosh(±oo) = oo (even, both directions grow without bound).
    if (arg->type_id() == TypeId::Infinity
        || arg->type_id() == TypeId::NegativeInfinity) {
        return S::Infinity();
    }
    if (arg->type_id() == TypeId::Float) {
        return unary_evalf(mpfr_cosh, arg);
    }
    // cosh(acosh(x)) = x.
    if (auto v = arg_of(arg, FunctionId::Acosh); v.has_value()) return *v;
    // cosh(asinh(x)) = √(x² + 1).
    if (auto v = arg_of(arg, FunctionId::Asinh); v.has_value()) {
        return sqrtp(x2_plus_1(*v));
    }
    // cosh(atanh(x)) = 1 / √(1 − x²).
    if (auto v = arg_of(arg, FunctionId::Atanh); v.has_value()) {
        return invsqrt(one_minus_x2(*v));
    }
    if (auto pos = strip_neg(arg); pos.has_value()) {
        return make<Cosh>(*pos);  // even
    }
    return make<Cosh>(arg);
}

Expr tanh(const Expr& arg) {
    if (arg == S::Zero()) return S::Zero();
    // tanh(±oo) = ±1.
    if (arg->type_id() == TypeId::Infinity) return S::One();
    if (arg->type_id() == TypeId::NegativeInfinity) return S::NegativeOne();
    if (arg->type_id() == TypeId::Float) {
        return unary_evalf(mpfr_tanh, arg);
    }
    // tanh(atanh(x)) = x.
    if (auto v = arg_of(arg, FunctionId::Atanh); v.has_value()) return *v;
    // tanh(asinh(x)) = x / √(x² + 1).
    if (auto v = arg_of(arg, FunctionId::Asinh); v.has_value()) {
        return mul(*v, invsqrt(x2_plus_1(*v)));
    }
    // tanh(acosh(x)) = √(x−1)·√(x+1) / x.
    if (auto v = arg_of(arg, FunctionId::Acosh); v.has_value()) {
        return mul(sqrt_x2_minus_1(*v), pow(*v, integer(-1)));
    }
    if (auto pos = strip_neg(arg); pos.has_value()) {
        return mul(S::NegativeOne(), make<Tanh>(*pos));
    }
    return make<Tanh>(arg);
}

Expr coth(const Expr& arg) {
    // coth = cosh/sinh, odd. coth(0)=zoo (pole), coth(±oo)=±1.
    if (arg == S::Zero()) return S::ComplexInfinity();
    if (arg->type_id() == TypeId::Infinity) return S::One();
    if (arg->type_id() == TypeId::NegativeInfinity) return S::NegativeOne();
    if (arg->type_id() == TypeId::Float) {
        return unary_evalf(mpfr_coth, arg);
    }
    // coth(atanh(x)) = 1/x.
    if (auto v = arg_of(arg, FunctionId::Atanh); v.has_value()) {
        return pow(*v, integer(-1));
    }
    if (auto pos = strip_neg(arg); pos.has_value()) {
        return mul(S::NegativeOne(), make<Coth>(*pos));
    }
    return make<Coth>(arg);
}

Expr sech(const Expr& arg) {
    // sech = 1/cosh, even. sech(0)=1, sech(±oo)=0.
    if (arg == S::Zero()) return S::One();
    if (arg->type_id() == TypeId::Infinity
        || arg->type_id() == TypeId::NegativeInfinity) {
        return S::Zero();
    }
    if (arg->type_id() == TypeId::Float) {
        return unary_evalf(mpfr_sech, arg);
    }
    // sech(acosh(x)) = 1/x.
    if (auto v = arg_of(arg, FunctionId::Acosh); v.has_value()) {
        return pow(*v, integer(-1));
    }
    if (auto pos = strip_neg(arg); pos.has_value()) {
        return make<Sech>(*pos);  // even
    }
    return make<Sech>(arg);
}

Expr csch(const Expr& arg) {
    // csch = 1/sinh, odd. csch(0)=zoo (pole), csch(±oo)=0.
    if (arg == S::Zero()) return S::ComplexInfinity();
    if (arg->type_id() == TypeId::Infinity
        || arg->type_id() == TypeId::NegativeInfinity) {
        return S::Zero();
    }
    if (arg->type_id() == TypeId::Float) {
        return unary_evalf(mpfr_csch, arg);
    }
    // csch(asinh(x)) = 1/x.
    if (auto v = arg_of(arg, FunctionId::Asinh); v.has_value()) {
        return pow(*v, integer(-1));
    }
    if (auto pos = strip_neg(arg); pos.has_value()) {
        return mul(S::NegativeOne(), make<Csch>(*pos));
    }
    return make<Csch>(arg);
}

Expr asinh(const Expr& arg) {
    if (arg == S::Zero()) return S::Zero();
    // asinh(±oo) = ±oo.
    if (arg->type_id() == TypeId::Infinity) return S::Infinity();
    if (arg->type_id() == TypeId::NegativeInfinity) return S::NegativeInfinity();
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
    // acosh(±oo) = oo.
    if (arg->type_id() == TypeId::Infinity
        || arg->type_id() == TypeId::NegativeInfinity) {
        return S::Infinity();
    }
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

// acoth(x) = atanh(1/x) [odd]. asech(x) = acosh(1/x) [asech(0)=oo]. acsch(x) =
// asinh(1/x) [acsch(0)=zoo, odd]. Each computes via the reciprocal-argument
// identity and keeps its own node when the inner inverse stays unevaluated.

Expr acoth(const Expr& arg) {
    if (auto pos = strip_neg(arg); pos.has_value()) {       // odd
        return mul(S::NegativeOne(), acoth(*pos));
    }
    Expr t = atanh(pow(arg, S::NegativeOne()));
    if (!is_bare_inverse(t, FunctionId::Atanh)) return t;
    return make<Acoth>(arg);
}

Expr asech(const Expr& arg) {
    if (arg == S::Zero()) return S::Infinity();             // asech(0) = oo
    Expr c = acosh(pow(arg, S::NegativeOne()));
    if (!is_bare_inverse(c, FunctionId::Acosh)) return c;
    return make<Asech>(arg);
}

Expr acsch(const Expr& arg) {
    if (arg == S::Zero()) return S::ComplexInfinity();      // acsch(0) = zoo
    if (auto pos = strip_neg(arg); pos.has_value()) {       // odd
        return mul(S::NegativeOne(), acsch(*pos));
    }
    Expr s = asinh(pow(arg, S::NegativeOne()));
    if (!is_bare_inverse(s, FunctionId::Asinh)) return s;
    return make<Acsch>(arg);
}

// ----- Derivatives ----------------------------------------------------------

Expr Sinh::diff_arg(std::size_t /*i*/) const {
    return cosh(args_[0]);
}
Expr Cosh::diff_arg(std::size_t /*i*/) const {
    return sinh(args_[0]);
}
Expr Tanh::diff_arg(std::size_t /*i*/) const {
    // 1 - tanh(x)^2
    auto t2 = pow(tanh(args_[0]), integer(2));
    return add(S::One(), mul(S::NegativeOne(), t2));
}
Expr Coth::diff_arg(std::size_t /*i*/) const {
    // d/dx coth(x) = -csch²(x).
    return mul(S::NegativeOne(), pow(csch(args_[0]), integer(2)));
}
Expr Sech::diff_arg(std::size_t /*i*/) const {
    // d/dx sech(x) = -sech(x)·tanh(x).
    return mul(S::NegativeOne(), mul(sech(args_[0]), tanh(args_[0])));
}
Expr Csch::diff_arg(std::size_t /*i*/) const {
    // d/dx csch(x) = -csch(x)·coth(x).
    return mul(S::NegativeOne(), mul(csch(args_[0]), coth(args_[0])));
}
Expr Asinh::diff_arg(std::size_t /*i*/) const {
    // 1 / sqrt(1 + x^2)
    auto one_plus_xsq = add(S::One(), pow(args_[0], integer(2)));
    return pow(one_plus_xsq, rational(-1, 2));
}
Expr Acosh::diff_arg(std::size_t /*i*/) const {
    // 1 / sqrt(x^2 - 1)
    auto xsq_minus_one = add(pow(args_[0], integer(2)), S::NegativeOne());
    return pow(xsq_minus_one, rational(-1, 2));
}
Expr Atanh::diff_arg(std::size_t /*i*/) const {
    // 1 / (1 - x^2)
    auto one_minus_xsq = add(S::One(), mul(S::NegativeOne(), pow(args_[0], integer(2))));
    return pow(one_minus_xsq, S::NegativeOne());
}
Expr Acoth::diff_arg(std::size_t /*i*/) const {
    // d/dx acoth(x) = 1/(1 - x²)  (same form as atanh').
    auto one_minus_xsq = add(S::One(), mul(S::NegativeOne(), pow(args_[0], integer(2))));
    return pow(one_minus_xsq, S::NegativeOne());
}
Expr Asech::diff_arg(std::size_t /*i*/) const {
    // d/dx asech(x) = -1/(x·√(1 - x²)).
    auto one_minus_xsq = add(S::One(), mul(S::NegativeOne(), pow(args_[0], integer(2))));
    auto denom = mul(args_[0], pow(one_minus_xsq, rational(1, 2)));
    return mul(S::NegativeOne(), pow(denom, S::NegativeOne()));
}
Expr Acsch::diff_arg(std::size_t /*i*/) const {
    // d/dx acsch(x) = -1/(x²·√(1 + 1/x²)).
    auto inner = add(S::One(), pow(args_[0], integer(-2)));
    auto denom = mul(pow(args_[0], integer(2)), pow(inner, rational(1, 2)));
    return mul(S::NegativeOne(), pow(denom, S::NegativeOne()));
}

}  // namespace sympp
