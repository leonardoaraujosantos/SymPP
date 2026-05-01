#include <sympp/functions/trigonometric.hpp>

#include <algorithm>
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
#include <sympp/core/number_symbol.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

// Detect a leading minus sign on a Mul: i.e. -1 * rest. Returns the
// stripped tail if so, std::nullopt otherwise.
[[nodiscard]] std::optional<Expr> strip_neg(const Expr& e) {
    // Negative numeric literal: just check sign.
    if (is_number(e)) {
        const auto& n = static_cast<const Number&>(*e);
        if (n.is_negative()) {
            // Build the absolute via numeric mul by -1.
            // Easier path: use the Mul factory which handles this.
            return mul(S::NegativeOne(), e);
        }
        return std::nullopt;
    }
    // Mul with leading -1.
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

// Numeric evaluator: evalf the argument and apply the MPFR trig op.
// Used when the user passes an explicit Float, or wants numeric output.
[[nodiscard]] Expr trig_evalf(int (*op)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t),
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

// ----- Sin -------------------------------------------------------------------

Sin::Sin(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Sin);
}

Expr Sin::rebuild(std::vector<Expr> new_args) const {
    return sin(new_args[0]);
}

std::optional<bool> Sin::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Real:
            if (is_real(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Finite:
            if (is_real(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

// ----- Cos -------------------------------------------------------------------

Cos::Cos(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Cos);
}

Expr Cos::rebuild(std::vector<Expr> new_args) const {
    return cos(new_args[0]);
}

std::optional<bool> Cos::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Real:
            if (is_real(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Finite:
            if (is_real(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

// ----- Tan -------------------------------------------------------------------

Tan::Tan(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Tan);
}

Expr Tan::rebuild(std::vector<Expr> new_args) const {
    return tan(new_args[0]);
}

std::optional<bool> Tan::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Real:
            // tan is real wherever cos(arg) ≠ 0; we don't have that info in
            // general so only return true when the arg is itself real (and
            // accept the implicit cos(arg) ≠ 0 caveat as in SymPy).
            if (is_real(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

// ----- Factories -------------------------------------------------------------

Expr sin(const Expr& arg) {
    // sin(0) = 0
    if (arg == S::Zero()) return S::Zero();

    // Numeric Float -> evalf via mpfr_sin
    if (arg->type_id() == TypeId::Float) {
        return trig_evalf(mpfr_sin, arg);
    }

    // sin(pi) = 0
    if (arg == S::Pi()) return S::Zero();

    // sin(pi/2) = 1
    if (arg->type_id() == TypeId::Mul) {
        const auto& a = arg->args();
        if (a.size() == 2 && a[0] == S::Half() && a[1] == S::Pi()) {
            return S::One();
        }
    }

    // Odd: sin(-x) = -sin(x)
    if (auto pos = strip_neg(arg); pos.has_value()) {
        return mul(S::NegativeOne(), make<Sin>(*pos));
    }

    return make<Sin>(arg);
}

Expr cos(const Expr& arg) {
    // cos(0) = 1
    if (arg == S::Zero()) return S::One();

    // Numeric Float -> evalf via mpfr_cos
    if (arg->type_id() == TypeId::Float) {
        return trig_evalf(mpfr_cos, arg);
    }

    // cos(pi) = -1
    if (arg == S::Pi()) return S::NegativeOne();

    // cos(pi/2) = 0
    if (arg->type_id() == TypeId::Mul) {
        const auto& a = arg->args();
        if (a.size() == 2 && a[0] == S::Half() && a[1] == S::Pi()) {
            return S::Zero();
        }
    }

    // Even: cos(-x) = cos(x)
    if (auto pos = strip_neg(arg); pos.has_value()) {
        return make<Cos>(*pos);
    }

    return make<Cos>(arg);
}

Expr tan(const Expr& arg) {
    // tan(0) = 0
    if (arg == S::Zero()) return S::Zero();

    // Numeric Float -> evalf via mpfr_tan
    if (arg->type_id() == TypeId::Float) {
        return trig_evalf(mpfr_tan, arg);
    }

    // tan(pi) = 0
    if (arg == S::Pi()) return S::Zero();

    // Odd: tan(-x) = -tan(x)
    if (auto pos = strip_neg(arg); pos.has_value()) {
        return mul(S::NegativeOne(), make<Tan>(*pos));
    }

    return make<Tan>(arg);
}

// ============================================================================
// Inverse trigonometric
// ============================================================================

Asin::Asin(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Asin);
}
Expr Asin::rebuild(std::vector<Expr> new_args) const { return asin(new_args[0]); }
std::optional<bool> Asin::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Real:
            // Real iff arg is real and |arg| <= 1. We don't track |arg|<=1
            // precisely yet; conservatively flag real only when arg is itself
            // real (the principal-branch convention SymPy uses).
            if (is_real(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Acos::Acos(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Acos);
}
Expr Acos::rebuild(std::vector<Expr> new_args) const { return acos(new_args[0]); }
std::optional<bool> Acos::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Real:
            if (is_real(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Nonnegative:
            // acos image is [0, pi], always nonnegative for real inputs.
            if (is_real(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Atan::Atan(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Atan);
}
Expr Atan::rebuild(std::vector<Expr> new_args) const { return atan(new_args[0]); }
std::optional<bool> Atan::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Real:
        case AssumptionKey::Finite:
            if (is_real(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Atan2::Atan2(Expr y, Expr x)
    : Function(std::vector<Expr>{std::move(y), std::move(x)}) {
    compute_hash(FunctionId::Atan2);
}
Expr Atan2::rebuild(std::vector<Expr> new_args) const {
    return atan2(new_args[0], new_args[1]);
}
std::optional<bool> Atan2::ask(AssumptionKey k) const noexcept {
    const auto& y = args_[0];
    const auto& x = args_[1];
    switch (k) {
        case AssumptionKey::Real:
            if (is_real(y) == true && is_real(x) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

// ----- asin / acos / atan factories -----------------------------------------

namespace {

// Numeric Float arg evaluator for unary inverse trig. mpfr's asin/acos/atan
// silently return NaN on out-of-domain inputs; we forward that behavior.
[[nodiscard]] Expr inv_trig_evalf(int (*op)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t),
                                  const Expr& arg) {
    return trig_evalf(op, arg);
}

}  // namespace

Expr asin(const Expr& arg) {
    if (arg == S::Zero()) return S::Zero();
    if (arg == S::One()) return mul(S::Half(), S::Pi());
    if (arg == S::NegativeOne()) return mul(S::NegativeOne(), mul(S::Half(), S::Pi()));

    if (arg->type_id() == TypeId::Float) {
        return inv_trig_evalf(mpfr_asin, arg);
    }

    if (auto pos = strip_neg(arg); pos.has_value()) {
        return mul(S::NegativeOne(), make<Asin>(*pos));
    }

    return make<Asin>(arg);
}

Expr acos(const Expr& arg) {
    if (arg == S::Zero()) return mul(S::Half(), S::Pi());
    if (arg == S::One()) return S::Zero();
    if (arg == S::NegativeOne()) return S::Pi();

    if (arg->type_id() == TypeId::Float) {
        return inv_trig_evalf(mpfr_acos, arg);
    }

    return make<Acos>(arg);
}

Expr atan(const Expr& arg) {
    if (arg == S::Zero()) return S::Zero();
    if (arg == S::One()) return mul(rational(1, 4), S::Pi());
    if (arg == S::NegativeOne()) return mul(rational(-1, 4), S::Pi());

    if (arg->type_id() == TypeId::Float) {
        return inv_trig_evalf(mpfr_atan, arg);
    }

    if (auto pos = strip_neg(arg); pos.has_value()) {
        return mul(S::NegativeOne(), make<Atan>(*pos));
    }

    return make<Atan>(arg);
}

Expr atan2(const Expr& y, const Expr& x) {
    // Both numeric Floats -> mpfr_atan2.
    if (y->type_id() == TypeId::Float && x->type_id() == TypeId::Float) {
        const auto& fy = static_cast<const Float&>(*y);
        const auto& fx = static_cast<const Float&>(*x);
        int dps = std::max(fy.precision_dps(), fx.precision_dps());
        mpfr_t r;
        mpfr_init2(r, dps_to_prec(dps));
        mpfr_atan2(r, fy.value(), fx.value(), MPFR_RNDN);
        auto out = make<Float>(static_cast<mpfr_srcptr>(r), dps);
        mpfr_clear(r);
        return out;
    }

    // Special exact reductions.
    if (y == S::Zero() && x == S::Zero()) {
        // Conventionally undefined; SymPy keeps it unevaluated.
        return make<Atan2>(y, x);
    }
    if (y == S::Zero()) {
        // atan2(0, x) = 0 for x>0, pi for x<0
        if (is_positive(x) == true) return S::Zero();
        if (is_negative(x) == true) return S::Pi();
    }
    if (x == S::Zero()) {
        // atan2(y, 0) = pi/2 for y>0, -pi/2 for y<0
        if (is_positive(y) == true) return mul(S::Half(), S::Pi());
        if (is_negative(y) == true) {
            return mul(S::NegativeOne(), mul(S::Half(), S::Pi()));
        }
    }

    return make<Atan2>(y, x);
}

}  // namespace sympp
