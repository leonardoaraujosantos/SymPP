#include <sympp/functions/trigonometric.hpp>

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

}  // namespace sympp
