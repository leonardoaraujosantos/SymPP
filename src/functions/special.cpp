#include <sympp/functions/special.hpp>

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
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

[[nodiscard]] std::optional<Expr> strip_neg(const Expr& e) {
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

[[nodiscard]] Expr float_unary_op(int (*op)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t),
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

// ============================================================================
// erf
// ============================================================================

Erf::Erf(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Erf);
}
Expr Erf::rebuild(std::vector<Expr> new_args) const { return erf(new_args[0]); }
std::optional<bool> Erf::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Real:
        case AssumptionKey::Finite:
            if (is_real(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Positive:
            if (is_positive(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Negative:
            if (is_negative(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Zero:
            if (is_zero(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Expr erf(const Expr& arg) {
    if (arg == S::Zero()) return S::Zero();
    if (arg->type_id() == TypeId::Float) {
        return float_unary_op(mpfr_erf, arg);
    }
    // Odd: erf(-x) = -erf(x)
    if (auto p = strip_neg(arg); p.has_value()) {
        return mul(S::NegativeOne(), make<Erf>(*p));
    }
    return make<Erf>(arg);
}

// ============================================================================
// erfc
// ============================================================================

Erfc::Erfc(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Erfc);
}
Expr Erfc::rebuild(std::vector<Expr> new_args) const { return erfc(new_args[0]); }
std::optional<bool> Erfc::ask(AssumptionKey k) const noexcept {
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

Expr erfc(const Expr& arg) {
    if (arg == S::Zero()) return S::One();
    if (arg->type_id() == TypeId::Float) {
        return float_unary_op(mpfr_erfc, arg);
    }
    return make<Erfc>(arg);
}

// ============================================================================
// Heaviside
// ============================================================================

HeavisideFn::HeavisideFn(Expr arg)
    : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Heaviside);
}
Expr HeavisideFn::rebuild(std::vector<Expr> new_args) const {
    return heaviside(new_args[0]);
}
std::optional<bool> HeavisideFn::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Real:
        case AssumptionKey::Finite:
            if (is_real(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Nonnegative:
            return true;
        default:
            return std::nullopt;
    }
}

Expr heaviside(const Expr& arg) {
    // SymPy default H(0) = 1/2
    if (arg == S::Zero()) return S::Half();
    if (is_positive(arg) == true) return S::One();
    if (is_negative(arg) == true) return S::Zero();
    return make<HeavisideFn>(arg);
}

// ============================================================================
// DiracDelta
// ============================================================================

DiracDeltaFn::DiracDeltaFn(Expr arg)
    : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::DiracDelta);
}
Expr DiracDeltaFn::rebuild(std::vector<Expr> new_args) const {
    return dirac_delta(new_args[0]);
}
std::optional<bool> DiracDeltaFn::ask(AssumptionKey k) const noexcept {
    // DiracDelta is a distribution; the only safe predicate at the function-
    // value level is nonneg (it's never negative).
    switch (k) {
        case AssumptionKey::Nonnegative:
            return true;
        default:
            return std::nullopt;
    }
}

Expr dirac_delta(const Expr& arg) {
    if (is_nonzero(arg) == true) return S::Zero();
    return make<DiracDeltaFn>(arg);
}

}  // namespace sympp
