#include <sympp/functions/special.hpp>

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
#include <sympp/functions/exponential.hpp>
#include <sympp/functions/hyperbolic.hpp>
#include <sympp/functions/trigonometric.hpp>

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
    // erf(±oo) = ±1.
    if (arg->type_id() == TypeId::Infinity) return S::One();
    if (arg->type_id() == TypeId::NegativeInfinity) return S::NegativeOne();
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
    // erfc(oo) = 0, erfc(-oo) = 2.
    if (arg->type_id() == TypeId::Infinity) return S::Zero();
    if (arg->type_id() == TypeId::NegativeInfinity) return integer(2);
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

// ----- Derivatives -----------------------------------------------------------
// diff()'s chain rule supplies the arg' factor; diff_arg returns ∂f/∂(arg).

Expr Erf::diff_arg(std::size_t /*i*/) const {
    // d/dx erf(x) = 2·exp(−x²)/√π.
    const Expr& a = args_[0];
    return mul(integer(2),
               mul(exp(mul(S::NegativeOne(), pow(a, integer(2)))),
                   pow(S::Pi(), rational(-1, 2))));
}

Expr Erfc::diff_arg(std::size_t /*i*/) const {
    // d/dx erfc(x) = −2·exp(−x²)/√π.
    const Expr& a = args_[0];
    return mul(integer(-2),
               mul(exp(mul(S::NegativeOne(), pow(a, integer(2)))),
                   pow(S::Pi(), rational(-1, 2))));
}

Expr HeavisideFn::diff_arg(std::size_t /*i*/) const {
    // d/dx Heaviside(x) = DiracDelta(x).
    return dirac_delta(args_[0]);
}

// ----- Riemann zeta ----------------------------------------------------------

Zeta::Zeta(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Zeta);
}
Expr Zeta::rebuild(std::vector<Expr> new_args) const { return zeta(new_args[0]); }
std::optional<bool> Zeta::ask(AssumptionKey k) const noexcept {
    // zeta(s) is real for real s ≠ 1 (we accept the implicit pole caveat).
    if (k == AssumptionKey::Real && is_real(args_[0]) == true) return true;
    return std::nullopt;
}

Expr zeta(const Expr& s) {
    if (s->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*s);
        if (z.fits_long()) {
            long n = z.to_long();
            if (n == 1) return S::ComplexInfinity();           // pole
            if (n == 0) return rational(-1, 2);
            if (n > 0 && n % 2 == 0) {
                // Even positive: zeta(2m) = rₘ·π^(2m) (Basel and relatives).
                Expr coeff;
                switch (n) {
                    case 2:  coeff = rational(1, 6); break;
                    case 4:  coeff = rational(1, 90); break;
                    case 6:  coeff = rational(1, 945); break;
                    case 8:  coeff = rational(1, 9450); break;
                    case 10: coeff = rational(1, 93555); break;
                    case 12: coeff = rational(691, 638512875); break;
                    case 14: coeff = rational(2, 18243225); break;
                    default: break;
                }
                if (coeff) return mul(coeff, pow(S::Pi(), integer(n)));
            }
            if (n < 0) {
                // Negative integers: trivial zeros at even, rationals at odd
                // (zeta(-n) = -B_{n+1}/(n+1)).
                if (n % 2 == 0) return S::Zero();               // zeta(-2k) = 0
                switch (n) {
                    case -1: return rational(-1, 12);
                    case -3: return rational(1, 120);
                    case -5: return rational(-1, 252);
                    case -7: return rational(1, 240);
                    case -9: return rational(-1, 132);
                    default: break;                              // bigger → symbolic
                }
            }
        }
    }
    return make<Zeta>(s);
}

// ----- Lambert W (principal branch) ------------------------------------------

LambertWFn::LambertWFn(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::LambertW);
}
Expr LambertWFn::rebuild(std::vector<Expr> new_args) const {
    return lambertw(new_args[0]);
}
std::optional<bool> LambertWFn::ask(AssumptionKey k) const noexcept {
    // Real on [-1/e, ∞); we only assert it for a known non-negative argument.
    if (k == AssumptionKey::Real && is_nonnegative(args_[0]) == true) return true;
    return std::nullopt;
}
Expr LambertWFn::diff_arg(std::size_t /*i*/) const {
    // W'(x) = W(x) / (x·(1 + W(x))).
    const Expr& x = args_[0];
    Expr w = lambertw(x);
    return mul(w, pow(mul(x, add(S::One(), w)), S::NegativeOne()));
}

Expr lambertw(const Expr& arg) {
    if (arg == S::Zero()) return S::Zero();                       // W(0) = 0
    if (arg == S::E()) return S::One();                           // W(e) = 1
    // W(-1/e) = -1 (the branch point); -1/e is canonically -E^(-1).
    if (arg == mul(S::NegativeOne(), pow(S::E(), S::NegativeOne()))) {
        return S::NegativeOne();
    }
    if (arg->type_id() == TypeId::Infinity) return S::Infinity();  // W(oo) = oo
    return make<LambertWFn>(arg);
}

// ----- Exponential / sine / cosine integrals ---------------------------------

Ei::Ei(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Ei);
}
Expr Ei::rebuild(std::vector<Expr> new_args) const { return expint_ei(new_args[0]); }
std::optional<bool> Ei::ask(AssumptionKey k) const noexcept {
    if (k == AssumptionKey::Real && is_real(args_[0]) == true) return true;
    return std::nullopt;
}
Expr Ei::diff_arg(std::size_t /*i*/) const {
    // Ei'(x) = eˣ/x.
    return mul(exp(args_[0]), pow(args_[0], S::NegativeOne()));
}

Expr expint_ei(const Expr& arg) {
    if (arg == S::Zero()) return S::NegativeInfinity();   // Ei(0) = -oo
    if (arg->type_id() == TypeId::Infinity) return S::Infinity();
    return make<Ei>(arg);
}

Si::Si(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Si);
}
Expr Si::rebuild(std::vector<Expr> new_args) const { return sinint(new_args[0]); }
std::optional<bool> Si::ask(AssumptionKey k) const noexcept {
    if (k == AssumptionKey::Real && is_real(args_[0]) == true) return true;
    return std::nullopt;
}
Expr Si::diff_arg(std::size_t /*i*/) const {
    // Si'(x) = sin(x)/x.
    return mul(sin(args_[0]), pow(args_[0], S::NegativeOne()));
}

Expr sinint(const Expr& arg) {
    if (arg == S::Zero()) return S::Zero();               // Si(0) = 0
    if (arg->type_id() == TypeId::Infinity) return mul(S::Half(), S::Pi());
    if (arg->type_id() == TypeId::NegativeInfinity) {
        return mul(rational(-1, 2), S::Pi());
    }
    if (auto p = strip_neg(arg); p.has_value()) {          // odd: Si(-x) = -Si(x)
        return mul(S::NegativeOne(), sinint(*p));
    }
    return make<Si>(arg);
}

Ci::Ci(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Ci);
}
Expr Ci::rebuild(std::vector<Expr> new_args) const { return cosint(new_args[0]); }
std::optional<bool> Ci::ask(AssumptionKey k) const noexcept {
    if (k == AssumptionKey::Real && is_positive(args_[0]) == true) return true;
    return std::nullopt;
}
Expr Ci::diff_arg(std::size_t /*i*/) const {
    // Ci'(x) = cos(x)/x.
    return mul(cos(args_[0]), pow(args_[0], S::NegativeOne()));
}

Expr cosint(const Expr& arg) {
    if (arg->type_id() == TypeId::Infinity) return S::Zero();  // Ci(oo) = 0
    return make<Ci>(arg);
}

Shi::Shi(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Shi);
}
Expr Shi::rebuild(std::vector<Expr> new_args) const { return sinhint(new_args[0]); }
std::optional<bool> Shi::ask(AssumptionKey k) const noexcept {
    if (k == AssumptionKey::Real && is_real(args_[0]) == true) return true;
    return std::nullopt;
}
Expr Shi::diff_arg(std::size_t /*i*/) const {
    // Shi'(x) = sinh(x)/x.
    return mul(sinh(args_[0]), pow(args_[0], S::NegativeOne()));
}

Expr sinhint(const Expr& arg) {
    if (arg == S::Zero()) return S::Zero();                   // Shi(0) = 0
    if (arg->type_id() == TypeId::Infinity) return S::Infinity();
    if (arg->type_id() == TypeId::NegativeInfinity) {
        return S::NegativeInfinity();
    }
    if (auto p = strip_neg(arg); p.has_value()) {              // odd
        return mul(S::NegativeOne(), sinhint(*p));
    }
    return make<Shi>(arg);
}

Chi::Chi(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Chi);
}
Expr Chi::rebuild(std::vector<Expr> new_args) const { return coshint(new_args[0]); }
std::optional<bool> Chi::ask(AssumptionKey k) const noexcept {
    if (k == AssumptionKey::Real && is_positive(args_[0]) == true) return true;
    return std::nullopt;
}
Expr Chi::diff_arg(std::size_t /*i*/) const {
    // Chi'(x) = cosh(x)/x.
    return mul(cosh(args_[0]), pow(args_[0], S::NegativeOne()));
}

Expr coshint(const Expr& arg) {
    if (arg->type_id() == TypeId::Infinity) return S::Infinity();  // Chi(oo) = oo
    return make<Chi>(arg);
}

// ----- Polylogarithm ---------------------------------------------------------

Polylog::Polylog(Expr s, Expr z)
    : Function(std::vector<Expr>{std::move(s), std::move(z)}) {
    compute_hash(FunctionId::Polylog);
}
Expr Polylog::rebuild(std::vector<Expr> new_args) const {
    return polylog(new_args[0], new_args[1]);
}
std::optional<bool> Polylog::ask(AssumptionKey /*k*/) const noexcept {
    return std::nullopt;
}
Expr Polylog::diff_arg(std::size_t i) const {
    // d/dz Li_s(z) = Li_{s-1}(z)/z. The order-derivative (i==0) has no
    // elementary form; return 0 (consistent with the discrete-order convention).
    if (i == 1) {
        const Expr& s = args_[0];
        const Expr& z = args_[1];
        return mul(polylog(add(s, S::NegativeOne()), z),
                   pow(z, S::NegativeOne()));
    }
    return S::Zero();
}

Expr polylog(const Expr& s, const Expr& z) {
    if (z == S::Zero()) return S::Zero();          // Li_s(0) = 0
    if (z == S::One()) return zeta(s);             // Li_s(1) = ζ(s)
    // Li_2(-1) = -π²/12.
    if (s == integer(2) && z == S::NegativeOne()) {
        return mul(rational(-1, 12), pow(S::Pi(), integer(2)));
    }
    return make<Polylog>(s, z);
}

}  // namespace sympp
