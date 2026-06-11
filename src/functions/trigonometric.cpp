#include <sympp/functions/trigonometric.hpp>

#include <algorithm>
#include <optional>
#include <utility>
#include <vector>

#include <mpfr.h>

#include <sympp/core/add.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/number_symbol.hpp>
#include <sympp/core/pow.hpp>
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

// If `arg` is the single-argument application of function `id`, return its
// inner argument; else nullopt. Used for the f(f⁻¹(x)) = x simplifications.
[[nodiscard]] std::optional<Expr> arg_of(const Expr& arg, FunctionId id) {
    if (arg->type_id() != TypeId::Function) return std::nullopt;
    const auto& fn = static_cast<const Function&>(*arg);
    if (fn.function_id() == id && fn.args().size() == 1) return fn.args()[0];
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

// ----- Exact values at rational multiples of π -------------------------------
// If `arg` is r·π for a rational r (including π itself), return r; else nullopt.
[[nodiscard]] std::optional<mpq_class> pi_coefficient(const Expr& arg) {
    if (arg == S::Pi()) return mpq_class(1);
    if (arg->type_id() != TypeId::Mul) return std::nullopt;
    bool has_pi = false;
    Expr coeff = S::One();
    for (const auto& f : arg->args()) {
        if (f == S::Pi()) {
            if (has_pi) return std::nullopt;  // π² etc. — not linear in π
            has_pi = true;
        } else if (is_number(f)) {
            coeff = mul(coeff, f);
        } else {
            return std::nullopt;  // a non-numeric, non-π factor
        }
    }
    if (!has_pi) return std::nullopt;
    if (coeff->type_id() == TypeId::Integer) {
        return mpq_class(static_cast<const Integer&>(*coeff).value());
    }
    if (coeff->type_id() == TypeId::Rational) {
        return static_cast<const Rational&>(*coeff).value();
    }
    return std::nullopt;
}

// cos(r·π) for a reference angle r ∈ [0, 1/2] with denominator in {1,2,3,4,6}.
// Returns nullopt for any other denominator (e.g. π/12 — a nested radical).
[[nodiscard]] std::optional<Expr> base_cos_pi(const mpq_class& r) {
    const mpz_class& num = r.get_num();
    const mpz_class& den = r.get_den();  // canonical: den > 0, gcd(num,den)=1
    if (num == 0) return S::One();                 // cos(0)        = 1
    if (num == 1 && den == 2) return S::Zero();     // cos(π/2)      = 0
    if (num == 1 && den == 3) return S::Half();      // cos(π/3)      = 1/2
    if (num == 1 && den == 4) {                      // cos(π/4)      = √2/2
        return mul(pow(integer(2), rational(1, 2)), S::Half());
    }
    if (num == 1 && den == 6) {                      // cos(π/6)      = √3/2
        return mul(pow(integer(3), rational(1, 2)), S::Half());
    }
    return std::nullopt;
}

// cos(r·π) for any rational r, via period (2) and reflection symmetry.
[[nodiscard]] std::optional<Expr> cos_pi(mpq_class r) {
    // Reduce modulo 2 into [0, 2): r -= 2·⌊r/2⌋.
    mpq_class half_r = r / 2;
    mpz_class k;
    mpz_fdiv_q(k.get_mpz_t(), half_r.get_num_mpz_t(), half_r.get_den_mpz_t());
    r -= 2 * mpq_class(k);
    // Fold [1, 2) → [0, 1] via cos(rπ) = cos((2−r)π).
    if (r.get_num() > r.get_den()) r = 2 - r;
    // Fold (1/2, 1] → [0, 1/2) with a sign flip: cos(rπ) = −cos((1−r)π).
    int sign = 1;
    if (2 * r.get_num() > r.get_den()) {
        sign = -1;
        r = 1 - r;
    }
    auto base = base_cos_pi(r);
    if (!base) return std::nullopt;
    return sign < 0 ? mul(S::NegativeOne(), *base) : *base;
}

// sin(r·π) = cos((1/2 − r)·π).
[[nodiscard]] std::optional<Expr> sin_pi(const mpq_class& r) {
    mpq_class half(1, 2);
    half.canonicalize();
    return cos_pi(half - r);
}

// tan(r·π) for a reference angle r ∈ [0, 1/2), denominator in {1,3,4,6}.
// Computed from a dedicated table (rather than sin/cos) for a clean result.
[[nodiscard]] std::optional<Expr> base_tan_pi(const mpq_class& r) {
    const mpz_class& num = r.get_num();
    const mpz_class& den = r.get_den();
    if (num == 0) return S::Zero();              // tan(0)   = 0
    if (num == 1 && den == 6) {                  // tan(π/6) = √3/3
        return mul(pow(integer(3), rational(1, 2)), rational(1, 3));
    }
    if (num == 1 && den == 4) return S::One();    // tan(π/4) = 1
    if (num == 1 && den == 3) {                  // tan(π/3) = √3
        return pow(integer(3), rational(1, 2));
    }
    return std::nullopt;
}

// tan(r·π) for any rational r, via period (1) and tan(π−x) = −tan(x). Returns
// nullopt at a pole (r ≡ 1/2 mod 1) or an out-of-table denominator.
[[nodiscard]] std::optional<Expr> tan_pi(mpq_class r) {
    mpz_class k;
    mpz_fdiv_q(k.get_mpz_t(), r.get_num_mpz_t(), r.get_den_mpz_t());
    r -= mpq_class(k);  // reduce modulo 1 into [0, 1)
    if (2 * r.get_num() == r.get_den()) return std::nullopt;  // pole at π/2
    int sign = 1;
    if (2 * r.get_num() > r.get_den()) {  // r > 1/2: tan(rπ) = −tan((1−r)π)
        sign = -1;
        r = 1 - r;
    }
    auto base = base_tan_pi(r);
    if (!base) return std::nullopt;
    return sign < 0 ? mul(S::NegativeOne(), *base) : *base;
}

// Split an additive argument `arg = rest + C·π`: returns the total rational
// π-coefficient C and the π-free remainder `rest`. Only fires for a true Add
// (a pure r·π Mul is left to pi_coefficient / the special-value tables).
[[nodiscard]] std::pair<mpq_class, Expr> split_pi_term(const Expr& arg) {
    if (arg->type_id() != TypeId::Add) return {mpq_class(0), arg};
    mpq_class total(0);
    std::vector<Expr> rest;
    for (const auto& term : arg->args()) {
        if (auto c = pi_coefficient(term); c.has_value()) {
            total += *c;
        } else {
            rest.push_back(term);
        }
    }
    if (rest.empty()) return {total, S::Zero()};
    return {total, add(std::move(rest))};
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

    // sin(asin(x)) = x.
    if (auto v = arg_of(arg, FunctionId::Asin); v.has_value()) return *v;

    // Exact value at a rational multiple of π (covers 0, π/6, π/4, π/3, π/2,
    // π and all their quadrant images).
    if (auto r = pi_coefficient(arg); r.has_value()) {
        if (auto v = sin_pi(*r); v.has_value()) return *v;
    }

    // Argument reduction by π multiples of the additive part:
    //   integer k:    sin(rest + k·π)     = (−1)^k·sin(rest)
    //   half-integer: sin(rest + (m/2)·π) = ±cos(rest)  (m odd, co-function)
    if (auto [c, rest] = split_pi_term(arg); c != 0) {
        if (c.get_den() == 1) {
            bool odd = mpz_odd_p(c.get_num_mpz_t());
            return odd ? mul(S::NegativeOne(), sin(rest)) : sin(rest);
        }
        if (c.get_den() == 2) {  // m odd; sin(mπ/2) = +1 (m≡1) / −1 (m≡3 mod 4)
            mpz_class m4 = ((c.get_num() % 4) + 4) % 4;
            return (m4 == 1) ? cos(rest) : mul(S::NegativeOne(), cos(rest));
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

    // cos(acos(x)) = x.
    if (auto v = arg_of(arg, FunctionId::Acos); v.has_value()) return *v;

    // Exact value at a rational multiple of π (covers 0, π/6, π/4, π/3, π/2,
    // π and all their quadrant images).
    if (auto r = pi_coefficient(arg); r.has_value()) {
        if (auto v = cos_pi(*r); v.has_value()) return *v;
    }

    // Argument reduction by π multiples of the additive part:
    //   integer k:    cos(rest + k·π)     = (−1)^k·cos(rest)
    //   half-integer: cos(rest + (m/2)·π) = ∓sin(rest)  (m odd, co-function)
    if (auto [c, rest] = split_pi_term(arg); c != 0) {
        if (c.get_den() == 1) {
            bool odd = mpz_odd_p(c.get_num_mpz_t());
            return odd ? mul(S::NegativeOne(), cos(rest)) : cos(rest);
        }
        if (c.get_den() == 2) {  // cos(rest + mπ/2) = −sin(mπ/2)·sin(rest)
            mpz_class m4 = ((c.get_num() % 4) + 4) % 4;
            return (m4 == 1) ? mul(S::NegativeOne(), sin(rest)) : sin(rest);
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

    // tan(atan(x)) = x.
    if (auto v = arg_of(arg, FunctionId::Atan); v.has_value()) return *v;

    // Exact value at a rational multiple of π. Poles (π/2 + kπ) and
    // out-of-table denominators are left unevaluated.
    if (auto r = pi_coefficient(arg); r.has_value()) {
        if (auto v = tan_pi(*r); v.has_value()) return *v;
    }

    // Period π: tan(rest + k·π) = tan(rest) for any integer k.
    if (auto [c, rest] = split_pi_term(arg); c.get_den() == 1 && c != 0) {
        return tan(rest);
    }

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

// √2/2 and √3/2 — built with the same constructors the forward trig table
// emits, so structural equality matches an incoming argument.
[[nodiscard]] Expr sqrt2_over_2() {
    return mul(pow(integer(2), rational(1, 2)), S::Half());
}
[[nodiscard]] Expr sqrt3_over_2() {
    return mul(pow(integer(3), rational(1, 2)), S::Half());
}

// asin of a non-negative special value → angle in [0, π/2], else nullopt.
[[nodiscard]] std::optional<Expr> asin_special(const Expr& arg) {
    if (arg == S::Zero()) return S::Zero();
    if (arg == S::Half()) return mul(rational(1, 6), S::Pi());
    if (arg == sqrt2_over_2()) return mul(rational(1, 4), S::Pi());
    if (arg == sqrt3_over_2()) return mul(rational(1, 3), S::Pi());
    if (arg == S::One()) return mul(S::Half(), S::Pi());
    return std::nullopt;
}

// atan of a non-negative special value → angle in [0, π/2), else nullopt.
[[nodiscard]] std::optional<Expr> atan_special(const Expr& arg) {
    if (arg == S::Zero()) return S::Zero();
    // √3/3 = tan(π/6)
    if (arg == mul(pow(integer(3), rational(1, 2)), rational(1, 3))) {
        return mul(rational(1, 6), S::Pi());
    }
    if (arg == S::One()) return mul(rational(1, 4), S::Pi());
    if (arg == pow(integer(3), rational(1, 2))) {  // √3 = tan(π/3)
        return mul(rational(1, 3), S::Pi());
    }
    return std::nullopt;
}

// True if `e` is an unevaluated Asin(...) or −Asin(...) — used so acos only
// folds via π/2 − asin when asin actually produced an exact angle.
[[nodiscard]] bool is_unevaluated_asin(const Expr& e) {
    auto is_asin = [](const Expr& f) {
        return f->type_id() == TypeId::Function
               && static_cast<const Function&>(*f).function_id()
                      == FunctionId::Asin;
    };
    if (is_asin(e)) return true;
    if (e->type_id() == TypeId::Mul) {
        for (const auto& f : e->args()) {
            if (is_asin(f)) return true;
        }
    }
    return false;
}

}  // namespace

Expr asin(const Expr& arg) {
    if (arg->type_id() == TypeId::Float) {
        return inv_trig_evalf(mpfr_asin, arg);
    }
    // Odd: asin(-x) = -asin(x). Recurse through the factory so the positive
    // part still gets the special-value table.
    if (auto pos = strip_neg(arg); pos.has_value()) {
        return mul(S::NegativeOne(), asin(*pos));
    }
    if (auto v = asin_special(arg); v.has_value()) return *v;
    return make<Asin>(arg);
}

Expr acos(const Expr& arg) {
    if (arg->type_id() == TypeId::Float) {
        return inv_trig_evalf(mpfr_acos, arg);
    }
    // acos(x) = π/2 − asin(x), but only adopt it when asin produced an exact
    // angle (otherwise keep acos(x) unevaluated, as SymPy does).
    Expr a = asin(arg);
    if (!is_unevaluated_asin(a)) {
        return add(mul(S::Half(), S::Pi()), mul(S::NegativeOne(), a));
    }
    return make<Acos>(arg);
}

Expr atan(const Expr& arg) {
    if (arg->type_id() == TypeId::Float) {
        return inv_trig_evalf(mpfr_atan, arg);
    }
    // Odd: atan(-x) = -atan(x).
    if (auto pos = strip_neg(arg); pos.has_value()) {
        return mul(S::NegativeOne(), atan(*pos));
    }
    if (auto v = atan_special(arg); v.has_value()) return *v;
    return make<Atan>(arg);
}

// ----- Derivatives ----------------------------------------------------------

Expr Sin::diff_arg(std::size_t /*i*/) const {
    return cos(args_[0]);
}
Expr Cos::diff_arg(std::size_t /*i*/) const {
    return mul(S::NegativeOne(), sin(args_[0]));
}
Expr Tan::diff_arg(std::size_t /*i*/) const {
    return add(S::One(), pow(tan(args_[0]), integer(2)));
}

Expr Asin::diff_arg(std::size_t /*i*/) const {
    auto one_minus_xsq = add(S::One(), mul(S::NegativeOne(), pow(args_[0], integer(2))));
    return pow(one_minus_xsq, rational(-1, 2));
}
Expr Acos::diff_arg(std::size_t /*i*/) const {
    auto one_minus_xsq = add(S::One(), mul(S::NegativeOne(), pow(args_[0], integer(2))));
    return mul(S::NegativeOne(), pow(one_minus_xsq, rational(-1, 2)));
}
Expr Atan::diff_arg(std::size_t /*i*/) const {
    auto one_plus_xsq = add(S::One(), pow(args_[0], integer(2)));
    return pow(one_plus_xsq, S::NegativeOne());
}
Expr Atan2::diff_arg(std::size_t i) const {
    auto y = args_[0];
    auto x = args_[1];
    auto denom = add(pow(x, integer(2)), pow(y, integer(2)));
    if (i == 0) return mul(x, pow(denom, S::NegativeOne()));
    return mul(mul(S::NegativeOne(), y), pow(denom, S::NegativeOne()));
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
