#include <sympp/functions/integers.hpp>

#include <optional>
#include <utility>
#include <vector>

#include <gmpxx.h>
#include <mpfr.h>

#include <sympp/core/add.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

// floor(p/q): mpz_fdiv_q rounds toward -inf; matches mathematical floor.
[[nodiscard]] Expr floor_rational(const Rational& q) {
    mpz_class r;
    mpz_fdiv_q(r.get_mpz_t(),
               mpq_numref(q.value().get_mpq_t()),
               mpq_denref(q.value().get_mpq_t()));
    return make<Integer>(std::move(r));
}

[[nodiscard]] Expr ceiling_rational(const Rational& q) {
    mpz_class r;
    mpz_cdiv_q(r.get_mpz_t(),
               mpq_numref(q.value().get_mpq_t()),
               mpq_denref(q.value().get_mpq_t()));
    return make<Integer>(std::move(r));
}

[[nodiscard]] Expr floor_float(const Float& f) {
    mpfr_t r;
    mpfr_init2(r, dps_to_prec(f.precision_dps()));
    mpfr_floor(r, f.value());
    auto out = make<Float>(static_cast<mpfr_srcptr>(r), f.precision_dps());
    mpfr_clear(r);
    return out;
}

[[nodiscard]] Expr ceiling_float(const Float& f) {
    mpfr_t r;
    mpfr_init2(r, dps_to_prec(f.precision_dps()));
    mpfr_ceil(r, f.value());
    auto out = make<Float>(static_cast<mpfr_srcptr>(r), f.precision_dps());
    mpfr_clear(r);
    return out;
}

// floor/ceiling of a constant real expression (pi, E, 2*pi, EulerGamma, …):
// evalf at high precision and round the resulting Float to an exact Integer.
// A boundary guard refuses to fold when the value sits within ~1e-40 of an
// integer (so a constant that is really an integer in disguise can't be
// mis-rounded); such cases stay symbolic.
[[nodiscard]] std::optional<Expr> constant_floor_ceiling(const Expr& arg,
                                                         bool is_ceiling) {
    if (!free_symbols(arg).empty()) return std::nullopt;
    // A real numeric constant evalfs to a Float; a complex or infinite value
    // does not, so the Float check both confirms reality and excludes I/oo.
    Expr v = evalf(arg, 60);
    if (v->type_id() != TypeId::Float) return std::nullopt;
    mpfr_srcptr val = static_cast<const Float&>(*v).value();
    if (mpfr_number_p(val) == 0) return std::nullopt;  // inf/nan guard

    // Distance to the nearest integer; bail if too close to call.
    mpfr_t nearest, diff;
    mpfr_init2(nearest, dps_to_prec(60));
    mpfr_init2(diff, dps_to_prec(60));
    mpfr_round(nearest, val);
    mpfr_sub(diff, val, nearest, MPFR_RNDN);
    mpfr_abs(diff, diff, MPFR_RNDN);
    const bool on_boundary = mpfr_cmp_d(diff, 1e-40) < 0;
    mpfr_clears(nearest, diff, static_cast<mpfr_ptr>(nullptr));
    if (on_boundary) return std::nullopt;

    mpz_class r;
    mpfr_get_z(r.get_mpz_t(), val, is_ceiling ? MPFR_RNDU : MPFR_RNDD);
    return make<Integer>(std::move(r));
}

[[nodiscard]] std::optional<bool> int_ask(AssumptionKey k, const Expr& a) noexcept {
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Real:
        case AssumptionKey::Integer:
        case AssumptionKey::Rational:
            if (is_real(a) == true) return true;
            return std::nullopt;
        case AssumptionKey::Finite:
            if (is_real(a) == true && is_finite(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

}  // namespace

// ----- Floor -----------------------------------------------------------------

Floor::Floor(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Floor);
}
Expr Floor::rebuild(std::vector<Expr> new_args) const { return floor(new_args[0]); }
std::optional<bool> Floor::ask(AssumptionKey k) const noexcept {
    return int_ask(k, args_[0]);
}

// ----- Ceiling ---------------------------------------------------------------

Ceiling::Ceiling(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Ceiling);
}
Expr Ceiling::rebuild(std::vector<Expr> new_args) const { return ceiling(new_args[0]); }
std::optional<bool> Ceiling::ask(AssumptionKey k) const noexcept {
    return int_ask(k, args_[0]);
}

// ----- Factories -------------------------------------------------------------

// If `arg` is a sum with one or more provably-integer terms and a non-integer
// remainder, returns (Σ integer terms, remaining sum). floor/ceiling are
// integer-shift invariant, so the integer part can be pulled out of the
// function. Returns nullopt when there's nothing to split.
[[nodiscard]] std::optional<std::pair<Expr, Expr>> pull_integer_shift(
    const Expr& arg) {
    if (arg->type_id() != TypeId::Add) return std::nullopt;
    std::vector<Expr> int_parts;
    std::vector<Expr> rest;
    for (const auto& term : arg->args()) {
        if (is_integer(term) == true) {
            int_parts.push_back(term);
        } else {
            rest.push_back(term);
        }
    }
    if (int_parts.empty() || rest.empty()) return std::nullopt;
    return std::make_pair(add(std::move(int_parts)), add(std::move(rest)));
}

Expr floor(const Expr& arg) {
    if (arg->type_id() == TypeId::Integer) return arg;
    if (arg->type_id() == TypeId::Rational) {
        return floor_rational(static_cast<const Rational&>(*arg));
    }
    if (arg->type_id() == TypeId::Float) {
        return floor_float(static_cast<const Float&>(*arg));
    }
    // Symbol with integer assumption -> identity.
    if (is_integer(arg) == true) return arg;
    // Integer-shift invariance: floor(n + x) = n + floor(x) when n is a sum of
    // provably-integer terms (e.g. floor(n + 1/2) = n).
    if (auto s = pull_integer_shift(arg); s.has_value()) {
        return add(s->first, floor(s->second));
    }
    // Constant real (pi, E, 2*pi, …): evaluate numerically.
    if (auto v = constant_floor_ceiling(arg, /*is_ceiling=*/false)) return *v;
    return make<Floor>(arg);
}

Expr ceiling(const Expr& arg) {
    if (arg->type_id() == TypeId::Integer) return arg;
    if (arg->type_id() == TypeId::Rational) {
        return ceiling_rational(static_cast<const Rational&>(*arg));
    }
    if (arg->type_id() == TypeId::Float) {
        return ceiling_float(static_cast<const Float&>(*arg));
    }
    if (is_integer(arg) == true) return arg;
    // Integer-shift invariance: ceiling(n + x) = n + ceiling(x).
    if (auto s = pull_integer_shift(arg); s.has_value()) {
        return add(s->first, ceiling(s->second));
    }
    // Constant real (pi, E, 2*pi, …): evaluate numerically.
    if (auto v = constant_floor_ceiling(arg, /*is_ceiling=*/true)) return *v;
    return make<Ceiling>(arg);
}

// ----- Frac (fractional part) ------------------------------------------------

Frac::Frac(Expr arg) : Function(std::vector<Expr>{std::move(arg)}) {
    compute_hash(FunctionId::Frac);
}
Expr Frac::rebuild(std::vector<Expr> new_args) const { return frac(new_args[0]); }
std::optional<bool> Frac::ask(AssumptionKey k) const noexcept {
    const auto& a = args_[0];
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Real:
        case AssumptionKey::Nonnegative:
            if (is_real(a) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Expr frac(const Expr& arg) {
    // frac(x) = x − floor(x), always in [0, 1). Reuse floor's numeric/constant
    // folding: when floor evaluates, return the difference; otherwise keep Frac.
    Expr fl = floor(arg);
    const bool fl_unevaluated =
        fl->type_id() == TypeId::Function
        && static_cast<const Function&>(*fl).function_id() == FunctionId::Floor;
    if (!fl_unevaluated) {
        return add(arg, mul(S::NegativeOne(), fl));
    }
    return make<Frac>(arg);
}

// ----- Mod (floored modulo) --------------------------------------------------

Mod::Mod(Expr p, Expr q)
    : Function(std::vector<Expr>{std::move(p), std::move(q)}) {
    compute_hash(FunctionId::Mod);
}
Expr Mod::rebuild(std::vector<Expr> new_args) const {
    return mod(new_args[0], new_args[1]);
}
std::optional<bool> Mod::ask(AssumptionKey k) const noexcept {
    const auto& p = args_[0];
    const auto& q = args_[1];
    switch (k) {
        case AssumptionKey::Complex:
        case AssumptionKey::Imaginary:
            return std::nullopt;  // derived by the generic ask() layer

        case AssumptionKey::Integer:
        case AssumptionKey::Real:
        case AssumptionKey::Rational:
            if (is_integer(p) == true && is_integer(q) == true) return true;
            return std::nullopt;
        default:
            return std::nullopt;
    }
}

Expr mod(const Expr& p, const Expr& q) {
    // q == 0 is undefined (SymPy raises); keep it unevaluated rather than throw.
    if (q == S::Zero()) return make<Mod>(p, q);
    // Structural identities valid for any p, q: mod(0,q)=0, mod(x,x)=0.
    if (p == S::Zero()) return S::Zero();
    if (p == q) return S::Zero();
    // Mod(integer, 1) = 0 (an integer leaves no remainder modulo 1).
    if (q == S::One() && is_integer(p) == true) return S::Zero();

    auto to_mpq = [](const Expr& e) -> std::optional<mpq_class> {
        if (e->type_id() == TypeId::Integer) {
            return mpq_class(static_cast<const Integer&>(*e).value());
        }
        if (e->type_id() == TypeId::Rational) {
            return static_cast<const Rational&>(*e).value();
        }
        return std::nullopt;
    };
    auto mp = to_mpq(p);
    auto mq = to_mpq(q);
    if (mp && mq) {
        // Floored modulo: r = p − q·⌊p/q⌋ (result takes the sign of q).
        mpq_class ratio = *mp / *mq;
        mpz_class fl;
        mpz_fdiv_q(fl.get_mpz_t(), mpq_numref(ratio.get_mpq_t()),
                   mpq_denref(ratio.get_mpq_t()));
        mpq_class r = *mp - *mq * mpq_class(fl);
        r.canonicalize();
        return rational(mpz_class(mpq_numref(r.get_mpq_t())),
                        mpz_class(mpq_denref(r.get_mpq_t())));
    }
    return make<Mod>(p, q);
}

}  // namespace sympp
