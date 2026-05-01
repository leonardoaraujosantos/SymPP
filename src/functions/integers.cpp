#include <sympp/functions/integers.hpp>

#include <optional>
#include <utility>
#include <vector>

#include <gmpxx.h>
#include <mpfr.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
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

[[nodiscard]] std::optional<bool> int_ask(AssumptionKey k, const Expr& a) noexcept {
    switch (k) {
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
    return make<Ceiling>(arg);
}

}  // namespace sympp
