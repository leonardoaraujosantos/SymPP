#include <sympp/core/number_arith.hpp>

#include <algorithm>
#include <optional>
#include <utility>

#include <gmpxx.h>
#include <mpfr.h>

#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

[[nodiscard]] bool is_int(const Number& n) noexcept { return n.type_id() == TypeId::Integer; }
[[nodiscard]] bool is_rat(const Number& n) noexcept { return n.type_id() == TypeId::Rational; }
[[nodiscard]] bool is_flt(const Number& n) noexcept { return n.type_id() == TypeId::Float; }

[[nodiscard]] const Integer& as_int(const Number& n) noexcept {
    return static_cast<const Integer&>(n);
}
[[nodiscard]] const Rational& as_rat(const Number& n) noexcept {
    return static_cast<const Rational&>(n);
}
[[nodiscard]] const Float& as_flt(const Number& n) noexcept {
    return static_cast<const Float&>(n);
}

// Largest exact integer/rational power we are willing to materialize, expressed
// as a result bit-length cap. base^ne needs roughly ne·bitlen(base) bits; GMP
// aborts with a hardware SIGFPE (via __gmp_overflow_in_mpz) when an mpz exceeds
// its internal size limit, so an unbounded exact power on input like
// (1000001/1000000)^(10^12) crashes the whole process. We cap the materialized
// size and defer (return nullopt) when exceeded, leaving the power symbolic /
// available to the float path instead of crashing. 2^27 bits ≈ 16 MB ≈ 40M
// decimal digits — far beyond any sane exact result, far below GMP's limit.
constexpr unsigned long kMaxPowResultBits = 1UL << 27;

// True when materializing base^|ne| as an exact integer would stay within the
// size cap. `base_bits` is the bit-length of the integer being raised.
[[nodiscard]] bool pow_result_within_cap(size_t base_bits, unsigned long ne) noexcept {
    if (base_bits == 0 || ne == 0) return true;  // 0/±1 base or exponent 0
    // Avoid overflow in the product: compare ne against cap / base_bits.
    return ne <= kMaxPowResultBits / base_bits;
}

// Normalize a binary numeric op so that if either operand is a Float, both are.
// Returns the working precision (max of the two when both Float; otherwise the
// Float's precision when only one is Float).
[[nodiscard]] int float_dps(const Number& a, const Number& b) noexcept {
    int dps = kDefaultDps;
    if (is_flt(a)) dps = std::max(dps, as_flt(a).precision_dps());
    if (is_flt(b)) dps = std::max(dps, as_flt(b).precision_dps());
    return dps;
}

}  // namespace

bool is_number(const Expr& e) noexcept {
    if (!e) return false;
    auto t = e->type_id();
    return t == TypeId::Integer || t == TypeId::Rational || t == TypeId::Float;
}

std::optional<Expr> number_add(const Number& a, const Number& b) {
    // --- exact arithmetic ---
    if (is_int(a) && is_int(b)) {
        mpz_class r = as_int(a).value() + as_int(b).value();
        return make<Integer>(std::move(r));
    }
    if ((is_int(a) || is_rat(a)) && (is_int(b) || is_rat(b))) {
        // Promote to Rational. mpq accepts mixed.
        mpq_class qa, qb;
        if (is_int(a)) qa = as_int(a).value();
        else qa = as_rat(a).value();
        if (is_int(b)) qb = as_int(b).value();
        else qb = as_rat(b).value();
        mpq_class r = qa + qb;
        r.canonicalize();
        if (r.get_den() == 1) {
            return make<Integer>(r.get_num());
        }
        return make<Rational>(std::move(r));
    }

    // --- floating point (Float involved) ---
    if (is_flt(a) || is_flt(b)) {
        int dps = float_dps(a, b);
        mpfr_t ax, bx, rx;
        mpfr_init2(ax, dps_to_prec(dps));
        mpfr_init2(bx, dps_to_prec(dps));
        mpfr_init2(rx, dps_to_prec(dps));

        if (is_flt(a)) mpfr_set(ax, as_flt(a).value(), MPFR_RNDN);
        else if (is_int(a)) mpfr_set_z(ax, as_int(a).value().get_mpz_t(), MPFR_RNDN);
        else /* rational */ mpfr_set_q(ax, as_rat(a).value().get_mpq_t(), MPFR_RNDN);

        if (is_flt(b)) mpfr_set(bx, as_flt(b).value(), MPFR_RNDN);
        else if (is_int(b)) mpfr_set_z(bx, as_int(b).value().get_mpz_t(), MPFR_RNDN);
        else mpfr_set_q(bx, as_rat(b).value().get_mpq_t(), MPFR_RNDN);

        mpfr_add(rx, ax, bx, MPFR_RNDN);
        auto result = make<Float>(static_cast<mpfr_srcptr>(rx), dps);
        mpfr_clear(ax);
        mpfr_clear(bx);
        mpfr_clear(rx);
        return result;
    }

    return std::nullopt;
}

std::optional<Expr> number_mul(const Number& a, const Number& b) {
    if (is_int(a) && is_int(b)) {
        mpz_class r = as_int(a).value() * as_int(b).value();
        return make<Integer>(std::move(r));
    }
    if ((is_int(a) || is_rat(a)) && (is_int(b) || is_rat(b))) {
        mpq_class qa, qb;
        if (is_int(a)) qa = as_int(a).value();
        else qa = as_rat(a).value();
        if (is_int(b)) qb = as_int(b).value();
        else qb = as_rat(b).value();
        mpq_class r = qa * qb;
        r.canonicalize();
        if (r.get_den() == 1) {
            return make<Integer>(r.get_num());
        }
        return make<Rational>(std::move(r));
    }

    if (is_flt(a) || is_flt(b)) {
        int dps = float_dps(a, b);
        mpfr_t ax, bx, rx;
        mpfr_init2(ax, dps_to_prec(dps));
        mpfr_init2(bx, dps_to_prec(dps));
        mpfr_init2(rx, dps_to_prec(dps));

        if (is_flt(a)) mpfr_set(ax, as_flt(a).value(), MPFR_RNDN);
        else if (is_int(a)) mpfr_set_z(ax, as_int(a).value().get_mpz_t(), MPFR_RNDN);
        else mpfr_set_q(ax, as_rat(a).value().get_mpq_t(), MPFR_RNDN);

        if (is_flt(b)) mpfr_set(bx, as_flt(b).value(), MPFR_RNDN);
        else if (is_int(b)) mpfr_set_z(bx, as_int(b).value().get_mpz_t(), MPFR_RNDN);
        else mpfr_set_q(bx, as_rat(b).value().get_mpq_t(), MPFR_RNDN);

        mpfr_mul(rx, ax, bx, MPFR_RNDN);
        auto result = make<Float>(static_cast<mpfr_srcptr>(rx), dps);
        mpfr_clear(ax);
        mpfr_clear(bx);
        mpfr_clear(rx);
        return result;
    }

    return std::nullopt;
}

std::optional<Expr> number_pow(const Number& base, const Number& exp) {
    // ---- universal identities ----
    if (exp.is_zero()) {
        // x^0 = 1 by SymPy convention (including 0^0).
        return make<Integer>(1);
    }
    if (exp.is_one()) {
        // x^1 = x. Caller must give back base shared_ptr, which we don't have
        // here — return nullopt and let the caller short-circuit.
        return std::nullopt;
    }
    if (base.is_one()) {
        return make<Integer>(1);
    }
    if (base.is_zero() && exp.is_positive()) {
        return make<Integer>(0);
    }

    // ---- Integer base, Integer exponent ----
    if (is_int(base) && is_int(exp)) {
        const auto& z = as_int(base);
        const auto& n = as_int(exp);
        if (!n.fits_long()) return std::nullopt;  // huge exponent — defer
        long ne = n.to_long();
        size_t base_bits = mpz_sizeinbase(z.value().get_mpz_t(), 2);
        unsigned long mag = static_cast<unsigned long>(ne >= 0 ? ne : -ne);
        if (!pow_result_within_cap(base_bits, mag)) return std::nullopt;  // too big — defer
        if (ne >= 0) {
            mpz_class r;
            mpz_pow_ui(r.get_mpz_t(), z.value().get_mpz_t(), static_cast<unsigned long>(ne));
            return make<Integer>(std::move(r));
        }
        // negative integer exponent: 1 / base^|ne|
        if (z.is_zero()) return std::nullopt;  // 0^negative — ComplexInfinity (Phase 1e+)
        mpz_class denom;
        mpz_pow_ui(denom.get_mpz_t(), z.value().get_mpz_t(), static_cast<unsigned long>(-ne));
        // Build Rational(1, denom). Sign normalization handled by canonicalize.
        mpq_class q;
        mpz_set_si(mpq_numref(q.get_mpq_t()), 1);
        mpz_set(mpq_denref(q.get_mpq_t()), denom.get_mpz_t());
        q.canonicalize();
        if (q.get_den() == 1) return make<Integer>(q.get_num());
        return make<Rational>(std::move(q));
    }

    // ---- Rational base, Integer exponent ----
    if (is_rat(base) && is_int(exp)) {
        const auto& q = as_rat(base);
        const auto& n = as_int(exp);
        if (!n.fits_long()) return std::nullopt;
        long ne = n.to_long();
        unsigned long mag = static_cast<unsigned long>(ne >= 0 ? ne : -ne);
        // Both numerator and denominator are raised to |ne|; bound by the larger.
        size_t num_bits = mpz_sizeinbase(mpq_numref(q.value().get_mpq_t()), 2);
        size_t den_bits = mpz_sizeinbase(mpq_denref(q.value().get_mpq_t()), 2);
        if (!pow_result_within_cap(std::max(num_bits, den_bits), mag)) {
            return std::nullopt;  // too big — defer
        }
        mpq_class r;
        if (ne >= 0) {
            mpz_pow_ui(mpq_numref(r.get_mpq_t()), mpq_numref(q.value().get_mpq_t()),
                       static_cast<unsigned long>(ne));
            mpz_pow_ui(mpq_denref(r.get_mpq_t()), mpq_denref(q.value().get_mpq_t()),
                       static_cast<unsigned long>(ne));
        } else {
            // (p/q)^-n = (q/p)^n
            mpz_pow_ui(mpq_numref(r.get_mpq_t()), mpq_denref(q.value().get_mpq_t()),
                       static_cast<unsigned long>(-ne));
            mpz_pow_ui(mpq_denref(r.get_mpq_t()), mpq_numref(q.value().get_mpq_t()),
                       static_cast<unsigned long>(-ne));
        }
        r.canonicalize();
        if (r.get_den() == 1) return make<Integer>(r.get_num());
        return make<Rational>(std::move(r));
    }

    // ---- Float involved ----
    if (is_flt(base) && is_int(exp)) {
        const auto& f = as_flt(base);
        const auto& n = as_int(exp);
        if (!n.fits_long()) return std::nullopt;
        mpfr_t r;
        mpfr_init2(r, dps_to_prec(f.precision_dps()));
        mpfr_pow_si(r, f.value(), n.to_long(), MPFR_RNDN);
        auto result = make<Float>(static_cast<mpfr_srcptr>(r), f.precision_dps());
        mpfr_clear(r);
        return result;
    }
    if ((is_flt(base) || is_flt(exp)) && (is_int(exp) || is_rat(exp) || is_flt(exp))) {
        int dps = float_dps(base, exp);
        mpfr_t bx, ex, rx;
        mpfr_init2(bx, dps_to_prec(dps));
        mpfr_init2(ex, dps_to_prec(dps));
        mpfr_init2(rx, dps_to_prec(dps));

        if (is_flt(base)) mpfr_set(bx, as_flt(base).value(), MPFR_RNDN);
        else if (is_int(base)) mpfr_set_z(bx, as_int(base).value().get_mpz_t(), MPFR_RNDN);
        else mpfr_set_q(bx, as_rat(base).value().get_mpq_t(), MPFR_RNDN);

        if (is_flt(exp)) mpfr_set(ex, as_flt(exp).value(), MPFR_RNDN);
        else if (is_int(exp)) mpfr_set_z(ex, as_int(exp).value().get_mpz_t(), MPFR_RNDN);
        else mpfr_set_q(ex, as_rat(exp).value().get_mpq_t(), MPFR_RNDN);

        mpfr_pow(rx, bx, ex, MPFR_RNDN);
        auto result = make<Float>(static_cast<mpfr_srcptr>(rx), dps);
        mpfr_clear(bx);
        mpfr_clear(ex);
        mpfr_clear(rx);
        return result;
    }

    return std::nullopt;
}

}  // namespace sympp
