#include <sympp/core/pow.hpp>

#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <gmpxx.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/imaginary_unit.hpp>
#include <sympp/core/infinity.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

constexpr std::size_t kPowHashSeed = 0xF00D'FACE'7777'8888ULL;

[[nodiscard]] std::size_t hash_combine(std::size_t a, std::size_t b) noexcept {
    return a ^ (b + 0x9E37'79B9'7F4A'7C15ULL + (a << 6) + (a >> 2));
}

[[nodiscard]] bool needs_parens_in_pow_base(const Expr& e) noexcept {
    if (!e) return false;
    auto t = e->type_id();
    if (t == TypeId::Add || t == TypeId::Mul || t == TypeId::Pow) return true;
    // Negative numeric base needs parens or `-12**(1/2)` parses as
    // `-(12**(1/2))` due to operator precedence. (-12)**(1/2) is correct.
    if (t == TypeId::Integer) {
        return static_cast<const Integer&>(*e).is_negative();
    }
    // Rational p/q always needs parens — `p/q**e` parses as `p / q**e`.
    if (t == TypeId::Rational) return true;
    if (t == TypeId::Float) {
        return static_cast<const Float&>(*e).is_negative();
    }
    return false;
}

[[nodiscard]] bool needs_parens_in_pow_exp(const Expr& e) noexcept {
    if (!e) return false;
    auto t = e->type_id();
    if (t == TypeId::Add || t == TypeId::Mul) return true;
    // Rationals (1/2) print with a slash that would associate wrong without
    // parens. Negative integer exponents likewise.
    if (t == TypeId::Rational) return true;
    if (t == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*e);
        return z.is_negative();
    }
    if (t == TypeId::Float) {
        const auto& f = static_cast<const Float&>(*e);
        return f.is_negative();
    }
    return false;
}

}  // namespace

Pow::Pow(Expr base, Expr exp) : args_{std::move(base), std::move(exp)} {
    hash_ = hash_combine(kPowHashSeed,
                        hash_combine(args_[0]->hash(), args_[1]->hash()));
}

bool Pow::equals(const Basic& other) const noexcept {
    if (other.type_id() != TypeId::Pow) return false;
    const auto& o = static_cast<const Pow&>(other);
    return args_[0] == o.args_[0] && args_[1] == o.args_[1];
}

std::optional<bool> Pow::ask(AssumptionKey k) const noexcept {
    const auto& base = args_[0];
    const auto& exp = args_[1];
    // A positive even integer exponent: x^(2m) ≥ 0 for any real x (and > 0 when
    // x ≠ 0). Used by the sign queries below.
    const bool even_pos_int_exp =
        exp->type_id() == TypeId::Integer
        && static_cast<const Integer&>(*exp).is_positive()
        && static_cast<const Integer&>(*exp).fits_long()
        && (static_cast<const Integer&>(*exp).to_long() % 2 == 0);
    // A Pow that survives construction with a positive rational base ≠ 1 and a
    // non-integer rational exponent is a real algebraic irrational: the factory
    // has already reduced every perfect rational root (√4 → 2, √(9/4) → 3/2) and
    // rationalized rational bases, so what remains (√2, 2^(1/3), 6^(1/2)) is not
    // rational. (A Rational-typed exponent is necessarily non-integer — an
    // integer-valued ratio normalizes to the Integer type.)
    const bool surviving_algebraic_root =
        (base->type_id() == TypeId::Integer
         || base->type_id() == TypeId::Rational)
        && !(base == S::One())
        && direct_ask(base, AssumptionKey::Positive) == true
        && exp->type_id() == TypeId::Rational;
    switch (k) {
        case AssumptionKey::Complex:
            if (direct_ask(base, AssumptionKey::Complex) == true
                && direct_ask(exp, AssumptionKey::Integer) == true
                && direct_ask(base, AssumptionKey::Nonzero) == true) return true;
            return std::nullopt;
        case AssumptionKey::Imaginary: {
            // (imaginary)^n: odd positive integer n → imaginary (i³ = −i),
            // even → real. Requires a nonzero base.
            const bool pos_int =
                exp->type_id() == TypeId::Integer
                && static_cast<const Integer&>(*exp).is_positive()
                && static_cast<const Integer&>(*exp).fits_long();
            if (pos_int && direct_ask(base, AssumptionKey::Imaginary) == true) {
                long n = static_cast<const Integer&>(*exp).to_long();
                return (n % 2 == 1);
            }
            return std::nullopt;
        }

        case AssumptionKey::Real: {
            // positive base + real exp → real
            if (direct_ask(base, AssumptionKey::Positive) == true
                && direct_ask(exp, AssumptionKey::Real) == true) return true;
            // (imaginary)^(even positive integer) → real (i² = −1).
            if (exp->type_id() == TypeId::Integer
                && static_cast<const Integer&>(*exp).is_positive()
                && static_cast<const Integer&>(*exp).fits_long()
                && (static_cast<const Integer&>(*exp).to_long() % 2 == 0)
                && direct_ask(base, AssumptionKey::Imaginary) == true) {
                return true;
            }
            // real base + integer exp → real (when 0^neg can't occur)
            if (direct_ask(base, AssumptionKey::Real) == true
                && direct_ask(exp, AssumptionKey::Integer) == true
                && (direct_ask(base, AssumptionKey::Nonzero) == true
                    || direct_ask(exp, AssumptionKey::Nonnegative) == true)) {
                return true;
            }
            return std::nullopt;
        }
        case AssumptionKey::Positive:
            // Any real-exp power of a positive is positive.
            if (direct_ask(base, AssumptionKey::Positive) == true) return true;
            // Even power of a nonzero real base: x^(2m) > 0.
            if (even_pos_int_exp && direct_ask(base, AssumptionKey::Real) == true
                && direct_ask(base, AssumptionKey::Nonzero) == true) {
                return true;
            }
            return std::nullopt;
        case AssumptionKey::Negative:
            // A negative real base raised to an odd integer power stays negative
            // (an even power turns positive — handled by the Positive case).
            if (direct_ask(base, AssumptionKey::Negative) == true
                && exp->type_id() == TypeId::Integer) {
                const auto& z = static_cast<const Integer&>(*exp);
                if (z.fits_long()) return z.to_long() % 2 != 0;
            }
            return std::nullopt;
        case AssumptionKey::Zero:
            if (direct_ask(base, AssumptionKey::Zero) == true
                && direct_ask(exp, AssumptionKey::Positive) == true) return true;
            if (direct_ask(base, AssumptionKey::Nonzero) == true) return false;
            return std::nullopt;
        case AssumptionKey::Nonzero:
            if (direct_ask(base, AssumptionKey::Nonzero) == true) return true;
            return std::nullopt;
        case AssumptionKey::Integer:
            // integer base + nonneg integer exp → integer
            if (direct_ask(base, AssumptionKey::Integer) == true
                && direct_ask(exp, AssumptionKey::Integer) == true
                && direct_ask(exp, AssumptionKey::Nonnegative) == true) return true;
            return std::nullopt;
        case AssumptionKey::Rational:
            // rational base + integer exp → rational
            if (direct_ask(base, AssumptionKey::Rational) == true
                && direct_ask(exp, AssumptionKey::Integer) == true) return true;
            // A surviving root: the factory already reduces every perfect rational
            // root (√4 → 2, √(9/4) → 3/2) and rationalizes rational bases, so a Pow
            // that persists with a positive rational base ≠ 1 and a non-integer
            // rational exponent (√2, 2^(1/3), 6^(1/2)) is a genuine algebraic
            // irrational — hence not rational.
            if (surviving_algebraic_root) return false;
            return std::nullopt;
        case AssumptionKey::Finite:
            // finite base + finite nonneg exp → finite (skip 0^neg edge case
            // by requiring nonneg or nonzero base)
            if (direct_ask(base, AssumptionKey::Finite) == true
                && direct_ask(exp, AssumptionKey::Finite) == true
                && (direct_ask(base, AssumptionKey::Nonzero) == true
                    || direct_ask(exp, AssumptionKey::Nonnegative) == true)) {
                return true;
            }
            return std::nullopt;
        case AssumptionKey::Nonnegative:
        case AssumptionKey::Nonpositive:
            // Could derive in some special cases (e.g. positive base) but
            // those reduce to is_positive which is covered above.
            if (direct_ask(base, AssumptionKey::Positive) == true) {
                if (k == AssumptionKey::Nonnegative) return true;
                if (k == AssumptionKey::Nonpositive) return false;
            }
            // Even power of a real base: x^(2m) ≥ 0 (nonpositive only if x = 0).
            if (even_pos_int_exp && direct_ask(base, AssumptionKey::Real) == true) {
                if (k == AssumptionKey::Nonnegative) return true;
                if (k == AssumptionKey::Nonpositive
                    && direct_ask(base, AssumptionKey::Nonzero) == true) {
                    return false;
                }
            }
            return std::nullopt;
        case AssumptionKey::Even:
        case AssumptionKey::Odd: {
            // base^n for a positive integer n preserves the base's parity
            // (oddⁿ = odd, evenⁿ = even). Non-positive / non-integer exponents
            // are left Unknown (e.g. x^(-1) isn't generally an integer).
            const bool pos_int_exp =
                exp->type_id() == TypeId::Integer
                && static_cast<const Integer&>(*exp).is_positive();
            if (!pos_int_exp) return std::nullopt;
            if (base->ask(k) == true) return true;
            if (k == AssumptionKey::Even
                && direct_ask(base, AssumptionKey::Odd) == true) {
                return false;
            }
            if (k == AssumptionKey::Odd
                && direct_ask(base, AssumptionKey::Even) == true) {
                return false;
            }
            return std::nullopt;
        }
        case AssumptionKey::Prime:
        case AssumptionKey::Composite:
            // Primality / compositeness of a symbolic power isn't decided
            // structurally.
            return std::nullopt;
        case AssumptionKey::Algebraic: {
            // algebraic base ^ rational exponent ⇒ algebraic (radicals, rational
            // powers). Guard 0^negative: require a nonnegative exponent or a
            // nonzero base.
            if (direct_ask(base, AssumptionKey::Algebraic) == true
                && direct_ask(exp, AssumptionKey::Rational) == true
                && (direct_ask(exp, AssumptionKey::Nonnegative) == true
                    || direct_ask(base, AssumptionKey::Nonzero) == true)) {
                return true;
            }
            // transcendental base ^ nonzero rational exponent ⇒ transcendental.
            if (direct_ask(base, AssumptionKey::Transcendental) == true
                && direct_ask(exp, AssumptionKey::Rational) == true
                && direct_ask(exp, AssumptionKey::Nonzero) == true) {
                return false;
            }
            return std::nullopt;
        }
        case AssumptionKey::Transcendental: {
            if (direct_ask(base, AssumptionKey::Transcendental) == true
                && direct_ask(exp, AssumptionKey::Rational) == true
                && direct_ask(exp, AssumptionKey::Nonzero) == true) {
                return true;
            }
            return std::nullopt;
        }
        case AssumptionKey::Irrational:
            // A surviving root of a positive rational base (√2, 2^(1/3)) is a real
            // algebraic irrational; see the Rational case above.
            if (surviving_algebraic_root) return true;
            return std::nullopt;
        case AssumptionKey::ExtendedReal:
        case AssumptionKey::Infinite:
            // Left to the generic derivation layer.
            return std::nullopt;
    }
    return std::nullopt;
}

std::string Pow::str() const {
    std::string b = args_[0]->str();
    std::string e = args_[1]->str();
    if (needs_parens_in_pow_base(args_[0])) b = "(" + b + ")";
    if (needs_parens_in_pow_exp(args_[1])) e = "(" + e + ")";
    return b + "**" + e;
}

namespace {

// Detect base^Rational(1, n) collapsing to a perfect n-th root, for a
// non-negative Integer or Rational base. For a rational p/q the numerator and
// denominator are rooted independently (so √(9/4) → 3/2). Returns the
// simplified Expr, or std::nullopt if no exact root exists.
//
// Reference: sympy/core/power.py — Pow.__new__ / _eval_power for unit
// numerator rational exponents on integer / rational bases.
[[nodiscard]] std::optional<Expr> try_perfect_root(const Expr& base,
                                                   const Expr& exp) {
    if (exp->type_id() != TypeId::Rational) return std::nullopt;
    const auto& q = static_cast<const Rational&>(*exp);
    auto den = q.denominator();
    if (!den.fits_ulong_p()) return std::nullopt;
    auto n = den.get_ui();
    if (n < 2) return std::nullopt;
    // Numerator p: base^(p/n) = (base^(1/n))^p when the n-th root is exact.
    const mpz_class& p = q.numerator();
    if (!p.fits_slong_p()) return std::nullopt;
    const long pl = p.get_si();

    // Exact n-th root of a non-negative integer, or nullopt.
    auto exact_root = [n](const mpz_class& z) -> std::optional<mpz_class> {
        mpz_class root, rem;
        mpz_rootrem(root.get_mpz_t(), rem.get_mpz_t(), z.get_mpz_t(), n);
        if (sgn(rem) != 0) return std::nullopt;
        return root;
    };

    if (base->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*base);
        if (z.is_negative()) return std::nullopt;  // branch handling — defer
        if (auto r = exact_root(z.value())) {
            // (root)^p — pow folds integer^integer (negative p → Rational).
            return pow(make<Integer>(std::move(*r)), make<Integer>(pl));
        }
        return std::nullopt;
    }
    if (base->type_id() == TypeId::Rational) {
        const auto& r = static_cast<const Rational&>(*base);
        if (r.is_negative()) return std::nullopt;  // branch handling — defer
        // Denominator is always positive (Rational invariant), numerator ≥ 0
        // here, so both roots are over non-negative integers.
        auto rn = exact_root(r.numerator());
        auto rd = exact_root(r.denominator());
        if (rn && rd) return pow(rational(*rn, *rd), make<Integer>(pl));
        return std::nullopt;
    }
    return std::nullopt;
}

// N-th-root factor extraction for a positive Integer or Rational radicand under
// a unit 1/n power (n ≥ 2): N^(1/n) = s · m^(1/n) with N = sⁿ · m. For a rational
// p/q it rationalises via (p/q)^(1/n) = (p·q^(n-1))^(1/n)/q — the n-th-root
// generalisation of √(p/q) = √(p·q)/q (so √(2/3) = √6/3, cbrt(16) = 2·cbrt(2),
// (2/3)^(1/3) = 18^(1/3)/3). Returns nullopt when nothing pulls out (s == 1,
// integer radicand) or the radicand is too large to factor cheaply.
[[nodiscard]] std::optional<Expr> try_nth_root_factor_extraction(const Expr& base,
                                                                 const Expr& exp) {
    if (exp->type_id() != TypeId::Rational) return std::nullopt;
    const auto& q = static_cast<const Rational&>(*exp);
    if (!(q.numerator() == 1)) return std::nullopt;
    const mpz_class& denq = q.denominator();
    if (!denq.fits_ulong_p()) return std::nullopt;
    const unsigned long n = denq.get_ui();
    if (n < 2) return std::nullopt;

    mpz_class num, den(1);
    if (base->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*base);
        if (sgn(z.value()) <= 0) return std::nullopt;
        num = z.value();
    } else if (base->type_id() == TypeId::Rational) {
        const auto& r = static_cast<const Rational&>(*base);
        if (r.is_negative()) return std::nullopt;
        num = r.numerator();
        den = r.denominator();
    } else {
        return std::nullopt;
    }

    // Radicand under a single n-th root: (num/den)^(1/n) = (num·den^(n-1))^(1/n)/den.
    mpz_class radicand = num;
    if (den != 1) {
        mpz_class denpow;
        mpz_pow_ui(denpow.get_mpz_t(), den.get_mpz_t(), n - 1);
        radicand *= denpow;
    }

    // Bound the trial division so a huge radicand can never hang the process.
    constexpr unsigned long kMaxRadicand = 1000000000000UL;  // 1e12 → ≤1e6 iters
    if (mpz_cmp_ui(radicand.get_mpz_t(), kMaxRadicand) > 0) return std::nullopt;

    // Pull out the largest n-th-power factor: radicand = sⁿ · m, m n-th-power-free.
    // d ranges over candidate factors; dⁿ grows fast, so the loop is bounded by
    // m^(1/n) (≤ 1e6 iterations for n=2, fewer for larger n).
    mpz_class s(1), m(radicand);
    for (mpz_class d(2);; ++d) {
        mpz_class dn;
        mpz_pow_ui(dn.get_mpz_t(), d.get_mpz_t(), n);
        if (dn > m) break;
        while (mpz_divisible_p(m.get_mpz_t(), dn.get_mpz_t())) {
            s *= d;
            m /= dn;
        }
    }
    // For an integer radicand (den == 1) with no factor (s == 1) there is nothing
    // to do — leave it symbolic. A rational radicand still gets rationalised even
    // when s == 1 (√(2/3) = √6/3), so only the integer case bails here.
    if (s == 1 && den == 1) return std::nullopt;

    // (s / den) · m^(1/n).  m is n-th-power-free, so the radical stays symbolic
    // and this bottoms out immediately.
    Expr radical = pow(make<Integer>(m), exp);
    return mul(rational(s, den), radical);
}

// Non-unit rational power of a positive integer base: bᵖ⸍۹ → coeff · ∏ pᵢ^(rᵢ/q),
// the n-th-root extraction generalised to numerator p ≥ 2. Prime-factorise
// b = ∏ pᵢ^aᵢ; each prime's exponent under the power is aᵢ·p/q, whose integer
// part floor(aᵢp/q) is pulled into the integer coefficient and whose remainder
// rᵢ = aᵢp mod q stays under a per-prime radical pᵢ^(rᵢ/q). Keeping primes
// separate matches SymPy's canonical form (16^(2/3) = 4·2^(2/3), not 4·4^(1/3);
// 2^(5/2) = 4√2). Returns nullopt when nothing pulls out (every aᵢp < q, e.g.
// 2^(2/3)) — leaving such already-reduced powers symbolic. p = 1 is NROOT-1 and
// the perfect case is try_perfect_root, both dispatched earlier.
// Pull a perfect-power factor out of a NEGATIVE integer base under a rational
// power, keeping the principal-branch sign on the residual: with a = |base| =
// sᵠ·m (m q-th-power-free) and exponent p/q, (−a)^(p/q) = sᵖ·(−m)^(p/q). So
// (−8)^(1/3) = 2·(−1)^(1/3), (−24)^(1/3) = 2·(−3)^(1/3), (−8)^(2/3) = 4·(−1)^(2/3).
// A base with no perfect-power factor (s = 1, e.g. (−2)^(1/3)) is left symbolic.
// The residual (−m)^(p/q) re-enters here with s = 1, so there is no recursion.
[[nodiscard]] std::optional<Expr> try_negative_root_factor_extraction(
    const Expr& base, const Expr& exp) {
    if (base->type_id() != TypeId::Integer) return std::nullopt;
    const auto& z = static_cast<const Integer&>(*base);
    if (!z.is_negative()) return std::nullopt;
    if (exp->type_id() != TypeId::Rational) return std::nullopt;
    const auto& q = static_cast<const Rational&>(*exp);
    const mpz_class& denq = q.denominator();
    if (!denq.fits_ulong_p()) return std::nullopt;
    const unsigned long n = denq.get_ui();  // n == q (the root order)
    if (n < 2) return std::nullopt;

    mpz_class a = -z.value();  // |base| > 0
    constexpr unsigned long kMaxRadicand = 1000000000000UL;  // bound trial division
    if (mpz_cmp_ui(a.get_mpz_t(), kMaxRadicand) > 0) return std::nullopt;

    // a = sⁿ · m, m n-th-power-free.
    mpz_class s(1), m(a);
    for (mpz_class d(2);; ++d) {
        mpz_class dn;
        mpz_pow_ui(dn.get_mpz_t(), d.get_mpz_t(), n);
        if (dn > m) break;
        while (mpz_divisible_p(m.get_mpz_t(), dn.get_mpz_t())) {
            s *= d;
            m /= dn;
        }
    }
    if (s == 1) return std::nullopt;  // nothing pulls out — keep (−a)^(p/q)

    // sᵖ · (−m)^(p/q).  −m is n-th-power-free, so the residual radical bottoms out.
    Expr coeff = pow(make<Integer>(std::move(s)), make<Integer>(q.numerator()));
    Expr radical = pow(make<Integer>(-m), exp);
    return mul(coeff, radical);
}

[[nodiscard]] std::optional<Expr> try_rational_power_extraction(const Expr& base,
                                                                const Expr& exp) {
    if (exp->type_id() != TypeId::Rational) return std::nullopt;
    const auto& q = static_cast<const Rational&>(*exp);
    const mpz_class& pnum = q.numerator();
    const mpz_class& qden = q.denominator();
    if (pnum < 2) return std::nullopt;                 // p = 1 → NROOT-1
    if (!qden.fits_ulong_p() || !pnum.fits_ulong_p()) return std::nullopt;
    const unsigned long qd = qden.get_ui();
    if (qd < 2) return std::nullopt;
    // Bound the numerator so the pulled-out coefficient can't blow up.
    if (mpz_cmp_ui(pnum.get_mpz_t(), 1000UL) > 0) return std::nullopt;

    // Integer base b → factorise b directly. Rational base a/b → use the
    // identity (a/b)^(P/q) = (a·b^(q-1))^(P/q) / b^P, reducing to the integer
    // case on a·b^(q-1) with a final division by den_pow = b^P (which
    // rationalises the denominator, e.g. (2/3)^(2/3) = 2^(2/3)·3^(1/3)/3).
    mpz_class m;        // integer base to factorise
    mpz_class den_pow(1);
    if (base->type_id() == TypeId::Integer) {
        const auto& zb = static_cast<const Integer&>(*base);
        if (sgn(zb.value()) <= 0) return std::nullopt;  // negative/zero → defer
        m = zb.value();
    } else if (base->type_id() == TypeId::Rational) {
        const auto& rb = static_cast<const Rational&>(*base);
        if (rb.is_negative()) return std::nullopt;
        mpz_class bq1;
        mpz_pow_ui(bq1.get_mpz_t(), rb.denominator().get_mpz_t(), qd - 1);
        m = rb.numerator() * bq1;
        mpz_pow_ui(den_pow.get_mpz_t(), rb.denominator().get_mpz_t(), pnum.get_ui());
    } else {
        return std::nullopt;
    }
    if (m == 1) return std::nullopt;
    constexpr unsigned long kMaxBase = 1000000000000UL;  // 1e12 → ≤1e6 trial iters
    if (mpz_cmp_ui(m.get_mpz_t(), kMaxBase) > 0) return std::nullopt;

    mpz_class coeff(1);
    std::vector<std::pair<mpz_class, unsigned long>> residuals;  // (prime, rᵢ)
    bool pulled = false;

    // Accumulate the integer coefficient and the per-prime residual exponents.
    // Crucially this builds NO Expr yet — the residual pow(pᵢ, rᵢ/q) factors
    // re-enter this function (e.g. building 2^(2/3) while reducing 16^(2/3)),
    // and a bare prime power must bail at the `!pulled` check below BEFORE any
    // Pow is constructed, or it would recurse without bound.
    auto handle_prime = [&](const mpz_class& prime, unsigned long a) {
        mpz_class ap = mpz_class(a) * pnum;
        mpz_class ipart = ap / qd;                 // floor (both positive)
        unsigned long r = mpz_class(ap % qd).get_ui();
        if (ipart > 0) {
            mpz_class pw;
            mpz_pow_ui(pw.get_mpz_t(), prime.get_mpz_t(), ipart.get_ui());
            coeff *= pw;
            pulled = true;
        }
        if (r > 0) residuals.emplace_back(prime, r);
    };

    // Trial-division prime factorisation of the base.
    for (mpz_class d(2); d * d <= m; ++d) {
        if (!mpz_divisible_p(m.get_mpz_t(), d.get_mpz_t())) continue;
        unsigned long a = 0;
        while (mpz_divisible_p(m.get_mpz_t(), d.get_mpz_t())) {
            m /= d;
            ++a;
        }
        handle_prime(d, a);
    }
    if (m > 1) handle_prime(m, 1);  // leftover prime factor > √base

    // For an integer base with nothing pulled out (e.g. 2^(2/3)) leave it
    // symbolic. A rational base always proceeds — the den_pow division
    // rationalises it even when no prime contributes an integer part.
    if (!pulled && den_pow == 1) return std::nullopt;
    Expr result = (den_pow == 1) ? make<Integer>(coeff)
                                 : rational(coeff, den_pow);
    for (const auto& [prime, r] : residuals) {
        result = mul(result,
                     pow(make<Integer>(prime), rational(mpz_class(r), mpz_class(qd))));
    }
    return result;
}

// Principal square root of a negative Integer/Rational: √(−a) = I·√a (a > 0).
// Restricted to the ½ power; other fractional powers of a negative base need
// full branch-cut handling and are left symbolic.
[[nodiscard]] std::optional<Expr> try_sqrt_of_negative(const Expr& base,
                                                       const Expr& exp) {
    if (!(exp == rational(1, 2))) return std::nullopt;
    bool negative = false;
    if (base->type_id() == TypeId::Integer) {
        negative = static_cast<const Integer&>(*base).is_negative();
    } else if (base->type_id() == TypeId::Rational) {
        negative = static_cast<const Rational&>(*base).is_negative();
    } else {
        return std::nullopt;
    }
    if (!negative) return std::nullopt;
    // √(−a) = I·√a; pow(|base|, ½) reuses the perfect-root / factor-extraction
    // paths, so the magnitude comes back fully reduced (√(−8) → 2·√2·I).
    Expr magnitude = mul(S::NegativeOne(), base);
    return mul(S::I(), pow(magnitude, exp));
}

// Rationalise an inverse square root: base^(-1/2) → sqrt of the reciprocal,
// which the positive-½ factor-extraction path then reduces. Examples:
// 3^(-1/2) = sqrt(1/3) = sqrt(3)/3, (2/3)^(-1/2) = sqrt(3/2) = sqrt(6)/2,
// 12^(-1/2) = sqrt(1/12) = sqrt(3)/6. Positive Integer/Rational bases only;
// a negative base needs branch-cut handling and is left symbolic.
[[nodiscard]] std::optional<Expr> try_inverse_sqrt(const Expr& base,
                                                   const Expr& exp) {
    if (!(exp == rational(-1, 2))) return std::nullopt;
    if (base->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*base);
        if (sgn(z.value()) <= 0) return std::nullopt;
        return pow(rational(mpz_class(1), z.value()), rational(1, 2));
    }
    if (base->type_id() == TypeId::Rational) {
        const auto& r = static_cast<const Rational&>(*base);
        if (r.is_negative() || r.numerator() == 0) return std::nullopt;
        // Reciprocal den/num, then the positive-½ path rationalises it.
        return pow(rational(r.denominator(), r.numerator()), rational(1, 2));
    }
    return std::nullopt;
}

}  // namespace

Expr pow(const Expr& base, const Expr& exp) {
    // ---- universal identity rules ----
    // x^1 → x
    if (exp == S::One()) return base;
    // x^0 → 1  (including 0^0 by SymPy convention)
    if (exp == S::Zero()) return S::One();
    // Infinity / nan in the base or exponent (1^oo = nan, so this must run
    // before the 1^x → 1 rule below).
    if (auto inf = pow_with_infinity(base, exp); inf.has_value()) return *inf;
    // 0^(negative) → zoo (the value of 1/0). 0^0 = 1 (handled above) and
    // 0^(positive) = 0 (number_pow) are fine; only division-by-zero needs
    // folding so it never escapes as the malformed 0**(-1). (Infinity exponents
    // are already resolved by pow_with_infinity above.)
    if (base == S::Zero() && is_negative(exp) == true) {
        return S::ComplexInfinity();
    }
    // 1^x → 1
    if (base == S::One()) return S::One();

    // ---- (a^m)^n → a^(m*n) when n is an Integer ----
    // Safe for integer outer exponent (avoids the (x²)^(1/2) ≠ x issue
    // for negative x). Common case: pow(pow(y, 2), -1) → pow(y, -2).
    if (exp->type_id() == TypeId::Integer
        && base->type_id() == TypeId::Pow) {
        const auto& inner_exp = base->args()[1];
        Expr new_exp = mul(inner_exp, exp);
        return pow(base->args()[0], new_exp);
    }

    // ---- Distribute an integer power over a symbol-free (numeric) Mul base ----
    // (∏ aᵢ)ⁿ = ∏ aᵢⁿ folds radical coefficients that the rationalised form leaves
    // split: (⅓·√3)⁻¹ = 3·3^(−½) = √3, so reciprocals of such products reduce (this
    // unblocks e.g. acot(√3/3) = atan(√3) = π/3). Restricted to a base with no free
    // symbols, so the compact form of (2·x)ⁿ is untouched.
    if (exp->type_id() == TypeId::Integer && base->type_id() == TypeId::Mul
        && free_symbols(base).empty()) {
        std::vector<Expr> factors;
        factors.reserve(base->args().size());
        for (const auto& f : base->args()) factors.push_back(pow(f, exp));
        return mul(std::move(factors));
    }

    // ---- I^Integer cycles through {1, I, -1, -I} ----
    if (base == S::I() && exp->type_id() == TypeId::Integer) {
        const auto& z = static_cast<const Integer&>(*exp);
        // Compute n mod 4 in [0, 3] without overflow concerns.
        mpz_class four(4);
        mpz_class r;
        mpz_mod(r.get_mpz_t(), z.value().get_mpz_t(), four.get_mpz_t());
        const long m = r.get_si();  // 0..3
        switch (m) {
            case 0: return S::One();
            case 1: return S::I();
            case 2: return S::NegativeOne();
            case 3: return mul(S::NegativeOne(), S::I());
        }
    }

    // ---- (−1)^k = 1 (k even) / −1 (k odd) for a provably-parity exponent ----
    // Closes (−1)^(2n) = 1, (−1)^(2n+1) = −1 for a known integer n.
    if (base == S::NegativeOne()) {
        if (is_provably_even(exp)) return S::One();
        if (is_provably_odd(exp)) return S::NegativeOne();
    }

    // ---- numeric base & numeric exp ----
    if (is_number(base) && is_number(exp)) {
        auto r = number_pow(static_cast<const Number&>(*base),
                           static_cast<const Number&>(*exp));
        if (r) return *r;
        // Fallback: Integer / Rational base ^ Rational(1, n) — perfect n-th root.
        if (auto root = try_perfect_root(base, exp); root.has_value()) {
            return *root;
        }
        // N^(1/n) where N is not a perfect n-th power — pull out the largest
        // n-th-power factor (cbrt(16) = 2·cbrt(2), √12 = 2√3).
        if (auto ext = try_nth_root_factor_extraction(base, exp); ext.has_value()) {
            return *ext;
        }
        // bᵖ⸍۹ with p ≥ 2 and a non-perfect base — pull out integer powers of each
        // prime factor (16^(2/3) = 4·2^(2/3), 2^(5/2) = 4√2).
        if (auto ext = try_rational_power_extraction(base, exp); ext.has_value()) {
            return *ext;
        }
        // (−a)^(p/q): pull a perfect q-th-power factor of a out, keeping the sign on
        // the residual — (−8)^(1/3) = 2·(−1)^(1/3). (Runs before the ½-only
        // try_sqrt_of_negative, which it subsumes for perfect-square radicands.)
        if (auto neg = try_negative_root_factor_extraction(base, exp);
            neg.has_value()) {
            return *neg;
        }
        // √(−a) = I·√a for a negative numeric base.
        if (auto neg = try_sqrt_of_negative(base, exp); neg.has_value()) {
            return *neg;
        }
        // base^(-1/2) → rationalised sqrt of the reciprocal (1/√3 = √3/3).
        if (auto inv = try_inverse_sqrt(base, exp); inv.has_value()) {
            return *inv;
        }
    }

    // 0^positive → 0  (handled by number_pow above for numeric exp; only need
    // the symbolic-exp fallback when exp is e.g. a Symbol with positive
    // assumption — defer to Phase 2)

    // Otherwise stay unevaluated.
    return make<Pow>(base, exp);
}

}  // namespace sympp
