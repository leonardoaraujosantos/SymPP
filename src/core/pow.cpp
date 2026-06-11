#include <sympp/core/pow.hpp>

#include <cstddef>
#include <optional>
#include <string>
#include <utility>

#include <gmpxx.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/imaginary_unit.hpp>
#include <sympp/core/infinity.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/number_arith.hpp>
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
    switch (k) {
        case AssumptionKey::Real: {
            // positive base + real exp → real
            if (base->ask(AssumptionKey::Positive) == true
                && exp->ask(AssumptionKey::Real) == true) return true;
            // real base + integer exp → real (when 0^neg can't occur)
            if (base->ask(AssumptionKey::Real) == true
                && exp->ask(AssumptionKey::Integer) == true
                && (base->ask(AssumptionKey::Nonzero) == true
                    || exp->ask(AssumptionKey::Nonnegative) == true)) {
                return true;
            }
            return std::nullopt;
        }
        case AssumptionKey::Positive:
            // Any real-exp power of a positive is positive.
            if (base->ask(AssumptionKey::Positive) == true) return true;
            return std::nullopt;
        case AssumptionKey::Negative:
            return std::nullopt;
        case AssumptionKey::Zero:
            if (base->ask(AssumptionKey::Zero) == true
                && exp->ask(AssumptionKey::Positive) == true) return true;
            if (base->ask(AssumptionKey::Nonzero) == true) return false;
            return std::nullopt;
        case AssumptionKey::Nonzero:
            if (base->ask(AssumptionKey::Nonzero) == true) return true;
            return std::nullopt;
        case AssumptionKey::Integer:
            // integer base + nonneg integer exp → integer
            if (base->ask(AssumptionKey::Integer) == true
                && exp->ask(AssumptionKey::Integer) == true
                && exp->ask(AssumptionKey::Nonnegative) == true) return true;
            return std::nullopt;
        case AssumptionKey::Rational:
            // rational base + integer exp → rational
            if (base->ask(AssumptionKey::Rational) == true
                && exp->ask(AssumptionKey::Integer) == true) return true;
            return std::nullopt;
        case AssumptionKey::Finite:
            // finite base + finite nonneg exp → finite (skip 0^neg edge case
            // by requiring nonneg or nonzero base)
            if (base->ask(AssumptionKey::Finite) == true
                && exp->ask(AssumptionKey::Finite) == true
                && (base->ask(AssumptionKey::Nonzero) == true
                    || exp->ask(AssumptionKey::Nonnegative) == true)) {
                return true;
            }
            return std::nullopt;
        case AssumptionKey::Nonnegative:
        case AssumptionKey::Nonpositive:
            // Could derive in some special cases (e.g. positive base) but
            // those reduce to is_positive which is covered above.
            if (base->ask(AssumptionKey::Positive) == true) {
                if (k == AssumptionKey::Nonnegative) return true;
                if (k == AssumptionKey::Nonpositive) return false;
            }
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
    if (q.numerator() != 1) return std::nullopt;  // only 1/n for now
    auto den = q.denominator();
    if (!den.fits_ulong_p()) return std::nullopt;
    auto n = den.get_ui();
    if (n < 2) return std::nullopt;

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
        if (auto r = exact_root(z.value())) return make<Integer>(std::move(*r));
        return std::nullopt;
    }
    if (base->type_id() == TypeId::Rational) {
        const auto& r = static_cast<const Rational&>(*base);
        if (r.is_negative()) return std::nullopt;  // branch handling — defer
        // Denominator is always positive (Rational invariant), numerator ≥ 0
        // here, so both roots are over non-negative integers.
        auto rn = exact_root(r.numerator());
        auto rd = exact_root(r.denominator());
        if (rn && rd) return rational(*rn, *rd);
        return std::nullopt;
    }
    return std::nullopt;
}

// Square-root factor extraction for a positive Integer or Rational radicand
// that is NOT a perfect square: √N = s·√m with N = s²·m. For a rational p/q it
// uses √(p/q) = √(p·q)/q, yielding SymPy's rationalised form (√(2/3) = √6/3,
// √(8/9) = 2√2/3). Scoped to square roots; returns nullopt when nothing pulls
// out (s == 1) or the radicand is too large to factor cheaply.
[[nodiscard]] std::optional<Expr> try_sqrt_factor_extraction(const Expr& base,
                                                             const Expr& exp) {
    if (exp->type_id() != TypeId::Rational) return std::nullopt;
    const auto& q = static_cast<const Rational&>(*exp);
    if (!(q.numerator() == 1 && q.denominator() == 2)) return std::nullopt;

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

    // Radicand under a single root: √(num/den) = √(num·den)/den.
    mpz_class radicand = num * den;

    // Bound the trial division so a huge radicand can never hang the process.
    constexpr unsigned long kMaxRadicand = 1000000000000UL;  // 1e12 → ≤1e6 iters
    if (mpz_cmp_ui(radicand.get_mpz_t(), kMaxRadicand) > 0) return std::nullopt;

    // Pull out the largest square factor: radicand = s² · m, m square-free.
    mpz_class s(1), m(radicand);
    for (mpz_class d(2); d * d <= m; ++d) {
        mpz_class d2 = d * d;
        while (mpz_divisible_p(m.get_mpz_t(), d2.get_mpz_t())) {
            s *= d;
            m /= d2;
        }
    }
    // For an integer radicand (den == 1) with no square factor (s == 1) there
    // is nothing to do — leave it symbolic. A rational radicand still gets
    // rationalised even when s == 1 (√(2/3) = √6/3), so only the integer case
    // bails here.
    if (s == 1 && den == 1) return std::nullopt;

    // (s / den) · √m.  den == 1 for an integer radicand → plain s·√m. √m stays
    // symbolic: m is square-free, so this recursion bottoms out immediately.
    Expr radical = pow(make<Integer>(m), exp);
    return mul(rational(s, den), radical);
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

    // ---- numeric base & numeric exp ----
    if (is_number(base) && is_number(exp)) {
        auto r = number_pow(static_cast<const Number&>(*base),
                           static_cast<const Number&>(*exp));
        if (r) return *r;
        // Fallback: Integer / Rational base ^ Rational(1, n) — perfect n-th root.
        if (auto root = try_perfect_root(base, exp); root.has_value()) {
            return *root;
        }
        // √N where N is not a perfect square — pull out the square factor.
        if (auto ext = try_sqrt_factor_extraction(base, exp); ext.has_value()) {
            return *ext;
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
