#include <sympp/core/pow.hpp>

#include <cstddef>
#include <optional>
#include <string>
#include <utility>

#include <gmpxx.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
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
    return t == TypeId::Add || t == TypeId::Mul || t == TypeId::Pow;
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

// Detect Integer^Rational(1, n) collapsing to a perfect n-th root. Returns
// the simplified Expr, or std::nullopt if no exact root exists.
//
// Reference: sympy/core/power.py — Pow.__new__ / _eval_power for unit
// numerator rational exponents on integer bases.
[[nodiscard]] std::optional<Expr> try_integer_perfect_root(const Expr& base,
                                                           const Expr& exp) {
    if (base->type_id() != TypeId::Integer) return std::nullopt;
    if (exp->type_id() != TypeId::Rational) return std::nullopt;

    const auto& z = static_cast<const Integer&>(*base);
    const auto& q = static_cast<const Rational&>(*exp);
    if (z.is_negative()) return std::nullopt;  // branch handling — defer
    if (q.numerator() != 1) return std::nullopt;  // only 1/n for now
    auto den = q.denominator();
    if (!den.fits_ulong_p()) return std::nullopt;
    auto n = den.get_ui();
    if (n < 2) return std::nullopt;

    mpz_class root;
    mpz_class rem;
    mpz_rootrem(root.get_mpz_t(), rem.get_mpz_t(), z.value().get_mpz_t(), n);
    if (sgn(rem) != 0) return std::nullopt;
    return make<Integer>(std::move(root));
}

}  // namespace

Expr pow(const Expr& base, const Expr& exp) {
    // ---- universal identity rules ----
    // x^1 → x
    if (exp == S::One()) return base;
    // x^0 → 1  (including 0^0 by SymPy convention)
    if (exp == S::Zero()) return S::One();
    // 1^x → 1
    if (base == S::One()) return S::One();

    // ---- numeric base & numeric exp ----
    if (is_number(base) && is_number(exp)) {
        auto r = number_pow(static_cast<const Number&>(*base),
                           static_cast<const Number&>(*exp));
        if (r) return *r;
        // Fallback: Integer base ^ Rational(1, n) — perfect n-th root.
        if (auto root = try_integer_perfect_root(base, exp); root.has_value()) {
            return *root;
        }
    }

    // 0^positive → 0  (handled by number_pow above for numeric exp; only need
    // the symbolic-exp fallback when exp is e.g. a Symbol with positive
    // assumption — defer to Phase 2)

    // Otherwise stay unevaluated.
    return make<Pow>(base, exp);
}

}  // namespace sympp
