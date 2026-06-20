#include <sympp/core/mul.hpp>

#include <algorithm>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <sympp/core/add.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/imaginary_unit.hpp>
#include <sympp/core/infinity.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>

#include "assumption_helpers.hpp"
#include "canonical_order.hpp"
#include "expr_map.hpp"

namespace sympp {

namespace {

// True if two of the factors are consecutive integers — they share the same
// non-constant part and their integer offsets differ by 1 (n and n+1, or 2n and
// 2n+1). One of two consecutive integers is even, so their product is even. The
// caller establishes that every factor is an integer first.
[[nodiscard]] bool has_consecutive_int_factors(
        const std::vector<Expr>& args) noexcept {
    // Decompose each factor into (rest, offset) with factor = rest + offset.
    std::vector<std::pair<Expr, long>> parts;
    parts.reserve(args.size());
    for (const auto& f : args) {
        if (f->type_id() == TypeId::Add) {
            std::vector<Expr> rest;
            long off = 0;
            bool got = false;
            bool ok = true;
            for (const auto& t : f->args()) {
                if (!got && t->type_id() == TypeId::Integer
                    && static_cast<const Integer&>(*t).fits_long()) {
                    off = static_cast<const Integer&>(*t).to_long();
                    got = true;
                } else {
                    rest.push_back(t);
                }
            }
            if (!ok || rest.empty()) continue;
            Expr r = rest.size() == 1 ? rest[0] : add(rest);
            parts.emplace_back(std::move(r), off);
        } else {
            parts.emplace_back(f, 0L);
        }
    }
    for (std::size_t i = 0; i < parts.size(); ++i) {
        for (std::size_t j = i + 1; j < parts.size(); ++j) {
            if (parts[i].first == parts[j].first) {
                const long d = parts[i].second - parts[j].second;
                if (d == 1 || d == -1) return true;
            }
        }
    }
    return false;
}

void flatten_into(const Expr& e, std::vector<Expr>& out) {
    if (e && e->type_id() == TypeId::Mul) {
        for (const auto& a : e->args()) {
            flatten_into(a, out);
        }
    } else {
        out.push_back(e);
    }
}

[[nodiscard]] std::size_t hash_combine(std::size_t a, std::size_t b) noexcept {
    return a ^ (b + 0x9E37'79B9'7F4A'7C15ULL + (a << 6) + (a >> 2));
}

constexpr std::size_t kMulHashSeed = 0xCAFE'BABE'5678'ABCDULL;

// Decompose `e` into (base, exponent) for Mul base collection.
//   x       -> (x, 1)
//   x**2    -> (x, 2)
//   x**y    -> (x, y)
//   y**(-1) -> (y, -1)
struct BaseExp {
    Expr base;
    Expr exp;
};

[[nodiscard]] BaseExp as_base_exp(const Expr& e) {
    if (e->type_id() == TypeId::Pow) {
        const auto& p = static_cast<const Pow&>(*e);
        return {p.base(), p.exp()};
    }
    return {e, S::One()};
}

}  // namespace

Mul::Mul(std::vector<Expr> args) : args_(std::move(args)) {
    compute_hash();
}

void Mul::compute_hash() noexcept {
    std::size_t h = kMulHashSeed;
    for (const auto& a : args_) {
        h = hash_combine(h, a->hash());
    }
    hash_ = h;
}

bool Mul::equals(const Basic& other) const noexcept {
    if (other.type_id() != TypeId::Mul) return false;
    const auto& o = static_cast<const Mul&>(other);
    if (args_.size() != o.args_.size()) return false;
    for (std::size_t i = 0; i < args_.size(); ++i) {
        if (!(args_[i] == o.args_[i])) return false;
    }
    return true;
}

std::optional<bool> Mul::ask(AssumptionKey k) const noexcept {
    using detail::all_args_have;
    using detail::any_arg_has;
    switch (k) {
        case AssumptionKey::ExtendedPositive:
        case AssumptionKey::ExtendedNegative:
        case AssumptionKey::ExtendedNonnegative:
        case AssumptionKey::ExtendedNonpositive:
        case AssumptionKey::Hermitian:
        case AssumptionKey::Antihermitian:
        case AssumptionKey::Commutative:
            return std::nullopt;  // derived by the generic ask() layer
        case AssumptionKey::Complex:
            // A product of finite complex factors is complex.
            if (all_args_have(args_, AssumptionKey::Complex, true)) return true;
            return std::nullopt;
        case AssumptionKey::Imaginary: {
            // Classify each factor as real or imaginary; the product is
            // imaginary iff an ODD number of factors are imaginary and the rest
            // are real, with every factor nonzero (I·real = imaginary,
            // I·I = real). A zero factor makes the product 0 (real, ¬imaginary).
            std::size_t imag = 0;
            bool all_classified = true;
            bool all_nonzero = true;
            for (const auto& f : args_) {
                if (direct_ask(f, AssumptionKey::Zero) == true) return false;
                if (direct_ask(f, AssumptionKey::Nonzero) != true) all_nonzero = false;
                if (direct_ask(f, AssumptionKey::Imaginary) == true) {
                    ++imag;
                } else if (direct_ask(f, AssumptionKey::Real) != true) {
                    all_classified = false;
                    break;
                }
            }
            if (!all_classified) return std::nullopt;
            if (imag % 2 == 0) return false;  // even # of i factors → real
            return all_nonzero ? std::optional<bool>{true} : std::nullopt;
        }

        case AssumptionKey::Real:
            if (all_args_have(args_, AssumptionKey::Real, true)) return true;
            // A product of an even number of imaginary factors (rest real) is
            // real: I·I = −1.
            {
                std::size_t imag = 0;
                bool ok = true;
                for (const auto& f : args_) {
                    if (direct_ask(f, AssumptionKey::Imaginary) == true) ++imag;
                    else if (direct_ask(f, AssumptionKey::Real) != true) { ok = false; break; }
                }
                if (ok && imag % 2 == 0) return true;
            }
            return std::nullopt;
        case AssumptionKey::Integer:
            if (all_args_have(args_, AssumptionKey::Integer, true)) return true;
            return std::nullopt;
        case AssumptionKey::Rational:
            if (all_args_have(args_, AssumptionKey::Rational, true)) return true;
            return std::nullopt;
        case AssumptionKey::Finite:
            if (all_args_have(args_, AssumptionKey::Finite, true)) return true;
            return std::nullopt;

        case AssumptionKey::Zero:
            // Any zero factor zeroes the product. (We can't claim the product
            // is *non-zero* without knowing every factor is finite, so that
            // case stays Unknown unless Nonzero handles it below.)
            if (any_arg_has(args_, AssumptionKey::Zero, true)) return true;
            if (all_args_have(args_, AssumptionKey::Nonzero, true)
                && all_args_have(args_, AssumptionKey::Finite, true)) return false;
            return std::nullopt;
        case AssumptionKey::Nonzero:
            if (all_args_have(args_, AssumptionKey::Nonzero, true)
                && all_args_have(args_, AssumptionKey::Finite, true)) return true;
            if (any_arg_has(args_, AssumptionKey::Zero, true)) return false;
            return std::nullopt;

        // Sign of product = parity of negative factors, given all are nonzero
        // and real.
        case AssumptionKey::Positive:
        case AssumptionKey::Negative: {
            if (!all_args_have(args_, AssumptionKey::Real, true)) return std::nullopt;
            // Each factor must have a definite sign (positive or negative); both
            // imply nonzero. Asking Positive as well as Negative lets a factor
            // like a² (positive for a ≠ 0, but whose Negative/Nonzero may be
            // unknown) count as a known-positive factor.
            int neg = 0;
            for (const auto& a : args_) {
                if (direct_ask(a, AssumptionKey::Positive) == std::optional<bool>{true}) {
                    continue;
                }
                if (direct_ask(a, AssumptionKey::Negative) == std::optional<bool>{true}) {
                    ++neg;
                    continue;
                }
                return std::nullopt;  // sign (or zeroness) of a factor unknown
            }
            const bool product_negative = (neg % 2) == 1;
            return k == AssumptionKey::Negative ? product_negative
                                                : !product_negative;
        }

        // Sign-direction of a product of real factors. Each factor must be known
        // ≥0 or ≤0; the product's direction is the parity of the ≤0 factors. A
        // provably-zero factor makes the product 0 (both ≥0 and ≤0). Conservative:
        // only the definite direction is reported (the unprovable opposite, which
        // would require ruling out a zero factor, stays nullopt). Closes e.g.
        // positive·nonnegative → nonnegative, nonpositive·nonpositive → nonnegative.
        case AssumptionKey::Nonnegative:
        case AssumptionKey::Nonpositive: {
            if (!all_args_have(args_, AssumptionKey::Real, true)) {
                return std::nullopt;
            }
            int neg_dir = 0;
            for (const auto& a : args_) {
                const bool ge0 =
                    direct_ask(a, AssumptionKey::Nonnegative) == std::optional<bool>{true};
                const bool le0 =
                    direct_ask(a, AssumptionKey::Nonpositive) == std::optional<bool>{true};
                if (ge0 && le0) return true;  // factor is 0 → product is 0
                if (ge0) continue;
                if (le0) { ++neg_dir; continue; }
                return std::nullopt;  // factor's sign direction unknown
            }
            const bool even = (neg_dir % 2) == 0;
            if (k == AssumptionKey::Nonnegative) {
                return even ? std::optional<bool>{true} : std::nullopt;
            }
            return even ? std::nullopt : std::optional<bool>{true};
        }

        // Parity of an integer product: even iff some factor is even; odd iff
        // every factor is odd. Requires all factors to be known integers.
        case AssumptionKey::Even:
            if (!all_args_have(args_, AssumptionKey::Integer, true)) {
                return std::nullopt;
            }
            if (any_arg_has(args_, AssumptionKey::Even, true)) return true;
            // Two consecutive integer factors ⇒ one is even ⇒ product even.
            if (has_consecutive_int_factors(args_)) return true;
            if (all_args_have(args_, AssumptionKey::Odd, true)) return false;
            return std::nullopt;
        case AssumptionKey::Odd:
            if (!all_args_have(args_, AssumptionKey::Integer, true)) {
                return std::nullopt;
            }
            if (all_args_have(args_, AssumptionKey::Odd, true)) return true;
            if (any_arg_has(args_, AssumptionKey::Even, true)) return false;
            // Two consecutive integer factors force an even product, hence not odd.
            if (has_consecutive_int_factors(args_)) return false;
            return std::nullopt;
        case AssumptionKey::Prime:
        case AssumptionKey::Composite:
            // Primality / compositeness of a symbolic product isn't decided
            // structurally.
            return std::nullopt;
        case AssumptionKey::Algebraic:
            // Algebraic numbers are closed under multiplication.
            if (all_args_have(args_, AssumptionKey::Algebraic, true)) return true;
            // Nonzero algebraics times exactly one transcendental ⇒ transcendental.
            if (detail::product_forces_transcendental(args_)) return false;
            return std::nullopt;
        case AssumptionKey::Transcendental:
            if (detail::product_forces_transcendental(args_)) return true;
            return std::nullopt;
        case AssumptionKey::Irrational:
            // A nonzero rational product times exactly one irrational factor is
            // irrational (2·π, √2/2, 3·√2).
            if (detail::product_forces_irrational(args_)) return true;
            return std::nullopt;
        case AssumptionKey::ExtendedReal:
        case AssumptionKey::Infinite:
            // Left to the generic derivation layer.
            return std::nullopt;
    }
    return std::nullopt;
}

std::string Mul::str() const {
    // Print "-a*b" instead of "(-1)*a*b" when the leading factor is -1.
    std::vector<std::string> pieces;
    pieces.reserve(args_.size());
    bool negate = false;
    for (std::size_t i = 0; i < args_.size(); ++i) {
        const auto& a = args_[i];
        if (i == 0 && a->type_id() == TypeId::Integer) {
            const auto& n = static_cast<const Integer&>(*a);
            if (n.value() == -1) {
                negate = true;
                continue;
            }
        }
        std::string s = a->str();
        // Wrap Add subexpressions in parens for readability.
        if (a->type_id() == TypeId::Add) {
            s = "(" + s + ")";
        }
        pieces.push_back(std::move(s));
    }
    std::string out;
    for (std::size_t i = 0; i < pieces.size(); ++i) {
        if (i > 0) out += "*";
        out += pieces[i];
    }
    if (negate) out = "-" + out;
    return out;
}

Expr mul(std::vector<Expr> args) {
    // ---- Step 1: flatten ----
    std::vector<Expr> flat;
    flat.reserve(args.size());
    for (auto& a : args) flatten_into(a, flat);

    // ---- Step 1b: infinity / nan pre-pass ----
    // An infinity factor makes the product an infinity whose sign is the
    // product of the signs of the finite/known factors (oo*pos=oo, oo*neg=-oo).
    // A literal zero factor gives nan (oo*0). zoo absorbs everything nonzero.
    // Unknown-sign symbolic factors coexist (oo*x). Numbers and known-sign
    // symbols are absorbed into the sign.
    if (std::any_of(flat.begin(), flat.end(), [](const Expr& a) {
            return is_infinity(a) || a->type_id() == TypeId::NaN;
        })) {
        bool has_zoo = false, has_zero = false;
        int real_inf = 0;
        int sign = 1;
        std::vector<Expr> keep;
        for (auto& a : flat) {
            switch (a->type_id()) {
                case TypeId::NaN: return S::NaN();
                case TypeId::Infinity: ++real_inf; break;
                case TypeId::NegativeInfinity: ++real_inf; sign = -sign; break;
                case TypeId::ComplexInfinity: has_zoo = true; break;
                default:
                    if (is_number(a)) {
                        const auto& n = static_cast<const Number&>(*a);
                        if (n.is_zero()) has_zero = true;
                        else if (n.is_negative()) sign = -sign;
                    } else if (is_negative(a) == true) {
                        sign = -sign;  // known-negative symbol absorbed
                    } else if (is_positive(a) == true) {
                        // known-positive symbol absorbed (sign unchanged)
                    } else {
                        keep.push_back(a);  // unknown sign coexists (oo*x)
                    }
            }
        }
        if (has_zero) return S::NaN();  // 0 * oo (or 0 * zoo) = nan
        Expr chosen = has_zoo ? S::ComplexInfinity()
                      : (sign > 0) ? S::Infinity()
                                   : S::NegativeInfinity();
        (void)real_inf;
        if (keep.empty()) return chosen;
        keep.push_back(std::move(chosen));
        flat = std::move(keep);
    }

    // ---- Step 2: combine numerics, short-circuit zero ----
    Expr running_prod;
    std::vector<Expr> non_numeric;
    non_numeric.reserve(flat.size());

    for (auto& a : flat) {
        if (is_number(a)) {
            if (static_cast<const Number&>(*a).is_zero()) {
                return S::Zero();
            }
            if (!running_prod) {
                running_prod = a;
            } else {
                auto combined = number_mul(static_cast<const Number&>(*running_prod),
                                           static_cast<const Number&>(*a));
                if (combined) {
                    running_prod = *combined;
                } else {
                    non_numeric.push_back(std::move(a));
                }
            }
        } else {
            non_numeric.push_back(std::move(a));
        }
    }

    // ---- Step 2b: fold ImaginaryUnit pairs (I * I = -1) ----
    // We do this before base collection so I^n in Pow doesn't have to
    // handle exponent summation; instead we count occurrences directly.
    {
        const Expr& I_atom = S::I();
        int i_count = 0;
        std::vector<Expr> filtered;
        filtered.reserve(non_numeric.size());
        for (auto& e : non_numeric) {
            if (e == I_atom) ++i_count;
            else filtered.push_back(std::move(e));
        }
        non_numeric = std::move(filtered);
        if (i_count > 0) {
            const int neg_pairs = i_count / 2;
            const bool extra_I = (i_count % 2) == 1;
            if (neg_pairs % 2 == 1) {
                // multiply running_prod by -1
                if (!running_prod) {
                    running_prod = S::NegativeOne();
                } else {
                    auto neg = number_mul(
                        static_cast<const Number&>(*running_prod),
                        static_cast<const Number&>(*S::NegativeOne()));
                    if (neg) running_prod = *neg;
                }
            }
            if (extra_I) non_numeric.push_back(I_atom);
        }
    }

    // ---- Step 3: base collection — group same base, sum exponents ----
    detail::ExprMap<Expr> base_to_exp;
    std::vector<Expr> base_order;
    base_order.reserve(non_numeric.size());

    for (auto& a : non_numeric) {
        auto be = as_base_exp(a);
        auto it = base_to_exp.find(be.base);
        if (it == base_to_exp.end()) {
            base_to_exp.emplace(be.base, be.exp);
            base_order.push_back(be.base);
        } else {
            // Sum the exponents — uses add() which canonicalizes.
            it->second = add(it->second, be.exp);
        }
    }

    // ---- Step 4: rebuild non-numeric args from the collected dict ----
    std::vector<Expr> rest;
    rest.reserve(base_order.size());
    for (const auto& base : base_order) {
        const auto& exp = base_to_exp.at(base);
        // pow() handles its own auto-eval (exp == 0 → 1, exp == 1 → base, etc.)
        auto p = pow(base, exp);
        if (p == S::One()) continue;  // x^0 dropped
        rest.push_back(std::move(p));
    }

    // ---- Step 4b: fold numeric results from base collection back into the
    // running product. pow(base, exp) in Step 4 can collapse a numeric base to
    // a Number (2^(1/2)·2^(1/2) → pow(2,1) = 2) or to a numeric·radical Mul
    // (2^(3/2) → 2·√2). Those numeric parts must merge into running_prod, or
    // they survive as un-collapsed factors (√2·√8 → 2·2 instead of 4, and the
    // analogous 1/2·2 in cot(π/4)).
    {
        std::vector<Expr> swept;
        swept.reserve(rest.size());
        auto fold = [&](const Expr& n) {
            if (!running_prod) {
                running_prod = n;
                return;
            }
            auto c = number_mul(static_cast<const Number&>(*running_prod),
                                static_cast<const Number&>(*n));
            if (c) running_prod = *c; else swept.push_back(n);
        };
        for (auto& p : rest) {
            if (is_number(p)) {
                fold(p);
            } else if (p->type_id() == TypeId::Mul) {
                for (const auto& f : p->args()) {
                    if (is_number(f)) fold(f);
                    else swept.push_back(f);
                }
            } else {
                swept.push_back(std::move(p));
            }
        }
        rest = std::move(swept);
    }

    // ---- Step 5: drop one factor, sort ----
    std::sort(rest.begin(), rest.end(), detail::canonical_less);

    if (running_prod && static_cast<const Number&>(*running_prod).is_one()) {
        running_prod.reset();
    }

    // ---- Step 6: assemble — number first per SymPy display convention ----
    std::vector<Expr> out;
    out.reserve(rest.size() + 1);
    if (running_prod) out.push_back(std::move(running_prod));
    for (auto& a : rest) out.push_back(std::move(a));

    if (out.empty()) return S::One();
    if (out.size() == 1) return std::move(out[0]);
    return make<Mul>(std::move(out));
}

Expr mul(const Expr& a, const Expr& b) {
    return mul(std::vector<Expr>{a, b});
}

}  // namespace sympp
