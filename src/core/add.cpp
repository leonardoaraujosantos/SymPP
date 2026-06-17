#include <sympp/core/add.hpp>

#include <algorithm>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <sympp/core/basic.hpp>
#include <sympp/core/infinity.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/traversal.hpp>
#include <sympp/core/type_id.hpp>

#include "assumption_helpers.hpp"
#include "canonical_order.hpp"
#include "expr_map.hpp"

namespace sympp {

namespace {

// Recursively flatten nested Adds into a flat arg list.
void flatten_into(const Expr& e, std::vector<Expr>& out) {
    if (e && e->type_id() == TypeId::Add) {
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

constexpr std::size_t kAddHashSeed = 0xADDA'DDAD'1234'5678ULL;

// Decompose `e` into (coefficient, term) for term collection.
//   x       -> (1, x)
//   2*x     -> (2, x)
//   2*x*y   -> (2, Mul(x, y))
//   -x      -> (-1, x)
//   x**2    -> (1, x**2)
// `coeff` is always a Number; `term` carries everything else.
struct CoeffTerm {
    Expr coeff;  // a Number (Integer/Rational/Float)
    Expr term;   // the non-numeric remainder
};

[[nodiscard]] CoeffTerm as_coeff_term(const Expr& e) {
    if (e->type_id() == TypeId::Mul) {
        const auto& args = e->args();
        // The Mul factory places the numeric coefficient first when present.
        if (!args.empty() && is_number(args[0])) {
            if (args.size() == 2) {
                return {args[0], args[1]};
            }
            // Rebuild a Mul of args[1..] as the term.
            std::vector<Expr> rest(args.begin() + 1, args.end());
            return {args[0], make<Mul>(std::move(rest))};
        }
        // Mul without a leading number — coefficient is 1, the whole Mul is the term.
        return {S::One(), e};
    }
    return {S::One(), e};
}

}  // namespace

Add::Add(std::vector<Expr> args) : args_(std::move(args)) {
    compute_hash();
}

void Add::compute_hash() noexcept {
    std::size_t h = kAddHashSeed;
    for (const auto& a : args_) {
        h = hash_combine(h, a->hash());
    }
    hash_ = h;
}

bool Add::equals(const Basic& other) const noexcept {
    if (other.type_id() != TypeId::Add) return false;
    const auto& o = static_cast<const Add&>(other);
    if (args_.size() != o.args_.size()) return false;
    for (std::size_t i = 0; i < args_.size(); ++i) {
        if (!(args_[i] == o.args_[i])) return false;
    }
    return true;
}

std::optional<bool> Add::ask(AssumptionKey k) const noexcept {
    using detail::all_args_have;
    using detail::any_arg_has;
    switch (k) {
        case AssumptionKey::Complex:
            // A sum of finite complex terms is complex.
            if (all_args_have(args_, AssumptionKey::Complex, true)) return true;
            return std::nullopt;
        case AssumptionKey::Imaginary: {
            // A sum is imaginary iff every term is imaginary (i + i = 2i). A mix
            // of a (nonzero) real term and an imaginary term is neither real nor
            // imaginary (1 + i). Zero terms are skipped.
            bool all_imag = true;
            bool any_imag = false;
            bool any_real_nonzero = false;
            for (const auto& a : args_) {
                if (a->ask(AssumptionKey::Zero) == true) continue;
                if (a->ask(AssumptionKey::Imaginary) == true) {
                    any_imag = true;
                } else {
                    all_imag = false;
                    if (a->ask(AssumptionKey::Real) == true
                        && a->ask(AssumptionKey::Nonzero) == true) {
                        any_real_nonzero = true;
                    }
                }
            }
            if (any_imag && all_imag) return true;
            if (any_imag && any_real_nonzero) return false;  // real + imaginary
            return std::nullopt;
        }

        // Closure under sum: if every term is real / integer / rational /
        // finite, the result is too.
        case AssumptionKey::Real:
            if (all_args_have(args_, AssumptionKey::Real, true)) return true;
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

        // Sign propagation. Sum of positives is positive; sum of nonneg with
        // at least one positive is positive (analogously for negative).
        case AssumptionKey::Positive:
            if (all_args_have(args_, AssumptionKey::Positive, true)) return true;
            if (all_args_have(args_, AssumptionKey::Nonnegative, true)
                && any_arg_has(args_, AssumptionKey::Positive, true)) return true;
            return std::nullopt;
        case AssumptionKey::Negative:
            if (all_args_have(args_, AssumptionKey::Negative, true)) return true;
            if (all_args_have(args_, AssumptionKey::Nonpositive, true)
                && any_arg_has(args_, AssumptionKey::Negative, true)) return true;
            return std::nullopt;
        case AssumptionKey::Nonnegative:
            if (all_args_have(args_, AssumptionKey::Nonnegative, true)) return true;
            return std::nullopt;
        case AssumptionKey::Nonpositive:
            if (all_args_have(args_, AssumptionKey::Nonpositive, true)) return true;
            return std::nullopt;

        // Without algebraic analysis we cannot rule out cancellation, so
        // defer Zero / Nonzero queries on Add.
        case AssumptionKey::Zero:
        case AssumptionKey::Nonzero:
            return std::nullopt;

        // Parity of a sum: fold the terms' parities (XOR). Every term must have
        // a known parity, else Unknown.
        case AssumptionKey::Even:
        case AssumptionKey::Odd: {
            int odd_terms = 0;
            for (const auto& a : args_) {
                if (a->ask(AssumptionKey::Even) == true) continue;
                if (a->ask(AssumptionKey::Odd) == true) { ++odd_terms; continue; }
                return std::nullopt;  // a term of unknown parity
            }
            const bool sum_even = (odd_terms % 2) == 0;
            return (k == AssumptionKey::Even) ? sum_even : !sum_even;
        }
        case AssumptionKey::Prime:
        case AssumptionKey::Composite:
            // Primality / compositeness of a symbolic sum isn't decided
            // structurally.
            return std::nullopt;
        case AssumptionKey::Irrational:
            // Left to the generic real ∧ ¬rational derivation.
            return std::nullopt;
    }
    return std::nullopt;
}

std::string Add::str() const {
    std::string out;
    bool first = true;
    for (const auto& a : args_) {
        std::string s = a->str();
        if (first) {
            out = s;
            first = false;
            continue;
        }
        // Pretty-print: write `a - b` instead of `a + -b` when the term has
        // a leading minus sign.
        if (!s.empty() && s[0] == '-') {
            out += " - ";
            out += s.substr(1);
        } else {
            out += " + ";
            out += s;
        }
    }
    return out;
}

Expr add(std::vector<Expr> args) {
    // ---- Step 1: flatten nested Adds ----
    std::vector<Expr> flat;
    flat.reserve(args.size());
    for (auto& a : args) flatten_into(a, flat);

    // ---- Step 1b: infinity / nan pre-pass ----
    // A single oo / -oo / zoo absorbs every finite numeric term; finite numbers
    // are dropped, symbolic terms coexist (oo + x). Contradictions (oo - oo,
    // zoo + oo, zoo + zoo, any nan) collapse to nan. When the result is just the
    // infinity it returns directly; otherwise the reduced list (symbolic
    // terms + one infinity) flows through the normal pipeline below.
    if (std::any_of(flat.begin(), flat.end(), [](const Expr& a) {
            return is_infinity(a) || a->type_id() == TypeId::NaN;
        })) {
        bool pinf = false, ninf = false;
        int zoo_count = 0;
        std::vector<Expr> keep;
        for (auto& a : flat) {
            switch (a->type_id()) {
                case TypeId::NaN: return S::NaN();
                case TypeId::Infinity: pinf = true; break;
                case TypeId::NegativeInfinity: ninf = true; break;
                case TypeId::ComplexInfinity: ++zoo_count; break;
                default:
                    // ±∞ absorbs every finite real term — not just numeric
                    // literals but also closed real constants like √2 or π
                    // (oo + √2 = oo). Symbolic terms (oo + x) are kept, since x
                    // is not known finite.
                    if (!(is_number(a)
                          || (free_symbols(a).empty() && is_real(a) == true))) {
                        keep.push_back(a);
                    }
            }
        }
        if ((pinf && ninf) || (zoo_count > 0 && (pinf || ninf))
            || zoo_count > 1) {
            return S::NaN();
        }
        Expr chosen = zoo_count    ? S::ComplexInfinity()
                      : pinf       ? S::Infinity()
                                   : S::NegativeInfinity();
        if (keep.empty()) return chosen;
        keep.push_back(std::move(chosen));
        flat = std::move(keep);
    }

    // ---- Step 2: separate numerics, accumulate numeric sum ----
    Expr running_sum;
    std::vector<Expr> non_numeric;
    non_numeric.reserve(flat.size());

    for (auto& a : flat) {
        if (is_number(a)) {
            if (!running_sum) {
                running_sum = a;
            } else {
                auto combined = number_add(static_cast<const Number&>(*running_sum),
                                           static_cast<const Number&>(*a));
                if (combined) {
                    running_sum = *combined;
                } else {
                    non_numeric.push_back(std::move(a));
                }
            }
        } else {
            non_numeric.push_back(std::move(a));
        }
    }

    // ---- Step 3: term collection — group same `term`, sum coefficients ----
    // Insertion order is preserved for first-seen terms so the result is
    // deterministic before we sort.
    detail::ExprMap<Expr> term_to_coeff;
    std::vector<Expr> term_order;
    term_order.reserve(non_numeric.size());

    for (auto& a : non_numeric) {
        auto ct = as_coeff_term(a);
        auto it = term_to_coeff.find(ct.term);
        if (it == term_to_coeff.end()) {
            term_to_coeff.emplace(ct.term, ct.coeff);
            term_order.push_back(ct.term);
        } else {
            // Sum the coefficients. Both are Numbers, so number_add must succeed.
            auto sum = number_add(static_cast<const Number&>(*it->second),
                                  static_cast<const Number&>(*ct.coeff));
            if (sum) {
                it->second = *sum;
            }
        }
    }

    // ---- Step 4: rebuild non-numeric args from the collected dict ----
    std::vector<Expr> rest;
    rest.reserve(term_order.size());
    for (const auto& term : term_order) {
        const auto& coeff = term_to_coeff.at(term);
        if (static_cast<const Number&>(*coeff).is_zero()) {
            continue;  // 0 * term drops
        }
        if (static_cast<const Number&>(*coeff).is_one()) {
            rest.push_back(term);
        } else {
            rest.push_back(mul(coeff, term));
        }
    }

    // ---- Step 5: drop a zero numeric, sort non-numerics ----
    std::sort(rest.begin(), rest.end(), detail::canonical_less);

    if (running_sum && static_cast<const Number&>(*running_sum).is_zero()) {
        running_sum.reset();
    }

    // ---- Step 6: assemble — number term last per SymPy display convention ----
    std::vector<Expr> out;
    out.reserve(rest.size() + 1);
    for (auto& a : rest) out.push_back(std::move(a));
    if (running_sum) out.push_back(std::move(running_sum));

    if (out.empty()) return S::Zero();
    if (out.size() == 1) return std::move(out[0]);
    return make<Add>(std::move(out));
}

Expr add(const Expr& a, const Expr& b) {
    return add(std::vector<Expr>{a, b});
}

}  // namespace sympp
