#include <sympp/core/add.hpp>

#include <algorithm>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <sympp/core/basic.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/mul.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>

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
