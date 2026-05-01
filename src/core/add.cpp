#include <sympp/core/add.hpp>

#include <algorithm>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <sympp/core/basic.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>

#include "canonical_order.hpp"

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

// Hash combine
[[nodiscard]] std::size_t hash_combine(std::size_t a, std::size_t b) noexcept {
    return a ^ (b + 0x9E37'79B9'7F4A'7C15ULL + (a << 6) + (a >> 2));
}

constexpr std::size_t kAddHashSeed = 0xADDA'DDAD'1234'5678ULL;

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
        // Pretty-print: write `a - b` instead of `a + -b` when the term has a
        // leading minus sign.
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

    // ---- Step 2: separate numerics from the rest, accumulate numeric sum ----
    Expr running_sum;  // null until we encounter the first numeric arg
    std::vector<Expr> rest;
    rest.reserve(flat.size());

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
                    rest.push_back(std::move(a));
                }
            }
        } else {
            rest.push_back(std::move(a));
        }
    }

    // ---- Step 3: drop a zero numeric, sort non-numerics ----
    std::sort(rest.begin(), rest.end(), detail::canonical_less);

    if (running_sum && static_cast<const Number&>(*running_sum).is_zero()) {
        running_sum.reset();
    }

    // ---- Step 4: assemble ----
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
