#include <sympp/core/mul.hpp>

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
        // Wrap Add and Pow-with-binary-op subexpressions in parens for readability.
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

    // ---- Step 2: combine numerics, short-circuit zero ----
    Expr running_prod;
    std::vector<Expr> rest;
    rest.reserve(flat.size());

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
                    rest.push_back(std::move(a));
                }
            }
        } else {
            rest.push_back(std::move(a));
        }
    }

    // ---- Step 3: drop one factor, sort ----
    std::sort(rest.begin(), rest.end(), detail::canonical_less);

    if (running_prod && static_cast<const Number&>(*running_prod).is_one()) {
        running_prod.reset();
    }

    // ---- Step 4: assemble — number first per SymPy display convention ----
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
