#include <sympp/core/piecewise.hpp>

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include <sympp/core/basic.hpp>
#include <sympp/core/boolean.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

constexpr std::size_t kPiecewiseHashSeed = 0xCAFE'BABE'1234'5678ULL;

[[nodiscard]] std::size_t hash_combine(std::size_t a, std::size_t b) noexcept {
    return a ^ (b + 0x9E37'79B9'7F4A'7C15ULL + (a << 6) + (a >> 2));
}

}  // namespace

Piecewise::Piecewise(std::vector<Expr> flat_args)
    : args_(std::move(flat_args)) {
    std::size_t h = kPiecewiseHashSeed;
    for (const auto& a : args_) h = hash_combine(h, a->hash());
    hash_ = h;
}

bool Piecewise::equals(const Basic& other) const noexcept {
    if (other.type_id() != TypeId::Piecewise) return false;
    const auto& o = static_cast<const Piecewise&>(other);
    if (args_.size() != o.args_.size()) return false;
    for (std::size_t i = 0; i < args_.size(); ++i) {
        if (!(args_[i] == o.args_[i])) return false;
    }
    return true;
}

std::string Piecewise::str() const {
    std::string out = "Piecewise(";
    for (std::size_t i = 0; i + 1 < args_.size(); i += 2) {
        if (i > 0) out += ", ";
        out += "(" + args_[i]->str() + ", " + args_[i + 1]->str() + ")";
    }
    out += ")";
    return out;
}

Expr piecewise(std::vector<PiecewiseBranch> branches) {
    std::vector<Expr> flat;
    flat.reserve(branches.size() * 2);
    for (auto& b : branches) {
        if (is_boolean_false(b.cond)) {
            continue;  // drop
        }
        if (is_boolean_true(b.cond)) {
            // First True condition: branch wins. Push it and stop — anything
            // after is unreachable.
            flat.push_back(std::move(b.value));
            flat.push_back(std::move(b.cond));
            break;
        }
        flat.push_back(std::move(b.value));
        flat.push_back(std::move(b.cond));
    }

    // No surviving branch — SymPy returns Undefined; we surface an empty
    // Piecewise that can stay symbolic.
    if (flat.empty()) {
        return make<Piecewise>(std::move(flat));
    }

    // Single branch with True condition: return the value directly.
    if (flat.size() == 2 && is_boolean_true(flat[1])) {
        return flat[0];
    }

    return make<Piecewise>(std::move(flat));
}

}  // namespace sympp
