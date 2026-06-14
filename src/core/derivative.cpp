#include <sympp/core/derivative.hpp>

#include <cstddef>
#include <string>
#include <utility>

#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

constexpr std::size_t kDerivativeHashSeed = 0xDE12'7A71'4E00'0000ULL;

[[nodiscard]] std::size_t hash_combine(std::size_t a, std::size_t b) noexcept {
    return a ^ (b + 0x9E37'79B9'7F4A'7C15ULL + (a << 6) + (a >> 2));
}

}  // namespace

Derivative::Derivative(Expr expr, Expr var, std::size_t order)
    : args_{std::move(expr), std::move(var)}, order_(order), hash_(0) {
    std::size_t h = kDerivativeHashSeed;
    h = hash_combine(h, args_[0]->hash());
    h = hash_combine(h, args_[1]->hash());
    h = hash_combine(h, order_);
    hash_ = h;
}

bool Derivative::equals(const Basic& other) const noexcept {
    if (other.type_id() != TypeId::Derivative) return false;
    const auto& o = static_cast<const Derivative&>(other);
    return order_ == o.order_ && args_[0] == o.args_[0] && args_[1] == o.args_[1];
}

std::string Derivative::str() const {
    std::string out = "Derivative(";
    out.append(args_[0]->str());
    out.append(", ");
    if (order_ == 1) {
        out.append(args_[1]->str());
    } else {
        // SymPy renders repeated single-variable derivatives as (var, n).
        out.push_back('(');
        out.append(args_[1]->str());
        out.append(", ");
        out.append(std::to_string(order_));
        out.push_back(')');
    }
    out.push_back(')');
    return out;
}

Expr derivative(const Expr& expr, const Expr& var, std::size_t order) {
    if (order == 0) return expr;
    // Fold d/dvar of an existing same-var Derivative into a single node.
    if (expr->type_id() == TypeId::Derivative) {
        const auto& d = static_cast<const Derivative&>(*expr);
        if (d.var() == var) {
            return make<Derivative>(d.expr(), var, d.order() + order);
        }
    }
    return make<Derivative>(expr, var, order);
}

}  // namespace sympp
