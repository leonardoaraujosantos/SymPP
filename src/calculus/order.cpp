#include <sympp/calculus/order.hpp>

#include <utility>

#include <sympp/core/basic.hpp>
#include <sympp/core/singletons.hpp>

namespace sympp {

namespace {

[[nodiscard]] std::size_t hash_combine(std::size_t a, std::size_t b) noexcept {
    return a ^ (b + 0x9e3779b9 + (a << 6) + (a >> 2));
}

}  // namespace

Order::Order(Expr expr, Expr var, Expr point)
    : args_{std::move(expr), std::move(var), std::move(point)} {
    hash_ = hash_combine(hash_combine(args_[0]->hash(), args_[1]->hash()),
                          args_[2]->hash());
    hash_ = hash_combine(hash_, 0x0BDE'1234'5678'9ABCULL);
}

bool Order::equals(const Basic& other) const noexcept {
    if (other.type_id() != TypeId::Order) return false;
    const auto& o = static_cast<const Order&>(other);
    return args_[0] == o.args_[0] && args_[1] == o.args_[1]
           && args_[2] == o.args_[2];
}

std::string Order::str() const {
    // SymPy print form: O(expr) when (var, point) is the default
    // (point=0); O(expr, (var, point)) otherwise.
    if (args_[2] == S::Zero()) {
        return "O(" + args_[0]->str() + ")";
    }
    return "O(" + args_[0]->str() + ", (" + args_[1]->str() + ", "
           + args_[2]->str() + "))";
}

std::optional<bool> Order::ask(AssumptionKey k) const noexcept {
    // Asymptotic remainders carry no definite finite-or-not answer in
    // general; report Unknown for everything except finiteness which is
    // also undecidable here.
    (void)k;
    return std::nullopt;
}

Expr order(const Expr& expr, const Expr& var) {
    return order(expr, var, S::Zero());
}

Expr order(const Expr& expr, const Expr& var, const Expr& point) {
    // O(0) → 0 — the trivial absorbing case.
    if (expr == S::Zero()) return S::Zero();
    return make<Order>(expr, var, point);
}

}  // namespace sympp
