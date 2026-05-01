#include <sympp/core/basic.hpp>

#include <span>
#include <utility>

#include <sympp/core/expr_cache.hpp>

namespace sympp {

namespace detail {

Expr intern_through_cache(Expr candidate) {
    return ExprCache::instance().intern(std::move(candidate));
}

}  // namespace detail

bool Basic::equals(const Basic& other) const noexcept {
    if (this == &other) return true;
    if (type_id() != other.type_id()) return false;
    if (hash() != other.hash()) return false;

    auto a = args();
    auto b = other.args();
    if (a.size() != b.size()) return false;

    for (std::size_t i = 0; i < a.size(); ++i) {
        if (!(a[i] == b[i])) return false;
    }
    return true;
}

std::optional<bool> Basic::ask(AssumptionKey /*k*/) const noexcept {
    return std::nullopt;
}

}  // namespace sympp
