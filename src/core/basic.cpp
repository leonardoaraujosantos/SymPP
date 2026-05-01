#include <sympp/core/basic.hpp>

#include <span>

namespace sympp {

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

}  // namespace sympp
