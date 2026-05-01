#include <sympp/core/symbol.hpp>

#include <functional>
#include <string>
#include <utility>

namespace sympp {

namespace {

// Combine the type tag into the name hash so Symbol("x") and a future
// Dummy("x") never collide.
constexpr std::size_t kSymbolHashSeed = 0x9E37'79B9'7F4A'7C15ULL;

[[nodiscard]] std::size_t hash_symbol(const std::string& name) noexcept {
    return std::hash<std::string>{}(name) ^ kSymbolHashSeed;
}

}  // namespace

Symbol::Symbol(std::string name) : name_(std::move(name)), hash_(hash_symbol(name_)) {}

bool Symbol::equals(const Basic& other) const noexcept {
    if (other.type_id() != TypeId::Symbol) return false;
    return name_ == static_cast<const Symbol&>(other).name_;
}

}  // namespace sympp
