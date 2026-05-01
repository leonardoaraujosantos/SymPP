#include <sympp/core/symbol.hpp>

#include <functional>
#include <string>
#include <utility>

#include <sympp/core/assumption_mask.hpp>

namespace sympp {

namespace {

// Combine the type tag into the name + assumption hash so a Symbol("x") and
// a future Dummy("x") never collide, and so two same-name Symbols with
// different assumption masks hash differently.
constexpr std::size_t kSymbolHashSeed = 0x9E37'79B9'7F4A'7C15ULL;

[[nodiscard]] std::size_t hash_symbol(const std::string& name,
                                      const AssumptionMask& assumptions) noexcept {
    std::size_t h = std::hash<std::string>{}(name) ^ kSymbolHashSeed;
    h ^= assumptions.hash() + 0x9E37'79B9'7F4A'7C15ULL + (h << 6) + (h >> 2);
    return h;
}

}  // namespace

Symbol::Symbol(std::string name)
    : name_(std::move(name)),
      assumptions_(),
      hash_(hash_symbol(name_, assumptions_)) {}

Symbol::Symbol(std::string name, AssumptionMask assumptions)
    : name_(std::move(name)),
      assumptions_(close_assumptions(std::move(assumptions))),
      hash_(hash_symbol(name_, assumptions_)) {}

bool Symbol::equals(const Basic& other) const noexcept {
    if (other.type_id() != TypeId::Symbol) return false;
    const auto& o = static_cast<const Symbol&>(other);
    return name_ == o.name_ && assumptions_ == o.assumptions_;
}

std::optional<bool> Symbol::ask(AssumptionKey k) const noexcept {
    return assumptions_.get(k);
}

}  // namespace sympp
