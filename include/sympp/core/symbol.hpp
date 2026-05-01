#pragma once

// Symbol — named symbolic variable. Atomic, immutable.
// Two Symbols with the same name compare equal. Assumptions (Phase 2) will
// extend this so e.g. Symbol("x", real=true) is distinct from Symbol("x").

#include <cstddef>
#include <span>
#include <string>
#include <string_view>

#include <sympp/core/api.hpp>
#include <sympp/core/assumption_key.hpp>
#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

class SYMPP_EXPORT Symbol final : public Basic {
public:
    explicit Symbol(std::string name);
    Symbol(std::string name, AssumptionMask assumptions);

    [[nodiscard]] TypeId type_id() const noexcept override { return TypeId::Symbol; }
    [[nodiscard]] std::size_t hash() const noexcept override { return hash_; }
    [[nodiscard]] std::span<const Expr> args() const noexcept override { return {}; }
    [[nodiscard]] bool equals(const Basic& other) const noexcept override;
    [[nodiscard]] std::string str() const override { return name_; }
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;

    [[nodiscard]] const std::string& name() const noexcept { return name_; }
    [[nodiscard]] const AssumptionMask& assumptions() const noexcept { return assumptions_; }

private:
    std::string name_;
    AssumptionMask assumptions_;
    std::size_t hash_;
};

// Factory. Two Symbols with the same name and same assumption mask are the
// same interned Expr; differing assumptions produce distinct Symbols (matches
// SymPy's `Symbol('x') != Symbol('x', positive=True)` semantics).
[[nodiscard]] inline Expr symbol(std::string_view name) {
    return make<Symbol>(std::string{name});
}
[[nodiscard]] inline Expr symbol(std::string_view name, AssumptionMask assumptions) {
    return make<Symbol>(std::string{name}, std::move(assumptions));
}

}  // namespace sympp
