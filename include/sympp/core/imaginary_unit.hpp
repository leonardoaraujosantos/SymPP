#pragma once

// ImaginaryUnit — the singleton I with I*I = -1.
//
// Modeled as a Basic atom (not a Number subclass), since it isn't a real
// number and shouldn't appear in the is_number(...) → is_real chain.
//
// Mul folds I*I to -1 via a dedicated pass (see src/core/mul.cpp).
// Pow folds I^n to {1, I, -1, -I} based on n mod 4 (see src/core/pow.cpp).
//
// Reference: sympy/core/numbers.py::ImaginaryUnit

#include <cstddef>
#include <span>
#include <string>

#include <sympp/core/api.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class SYMPP_EXPORT ImaginaryUnit final : public Basic {
public:
    ImaginaryUnit() = default;

    [[nodiscard]] TypeId type_id() const noexcept override {
        return TypeId::ImaginaryUnit;
    }
    [[nodiscard]] std::size_t hash() const noexcept override;
    [[nodiscard]] std::span<const Expr> args() const noexcept override { return {}; }
    [[nodiscard]] std::string str() const override { return "I"; }
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

}  // namespace sympp
