#pragma once

// Infinity atoms — +oo, -oo, zoo (complex/directed infinity) and nan.
//
// Modeled as Basic atoms (not Number subclasses), like ImaginaryUnit: they are
// not finite real numbers and must stay out of the is_number(...) → number
// arithmetic chain. Infinity arithmetic is handled by dedicated passes in
// add()/mul()/pow() (see src/core/infinity.cpp + the add/mul/pow factories).
//
// Reference: sympy/core/numbers.py::Infinity / NegativeInfinity /
//            ComplexInfinity / NaN

#include <cstddef>
#include <optional>
#include <span>
#include <string>

#include <sympp/core/api.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

// +oo — positive real infinity.
class SYMPP_EXPORT Infinity final : public Basic {
public:
    Infinity() = default;
    [[nodiscard]] TypeId type_id() const noexcept override { return TypeId::Infinity; }
    [[nodiscard]] std::size_t hash() const noexcept override;
    [[nodiscard]] std::span<const Expr> args() const noexcept override { return {}; }
    [[nodiscard]] std::string str() const override { return "oo"; }
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

// -oo — negative real infinity.
class SYMPP_EXPORT NegativeInfinity final : public Basic {
public:
    NegativeInfinity() = default;
    [[nodiscard]] TypeId type_id() const noexcept override {
        return TypeId::NegativeInfinity;
    }
    [[nodiscard]] std::size_t hash() const noexcept override;
    [[nodiscard]] std::span<const Expr> args() const noexcept override { return {}; }
    [[nodiscard]] std::string str() const override { return "-oo"; }
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

// zoo — complex (directed) infinity, the value of 1/0.
class SYMPP_EXPORT ComplexInfinity final : public Basic {
public:
    ComplexInfinity() = default;
    [[nodiscard]] TypeId type_id() const noexcept override {
        return TypeId::ComplexInfinity;
    }
    [[nodiscard]] std::size_t hash() const noexcept override;
    [[nodiscard]] std::span<const Expr> args() const noexcept override { return {}; }
    [[nodiscard]] std::string str() const override { return "zoo"; }
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

// nan — the indeterminate result (oo - oo, oo * 0, oo / oo, 1**oo, …).
class SYMPP_EXPORT NaN final : public Basic {
public:
    NaN() = default;
    [[nodiscard]] TypeId type_id() const noexcept override { return TypeId::NaN; }
    [[nodiscard]] std::size_t hash() const noexcept override;
    [[nodiscard]] std::span<const Expr> args() const noexcept override { return {}; }
    [[nodiscard]] std::string str() const override { return "nan"; }
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

// True if `e` is one of the infinity atoms (oo, -oo, zoo).
[[nodiscard]] SYMPP_EXPORT bool is_infinity(const Expr& e) noexcept;

// Evaluate base^exp when an infinity (or nan) is involved. Returns the folded
// value (oo / -oo / zoo / 0 / 1 / nan) for the determinate cases, std::nullopt
// when no infinity is involved or the result is left unevaluated.
[[nodiscard]] SYMPP_EXPORT std::optional<Expr> pow_with_infinity(
    const Expr& base, const Expr& exp);

}  // namespace sympp
