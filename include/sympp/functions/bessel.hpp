#pragma once

// Bessel functions of the first/second kind (J, Y) and the modified Bessel
// functions (I, K), as two-argument special-function classes besselj(ν, z),
// bessely(ν, z), besseli(ν, z), besselk(ν, z).
//
// Auto-evaluation:
//   * besselj(0, 0) = 1, besselj(ν, 0) = 0 for Re(ν) > 0 (likewise besseli)
//   * half-integer orders ±1/2 reduce to elementary closed forms:
//       J_{1/2}(z)  = √(2/πz)·sin z,   J_{−1/2}(z) = √(2/πz)·cos z
//       I_{1/2}(z)  = √(2/πz)·sinh z,  I_{−1/2}(z) = √(2/πz)·cosh z
//       K_{±1/2}(z) = √(π/2z)·e^{−z},  Y_{±1/2} via Y = (J cos νπ − J_{−ν})/sin νπ
//
// Exact z-derivatives use the standard recurrences; the ν-derivative is
// non-elementary and falls back to an unevaluated Derivative.
//
// Reference: sympy/functions/special/bessel.py, DLMF 10.

#include <cstddef>
#include <optional>
#include <string_view>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class SYMPP_EXPORT BesselJ final : public Function {
public:
    BesselJ(Expr nu, Expr z);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::BesselJ; }
    [[nodiscard]] std::string_view name() const noexcept override { return "besselj"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

class SYMPP_EXPORT BesselY final : public Function {
public:
    BesselY(Expr nu, Expr z);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::BesselY; }
    [[nodiscard]] std::string_view name() const noexcept override { return "bessely"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

class SYMPP_EXPORT BesselI final : public Function {
public:
    BesselI(Expr nu, Expr z);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::BesselI; }
    [[nodiscard]] std::string_view name() const noexcept override { return "besseli"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

class SYMPP_EXPORT BesselK final : public Function {
public:
    BesselK(Expr nu, Expr z);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::BesselK; }
    [[nodiscard]] std::string_view name() const noexcept override { return "besselk"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] Expr diff_arg(std::size_t i) const override;
};

[[nodiscard]] SYMPP_EXPORT Expr besselj(const Expr& nu, const Expr& z);
[[nodiscard]] SYMPP_EXPORT Expr bessely(const Expr& nu, const Expr& z);
[[nodiscard]] SYMPP_EXPORT Expr besseli(const Expr& nu, const Expr& z);
[[nodiscard]] SYMPP_EXPORT Expr besselk(const Expr& nu, const Expr& z);

}  // namespace sympp
