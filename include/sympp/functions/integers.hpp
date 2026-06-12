#pragma once

// Integer rounding functions: floor, ceiling; floored modulo: mod.
//
// Auto-evaluation:
//   * floor(integer)  -> integer (identity)
//   * ceiling(integer) -> integer (identity)
//   * floor(rational) / ceiling(rational) -> via mpz_fdiv_q / mpz_cdiv_q
//   * floor(float) / ceiling(float)       -> via mpfr_floor / mpfr_ceil
//   * floor / ceiling on a NumberSymbol (Pi, E, ...) — defer to evalf to
//     decide; without simplify support we keep it unevaluated.
//   * mod(p, q) = p - q*floor(p/q) — floored modulo (result takes the sign of
//     q): mod(7,3)=1, mod(-7,3)=2, mod(7,-3)=-2. Integer and rational p,q
//     evaluate; mod(0,q)=0, mod(p,p)=0; q==0 stays unevaluated.
//
// Reference: sympy/functions/elementary/integers.py (Mod is in core/mod.py)

#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class SYMPP_EXPORT Floor final : public Function {
public:
    explicit Floor(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Floor; }
    [[nodiscard]] std::string_view name() const noexcept override { return "floor"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT Ceiling final : public Function {
public:
    explicit Ceiling(Expr arg);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Ceiling; }
    [[nodiscard]] std::string_view name() const noexcept override { return "ceiling"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

class SYMPP_EXPORT Mod final : public Function {
public:
    Mod(Expr p, Expr q);
    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Mod; }
    [[nodiscard]] std::string_view name() const noexcept override { return "Mod"; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] std::optional<bool> ask(AssumptionKey k) const noexcept override;
};

[[nodiscard]] SYMPP_EXPORT Expr floor(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr ceiling(const Expr& arg);
[[nodiscard]] SYMPP_EXPORT Expr mod(const Expr& p, const Expr& q);

}  // namespace sympp
