#pragma once

// Derivative — an unevaluated derivative d^n/dvar^n (expr).
//
// Produced when diff() meets a function with no closed-form derivative rule:
// an undefined function f(x), or a special function (besselj, zeta, li, …)
// whose argument-derivative is not tabulated. Returning this node — instead of
// silently collapsing to 0 — keeps the result correct: SymPy does the same
// (`Derivative(f(x), x)`).
//
// Single-variable only: holds (expr, var, order). Higher order prints as
// Derivative(expr, (var, n)) to match SymPy.
//
// Reference: sympy/core/function.py::Derivative

#include <cstddef>
#include <span>
#include <string>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class SYMPP_EXPORT Derivative final : public Basic {
public:
    Derivative(Expr expr, Expr var, std::size_t order);

    [[nodiscard]] TypeId type_id() const noexcept override { return TypeId::Derivative; }
    [[nodiscard]] std::size_t hash() const noexcept override { return hash_; }
    [[nodiscard]] std::span<const Expr> args() const noexcept override { return args_; }
    [[nodiscard]] bool equals(const Basic& other) const noexcept override;
    [[nodiscard]] std::string str() const override;

    [[nodiscard]] const Expr& expr() const noexcept { return args_[0]; }
    [[nodiscard]] const Expr& var() const noexcept { return args_[1]; }
    [[nodiscard]] std::size_t order() const noexcept { return order_; }

private:
    std::vector<Expr> args_;  // {expr, var}
    std::size_t order_;
    std::size_t hash_;
};

// Factory: the unevaluated order-th derivative of `expr` w.r.t. `var`.
// Collapses a derivative-of-a-derivative w.r.t. the same var into one node with
// the orders summed (Derivative(Derivative(f, x), x) -> Derivative(f, (x, 2))).
[[nodiscard]] SYMPP_EXPORT Expr derivative(const Expr& expr, const Expr& var,
                                           std::size_t order = 1);

}  // namespace sympp
