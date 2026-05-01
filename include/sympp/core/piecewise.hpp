#pragma once

// Piecewise — case-based expression.
//
//     Piecewise((expr1, cond1), (expr2, cond2), ..., (exprN, condN))
//
// On construction:
//   * Drop branches whose condition is False.
//   * If the first branch's condition is True, return its value directly.
//   * If only one branch remains and it's True, return its value.
//   * Otherwise carry a flat args list: [expr1, cond1, expr2, cond2, ...].
//
// Conditions may be Boolean atoms (S::True()/S::False()) or Relational
// expressions. Symbolic conditions stay unevaluated.
//
// Reference: sympy/functions/elementary/piecewise.py

#include <cstddef>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

struct SYMPP_EXPORT PiecewiseBranch {
    Expr value;
    Expr cond;
};

class SYMPP_EXPORT Piecewise final : public Basic {
public:
    explicit Piecewise(std::vector<Expr> flat_args);

    [[nodiscard]] TypeId type_id() const noexcept override { return TypeId::Piecewise; }
    [[nodiscard]] std::size_t hash() const noexcept override { return hash_; }
    [[nodiscard]] std::span<const Expr> args() const noexcept override { return args_; }
    [[nodiscard]] bool equals(const Basic& other) const noexcept override;
    [[nodiscard]] std::string str() const override;

    // Number of branches.
    [[nodiscard]] std::size_t num_branches() const noexcept { return args_.size() / 2; }
    [[nodiscard]] PiecewiseBranch branch(std::size_t i) const noexcept {
        return {args_[2 * i], args_[2 * i + 1]};
    }

private:
    std::vector<Expr> args_;
    std::size_t hash_;
};

[[nodiscard]] SYMPP_EXPORT Expr piecewise(std::vector<PiecewiseBranch> branches);

}  // namespace sympp
