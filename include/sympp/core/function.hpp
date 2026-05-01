#pragma once

// Function — abstract base for every named function applied to arguments.
//
// Concrete subclasses (sin, cos, exp, log, gamma, ...) inherit this and
// declare their own FunctionId, name, evaluation rules, and chain-rule
// derivative. The default str/equals/hash work via name() + args.
//
// Reference: sympy/core/function.py::Function

#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/function_id.hpp>
#include <sympp/core/type_id.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class SYMPP_EXPORT Function : public Basic {
public:
    [[nodiscard]] TypeId type_id() const noexcept override { return TypeId::Function; }
    [[nodiscard]] std::span<const Expr> args() const noexcept override { return args_; }
    [[nodiscard]] std::size_t hash() const noexcept override { return hash_; }
    [[nodiscard]] bool equals(const Basic& other) const noexcept override;
    [[nodiscard]] std::string str() const override;

    // Subclass identity. Two Function exprs are equal iff their function_id
    // matches and their args are structurally equal.
    [[nodiscard]] virtual FunctionId function_id() const noexcept = 0;

    // Pretty-printed name (e.g. "sin", "log", "Abs"). For Undefined,
    // returns the user-supplied name.
    [[nodiscard]] virtual std::string_view name() const noexcept = 0;

    // Rebuild the same kind of function with new args. Used by traversal /
    // refine to descend through Function nodes.
    [[nodiscard]] virtual Expr rebuild(std::vector<Expr> new_args) const = 0;

    // Partial derivative with respect to the i-th argument, evaluated at the
    // current args. Default implementation returns an unevaluated Derivative
    // marker (a placeholder Expr); overrides supply closed-form rules.
    //
    // Reference: sympy/core/function.py::Function.fdiff (one-based there).
    [[nodiscard]] virtual Expr diff_arg(std::size_t i) const;

protected:
    explicit Function(std::vector<Expr> args);

    void compute_hash(FunctionId fid) noexcept;

    std::vector<Expr> args_;
    std::size_t hash_;
};

}  // namespace sympp
