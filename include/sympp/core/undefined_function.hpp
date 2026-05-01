#pragma once

// UndefinedFunction — a user-declared anonymous function `f` applied to args.
//
// Two UndefinedFunctions are equal iff they share a name and structurally-
// equal args. A SymPy-style usage:
//
//     auto f = function_symbol("f");          // declare
//     auto e = f(x) + integer(2);             // build expression
//
// The factory `function_symbol(name)` returns a callable adaptor whose
// operator() yields an Expr. This is the closest C++ idiom to SymPy's
// `f = Function('f')` then `f(x)`.

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

class SYMPP_EXPORT UndefinedFunction final : public Function {
public:
    UndefinedFunction(std::string name, std::vector<Expr> args);

    [[nodiscard]] FunctionId function_id() const noexcept override { return FunctionId::Undefined; }
    [[nodiscard]] std::string_view name() const noexcept override { return name_; }
    [[nodiscard]] Expr rebuild(std::vector<Expr> new_args) const override;
    [[nodiscard]] bool equals(const Basic& other) const noexcept override;

private:
    std::string name_;
};

// Lightweight callable bound to a name. Each call builds a fresh
// UndefinedFunction node (interned through the cache).
class SYMPP_EXPORT FunctionSymbol {
public:
    explicit FunctionSymbol(std::string name) : name_(std::move(name)) {}

    [[nodiscard]] const std::string& name() const noexcept { return name_; }

    [[nodiscard]] Expr operator()(Expr x) const;
    [[nodiscard]] Expr operator()(Expr x, Expr y) const;
    [[nodiscard]] Expr operator()(std::vector<Expr> args) const;

private:
    std::string name_;
};

[[nodiscard]] inline FunctionSymbol function_symbol(std::string_view name) {
    return FunctionSymbol{std::string{name}};
}

}  // namespace sympp
