#include <sympp/core/undefined_function.hpp>

#include <cstddef>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include <sympp/core/basic.hpp>
#include <sympp/core/function.hpp>
#include <sympp/core/function_id.hpp>
#include "string_hash.hpp"

namespace sympp {

namespace {

// Mix the user-supplied name into the hash so f(x) and g(x) collide-prevent.
// Function::compute_hash already mixes function_id; we extend with the name
// hash for UndefinedFunction specifically.
[[nodiscard]] std::size_t hash_combine(std::size_t a, std::size_t b) noexcept {
    return a ^ (b + 0x9E37'79B9'7F4A'7C15ULL + (a << 6) + (a >> 2));
}

}  // namespace

UndefinedFunction::UndefinedFunction(std::string name, std::vector<Expr> args)
    : Function(std::move(args)), name_(std::move(name)) {
    compute_hash(FunctionId::Undefined);
    // Mix the name into the existing hash — Function::compute_hash doesn't
    // know about the per-instance name.
    hash_ = hash_combine(hash_, sympp::detail::fnv1a_64(name_));
}

bool UndefinedFunction::equals(const Basic& other) const noexcept {
    if (!Function::equals(other)) return false;
    const auto& o = static_cast<const UndefinedFunction&>(other);
    return name_ == o.name_;
}

Expr UndefinedFunction::rebuild(std::vector<Expr> new_args) const {
    return make<UndefinedFunction>(name_, std::move(new_args));
}

Expr FunctionSymbol::operator()(Expr x) const {
    return make<UndefinedFunction>(name_, std::vector<Expr>{std::move(x)});
}

Expr FunctionSymbol::operator()(Expr x, Expr y) const {
    return make<UndefinedFunction>(name_, std::vector<Expr>{std::move(x), std::move(y)});
}

Expr FunctionSymbol::operator()(Expr x, Expr y, Expr z) const {
    return make<UndefinedFunction>(
        name_, std::vector<Expr>{std::move(x), std::move(y), std::move(z)});
}

Expr FunctionSymbol::operator()(std::vector<Expr> args) const {
    return make<UndefinedFunction>(name_, std::move(args));
}

}  // namespace sympp
