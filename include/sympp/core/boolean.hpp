#pragma once

// Minimal Boolean infrastructure for Piecewise conditions and assumption
// reasoning at the expression level.
//
//   * BooleanTrue / BooleanFalse — singleton atoms, accessible via
//     S::True() / S::False().
//   * Relational  — Eq, Ne, Lt, Le, Gt, Ge over two Exprs.
//   * BoolAnd / BoolOr / BoolNot — connectives.
//
// Auto-evaluation:
//   * Eq(x, x)              -> true
//   * Lt(numeric, numeric)  -> true/false by comparison
//   * BoolAnd(true, x)      -> x
//   * BoolAnd(false, _)     -> false
//   * BoolNot(true) -> false; BoolNot(BoolNot(x)) -> x
//
// For now, full SAT-based reasoning is out of scope (per Phase 2c notes).
// These nodes are mainly carriers for Piecewise's branch conditions.

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

class SYMPP_EXPORT BooleanTrue final : public Basic {
public:
    [[nodiscard]] TypeId type_id() const noexcept override { return TypeId::Boolean; }
    [[nodiscard]] std::size_t hash() const noexcept override;
    [[nodiscard]] std::span<const Expr> args() const noexcept override { return {}; }
    [[nodiscard]] bool equals(const Basic& other) const noexcept override;
    [[nodiscard]] std::string str() const override { return "True"; }
};

class SYMPP_EXPORT BooleanFalse final : public Basic {
public:
    [[nodiscard]] TypeId type_id() const noexcept override { return TypeId::Boolean; }
    [[nodiscard]] std::size_t hash() const noexcept override;
    [[nodiscard]] std::span<const Expr> args() const noexcept override { return {}; }
    [[nodiscard]] bool equals(const Basic& other) const noexcept override;
    [[nodiscard]] std::string str() const override { return "False"; }
};

enum class RelKind : std::uint8_t { Eq, Ne, Lt, Le, Gt, Ge };

class SYMPP_EXPORT Relational final : public Basic {
public:
    Relational(RelKind kind, Expr lhs, Expr rhs);

    [[nodiscard]] TypeId type_id() const noexcept override { return TypeId::Relational; }
    [[nodiscard]] std::size_t hash() const noexcept override { return hash_; }
    [[nodiscard]] std::span<const Expr> args() const noexcept override { return args_; }
    [[nodiscard]] bool equals(const Basic& other) const noexcept override;
    [[nodiscard]] std::string str() const override;

    [[nodiscard]] RelKind kind() const noexcept { return kind_; }
    [[nodiscard]] const Expr& lhs() const noexcept { return args_[0]; }
    [[nodiscard]] const Expr& rhs() const noexcept { return args_[1]; }

private:
    Expr args_[2];
    RelKind kind_;
    std::size_t hash_;
};

// ----- Factories -------------------------------------------------------------

[[nodiscard]] SYMPP_EXPORT Expr eq(const Expr& a, const Expr& b);
[[nodiscard]] SYMPP_EXPORT Expr ne(const Expr& a, const Expr& b);
[[nodiscard]] SYMPP_EXPORT Expr lt(const Expr& a, const Expr& b);
[[nodiscard]] SYMPP_EXPORT Expr le(const Expr& a, const Expr& b);
[[nodiscard]] SYMPP_EXPORT Expr gt(const Expr& a, const Expr& b);
[[nodiscard]] SYMPP_EXPORT Expr ge(const Expr& a, const Expr& b);

// True if expr is the singleton True or False atom.
[[nodiscard]] SYMPP_EXPORT bool is_boolean_true(const Expr& e) noexcept;
[[nodiscard]] SYMPP_EXPORT bool is_boolean_false(const Expr& e) noexcept;

}  // namespace sympp
