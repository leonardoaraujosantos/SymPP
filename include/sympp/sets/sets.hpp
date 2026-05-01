#pragma once

// Set algebra — Interval, FiniteSet, Union, Intersection, Complement,
// EmptySet, Reals. Foundation for solveset return types.
//
// The set hierarchy lives outside the Basic Expr tree because membership
// queries and operations like Union/Intersection are first-class
// algorithms over set objects, not expression rewrites.
//
// Reference: sympy/sets/{sets,fancysets}.py

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class Set;
using SetPtr = std::shared_ptr<const Set>;

enum class SetKind : std::uint8_t {
    Empty,
    Reals,
    Integers,
    Interval,
    FiniteSet,
    Union,
    Intersection,
    Complement,
    ConditionSet,   // {x ∈ S : φ(x)}
    ImageSet,       // {f(n) : n ∈ S}
};

class SYMPP_EXPORT Set {
public:
    virtual ~Set() = default;
    [[nodiscard]] virtual SetKind kind() const noexcept = 0;
    [[nodiscard]] virtual std::string str() const = 0;
    // contains(elem) — best-effort. Returns std::nullopt for symbolic
    // elements where membership can't be decided.
    [[nodiscard]] virtual std::optional<bool> contains(const Expr& elem) const = 0;
};

// Empty set.
class SYMPP_EXPORT EmptySet final : public Set {
public:
    [[nodiscard]] SetKind kind() const noexcept override { return SetKind::Empty; }
    [[nodiscard]] std::string str() const override { return "EmptySet"; }
    [[nodiscard]] std::optional<bool> contains(const Expr&) const override {
        return false;
    }
};

// All real numbers.
class SYMPP_EXPORT Reals final : public Set {
public:
    [[nodiscard]] SetKind kind() const noexcept override { return SetKind::Reals; }
    [[nodiscard]] std::string str() const override { return "Reals"; }
    [[nodiscard]] std::optional<bool> contains(const Expr& e) const override;
};

// All integers.
class SYMPP_EXPORT Integers final : public Set {
public:
    [[nodiscard]] SetKind kind() const noexcept override {
        return SetKind::Integers;
    }
    [[nodiscard]] std::string str() const override { return "Integers"; }
    [[nodiscard]] std::optional<bool> contains(const Expr& e) const override;
};

// [a, b] / (a, b) / [a, b) / (a, b].
class SYMPP_EXPORT Interval final : public Set {
public:
    Interval(Expr lo, Expr hi, bool left_open = false, bool right_open = false);
    [[nodiscard]] SetKind kind() const noexcept override { return SetKind::Interval; }
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] std::optional<bool> contains(const Expr& e) const override;
    [[nodiscard]] const Expr& lo() const { return lo_; }
    [[nodiscard]] const Expr& hi() const { return hi_; }
    [[nodiscard]] bool left_open() const noexcept { return left_open_; }
    [[nodiscard]] bool right_open() const noexcept { return right_open_; }

private:
    Expr lo_, hi_;
    bool left_open_, right_open_;
};

// Finite set of explicit elements.
class SYMPP_EXPORT FiniteSet final : public Set {
public:
    explicit FiniteSet(std::vector<Expr> elements);
    [[nodiscard]] SetKind kind() const noexcept override { return SetKind::FiniteSet; }
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] std::optional<bool> contains(const Expr& e) const override;
    [[nodiscard]] const std::vector<Expr>& elements() const { return elements_; }
    [[nodiscard]] std::size_t size() const noexcept { return elements_.size(); }

private:
    std::vector<Expr> elements_;
};

class SYMPP_EXPORT Union final : public Set {
public:
    Union(SetPtr a, SetPtr b);
    [[nodiscard]] SetKind kind() const noexcept override { return SetKind::Union; }
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] std::optional<bool> contains(const Expr& e) const override;
    [[nodiscard]] const SetPtr& lhs() const { return a_; }
    [[nodiscard]] const SetPtr& rhs() const { return b_; }

private:
    SetPtr a_, b_;
};

class SYMPP_EXPORT Intersection final : public Set {
public:
    Intersection(SetPtr a, SetPtr b);
    [[nodiscard]] SetKind kind() const noexcept override {
        return SetKind::Intersection;
    }
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] std::optional<bool> contains(const Expr& e) const override;
    [[nodiscard]] const SetPtr& lhs() const { return a_; }
    [[nodiscard]] const SetPtr& rhs() const { return b_; }

private:
    SetPtr a_, b_;
};

class SYMPP_EXPORT Complement final : public Set {
public:
    Complement(SetPtr universal, SetPtr removed);
    [[nodiscard]] SetKind kind() const noexcept override {
        return SetKind::Complement;
    }
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] std::optional<bool> contains(const Expr& e) const override;
    [[nodiscard]] const SetPtr& universal() const { return a_; }
    [[nodiscard]] const SetPtr& removed() const { return b_; }

private:
    SetPtr a_, b_;
};

// {var ∈ base_set : condition(var)} where condition is an Expr that
// evaluates to a Boolean (or Relational) when given a value for var.
class SYMPP_EXPORT ConditionSet final : public Set {
public:
    ConditionSet(Expr var, Expr condition, SetPtr base_set);
    [[nodiscard]] SetKind kind() const noexcept override {
        return SetKind::ConditionSet;
    }
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] std::optional<bool> contains(const Expr& e) const override;
    [[nodiscard]] const Expr& variable() const { return var_; }
    [[nodiscard]] const Expr& condition() const { return condition_; }
    [[nodiscard]] const SetPtr& base_set() const { return base_; }

private:
    Expr var_;
    Expr condition_;
    SetPtr base_;
};

// {expr(var) : var ∈ base_set} — image of base_set under the map
// var → expr. Used for periodic solution sets like {n·π : n ∈ ℤ}.
class SYMPP_EXPORT ImageSet final : public Set {
public:
    ImageSet(Expr var, Expr expr, SetPtr base_set);
    [[nodiscard]] SetKind kind() const noexcept override {
        return SetKind::ImageSet;
    }
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] std::optional<bool> contains(const Expr& e) const override;
    [[nodiscard]] const Expr& variable() const { return var_; }
    [[nodiscard]] const Expr& image_expr() const { return expr_; }
    [[nodiscard]] const SetPtr& base_set() const { return base_; }

private:
    Expr var_;
    Expr expr_;
    SetPtr base_;
};

// ---- Factories ------------------------------------------------------------

[[nodiscard]] SYMPP_EXPORT SetPtr empty_set();
[[nodiscard]] SYMPP_EXPORT SetPtr reals();
[[nodiscard]] SYMPP_EXPORT SetPtr integers();
[[nodiscard]] SYMPP_EXPORT SetPtr interval(const Expr& lo, const Expr& hi,
                                              bool left_open = false,
                                              bool right_open = false);
[[nodiscard]] SYMPP_EXPORT SetPtr finite_set(std::vector<Expr> elements);
[[nodiscard]] SYMPP_EXPORT SetPtr set_union(SetPtr a, SetPtr b);
[[nodiscard]] SYMPP_EXPORT SetPtr set_intersection(SetPtr a, SetPtr b);
[[nodiscard]] SYMPP_EXPORT SetPtr set_complement(SetPtr universal, SetPtr removed);
[[nodiscard]] SYMPP_EXPORT SetPtr condition_set(const Expr& var,
                                                  const Expr& condition,
                                                  SetPtr base_set);
[[nodiscard]] SYMPP_EXPORT SetPtr image_set(const Expr& var,
                                              const Expr& expr,
                                              SetPtr base_set);

}  // namespace sympp
