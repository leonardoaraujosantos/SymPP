#include <sympp/sets/sets.hpp>

#include <utility>

#include <sympp/core/basic.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/queries.hpp>

namespace sympp {

// ---- Reals ----------------------------------------------------------------

std::optional<bool> Reals::contains(const Expr& e) const {
    return is_real(e);
}

// ---- Interval -------------------------------------------------------------

Interval::Interval(Expr lo, Expr hi, bool left_open, bool right_open)
    : lo_(std::move(lo)), hi_(std::move(hi)),
      left_open_(left_open), right_open_(right_open) {}

std::string Interval::str() const {
    std::string s;
    s += left_open_ ? "(" : "[";
    s += lo_->str() + ", " + hi_->str();
    s += right_open_ ? ")" : "]";
    return s;
}

std::optional<bool> Interval::contains(const Expr& e) const {
    // Decide membership only when all three values are numeric.
    if (!is_number(e) || !is_number(lo_) || !is_number(hi_)) {
        return std::nullopt;
    }
    Expr lower_diff = e - lo_;       // ≥ 0 (or > 0 if left_open) for membership
    Expr upper_diff = hi_ - e;       // ≥ 0 (or > 0 if right_open) for membership

    auto sign_nonneg = [](const Expr& d, bool strict) -> std::optional<bool> {
        // Numeric e and bounds → e - lo is numeric. Use is_positive /
        // is_negative for a definite answer.
        if (auto p = is_positive(d); p && *p) return true;
        if (auto z = is_zero(d); z && *z) return !strict;
        if (auto n = is_negative(d); n && *n) return false;
        return std::nullopt;
    };

    auto lo_ok = sign_nonneg(lower_diff, left_open_);
    auto hi_ok = sign_nonneg(upper_diff, right_open_);
    if (!lo_ok || !hi_ok) return std::nullopt;
    return *lo_ok && *hi_ok;
}

// ---- FiniteSet ------------------------------------------------------------

FiniteSet::FiniteSet(std::vector<Expr> elements)
    : elements_(std::move(elements)) {}

std::string FiniteSet::str() const {
    std::string s = "{";
    for (std::size_t i = 0; i < elements_.size(); ++i) {
        if (i > 0) s += ", ";
        s += elements_[i]->str();
    }
    s += "}";
    return s;
}

std::optional<bool> FiniteSet::contains(const Expr& e) const {
    bool any_unknown = false;
    for (const auto& m : elements_) {
        if (m == e) return true;
        // Structural inequality could still hide algebraic equality;
        // remain conservative when symbols differ.
        if (m->hash() != e->hash()) continue;
        any_unknown = true;
    }
    if (any_unknown) return std::nullopt;
    return false;
}

// ---- Union ----------------------------------------------------------------

Union::Union(SetPtr a, SetPtr b) : a_(std::move(a)), b_(std::move(b)) {}

std::string Union::str() const {
    return a_->str() + " ∪ " + b_->str();
}

std::optional<bool> Union::contains(const Expr& e) const {
    auto la = a_->contains(e);
    auto lb = b_->contains(e);
    if (la == std::optional<bool>{true} || lb == std::optional<bool>{true}) {
        return true;
    }
    if (la == std::optional<bool>{false} && lb == std::optional<bool>{false}) {
        return false;
    }
    return std::nullopt;
}

// ---- Intersection ---------------------------------------------------------

Intersection::Intersection(SetPtr a, SetPtr b)
    : a_(std::move(a)), b_(std::move(b)) {}

std::string Intersection::str() const {
    return a_->str() + " ∩ " + b_->str();
}

std::optional<bool> Intersection::contains(const Expr& e) const {
    auto la = a_->contains(e);
    auto lb = b_->contains(e);
    if (la == std::optional<bool>{false} || lb == std::optional<bool>{false}) {
        return false;
    }
    if (la == std::optional<bool>{true} && lb == std::optional<bool>{true}) {
        return true;
    }
    return std::nullopt;
}

// ---- Complement -----------------------------------------------------------

Complement::Complement(SetPtr universal, SetPtr removed)
    : a_(std::move(universal)), b_(std::move(removed)) {}

std::string Complement::str() const {
    return a_->str() + " \\ " + b_->str();
}

std::optional<bool> Complement::contains(const Expr& e) const {
    auto la = a_->contains(e);
    auto lb = b_->contains(e);
    if (la == std::optional<bool>{true} && lb == std::optional<bool>{false}) {
        return true;
    }
    if (la == std::optional<bool>{false}
        || lb == std::optional<bool>{true}) {
        return false;
    }
    return std::nullopt;
}

// ---- Factories ------------------------------------------------------------

SetPtr empty_set() { return std::make_shared<EmptySet>(); }
SetPtr reals() { return std::make_shared<Reals>(); }
SetPtr interval(const Expr& lo, const Expr& hi, bool left_open, bool right_open) {
    return std::make_shared<Interval>(lo, hi, left_open, right_open);
}
SetPtr finite_set(std::vector<Expr> elements) {
    if (elements.empty()) return empty_set();
    return std::make_shared<FiniteSet>(std::move(elements));
}
SetPtr set_union(SetPtr a, SetPtr b) {
    if (a->kind() == SetKind::Empty) return b;
    if (b->kind() == SetKind::Empty) return a;
    return std::make_shared<Union>(std::move(a), std::move(b));
}
SetPtr set_intersection(SetPtr a, SetPtr b) {
    if (a->kind() == SetKind::Empty || b->kind() == SetKind::Empty) {
        return empty_set();
    }
    return std::make_shared<Intersection>(std::move(a), std::move(b));
}
SetPtr set_complement(SetPtr universal, SetPtr removed) {
    if (removed->kind() == SetKind::Empty) return universal;
    return std::make_shared<Complement>(std::move(universal), std::move(removed));
}

}  // namespace sympp
