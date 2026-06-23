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

// ---- Integers --------------------------------------------------------------

std::optional<bool> Integers::contains(const Expr& e) const {
    return is_integer(e);
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
    if (!is_number(e)) return std::nullopt;
    // A ±∞ endpoint is an (always-open) unbounded side: the bound on that side
    // is automatically satisfied for any finite e. The remaining (finite) bound
    // still has to be numeric to decide.
    const bool lo_inf = (lo_ == S::NegativeInfinity());
    const bool hi_inf = (hi_ == S::Infinity());
    if ((!lo_inf && !is_number(lo_)) || (!hi_inf && !is_number(hi_))) {
        return std::nullopt;
    }

    auto sign_nonneg = [](const Expr& d, bool strict) -> std::optional<bool> {
        if (auto p = is_positive(d); p && *p) return true;
        if (auto z = is_zero(d); z && *z) return !strict;
        if (auto n = is_negative(d); n && *n) return false;
        return std::nullopt;
    };

    std::optional<bool> lo_ok =
        lo_inf ? std::optional<bool>{true} : sign_nonneg(e - lo_, left_open_);
    std::optional<bool> hi_ok =
        hi_inf ? std::optional<bool>{true} : sign_nonneg(hi_ - e, right_open_);
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

// ---- ConditionSet ---------------------------------------------------------

ConditionSet::ConditionSet(Expr var, Expr condition, SetPtr base_set)
    : var_(std::move(var)), condition_(std::move(condition)),
      base_(std::move(base_set)) {}

std::string ConditionSet::str() const {
    return "{" + var_->str() + " ∈ " + base_->str()
           + " : " + condition_->str() + "}";
}

std::optional<bool> ConditionSet::contains(const Expr& e) const {
    auto in_base = base_->contains(e);
    if (in_base == std::optional<bool>{false}) return false;
    // We can't symbolically decide whether the condition holds for e
    // without an evaluator over Boolean expressions. Report Unknown.
    (void)in_base;
    return std::nullopt;
}

// ---- ImageSet -------------------------------------------------------------

ImageSet::ImageSet(Expr var, Expr expr, SetPtr base_set)
    : var_(std::move(var)), expr_(std::move(expr)),
      base_(std::move(base_set)) {}

std::string ImageSet::str() const {
    return "ImageSet(Lambda(" + var_->str() + ", " + expr_->str() + "), "
           + base_->str() + ")";
}

std::optional<bool> ImageSet::contains(const Expr& e) const {
    // Membership requires solving expr_(var) = e for var ∈ base_.
    // Beyond the minimal scope; report Unknown.
    (void)e;
    return std::nullopt;
}

// ---- Factories ------------------------------------------------------------

SetPtr empty_set() { return std::make_shared<EmptySet>(); }
SetPtr reals() { return std::make_shared<Reals>(); }
SetPtr complexes() { return std::make_shared<Complexes>(); }
SetPtr integers() { return std::make_shared<Integers>(); }
SetPtr interval(const Expr& lo, const Expr& hi, bool left_open, bool right_open) {
    return std::make_shared<Interval>(lo, hi, left_open, right_open);
}
SetPtr finite_set(std::vector<Expr> elements) {
    if (elements.empty()) return empty_set();
    return std::make_shared<FiniteSet>(std::move(elements));
}
namespace {
// Order two endpoints: 1 if a > b, −1 if a < b, 0 if equal, nullopt if the
// comparison can't be decided (symbolic bounds).
[[nodiscard]] std::optional<int> endpoint_cmp(const Expr& a, const Expr& b) {
    if (a == b) return 0;
    Expr d = a - b;
    if (is_positive(d) == std::optional<bool>{true}) return 1;
    if (is_negative(d) == std::optional<bool>{true}) return -1;
    if (d == S::Zero()) return 0;
    return std::nullopt;
}
}  // namespace

SetPtr set_union(SetPtr a, SetPtr b) {
    if (a->kind() == SetKind::Empty) return b;
    if (b->kind() == SetKind::Empty) return a;
    // Complexes is the universal domain: C ∪ X = C.
    if (a->kind() == SetKind::Complexes) return a;
    if (b->kind() == SetKind::Complexes) return b;
    // Merge two overlapping / adjacent real intervals into one.
    if (a->kind() == SetKind::Interval && b->kind() == SetKind::Interval) {
        const auto& ia = static_cast<const Interval&>(*a);
        const auto& ib = static_cast<const Interval&>(*b);
        // Overlap / touch test: ib.lo ≤ ia.hi and ia.lo ≤ ib.hi.
        auto c1 = endpoint_cmp(ib.lo(), ia.hi());
        auto c2 = endpoint_cmp(ia.lo(), ib.hi());
        if (c1 && c2 && *c1 <= 0 && *c2 <= 0) {
            auto clo = endpoint_cmp(ia.lo(), ib.lo());
            auto chi = endpoint_cmp(ia.hi(), ib.hi());
            if (clo && chi) {
                Expr lo;
                bool lo_open;
                if (*clo < 0) { lo = ia.lo(); lo_open = ia.left_open(); }
                else if (*clo > 0) { lo = ib.lo(); lo_open = ib.left_open(); }
                else { lo = ia.lo(); lo_open = ia.left_open() && ib.left_open(); }
                Expr hi;
                bool hi_open;
                if (*chi > 0) { hi = ia.hi(); hi_open = ia.right_open(); }
                else if (*chi < 0) { hi = ib.hi(); hi_open = ib.right_open(); }
                else { hi = ia.hi(); hi_open = ia.right_open() && ib.right_open(); }
                return interval(lo, hi, lo_open, hi_open);
            }
        }
    }
    return std::make_shared<Union>(std::move(a), std::move(b));
}
SetPtr set_intersection(SetPtr a, SetPtr b) {
    if (a->kind() == SetKind::Empty || b->kind() == SetKind::Empty) {
        return empty_set();
    }
    // Complexes is the universal domain: C ∩ X = X.
    if (a->kind() == SetKind::Complexes) return b;
    if (b->kind() == SetKind::Complexes) return a;
    // Intersect two real intervals: [max(los), min(his)].
    if (a->kind() == SetKind::Interval && b->kind() == SetKind::Interval) {
        const auto& ia = static_cast<const Interval&>(*a);
        const auto& ib = static_cast<const Interval&>(*b);
        auto clo = endpoint_cmp(ia.lo(), ib.lo());
        auto chi = endpoint_cmp(ia.hi(), ib.hi());
        if (clo && chi) {
            Expr lo;
            bool lo_open;
            if (*clo > 0) { lo = ia.lo(); lo_open = ia.left_open(); }
            else if (*clo < 0) { lo = ib.lo(); lo_open = ib.left_open(); }
            else { lo = ia.lo(); lo_open = ia.left_open() || ib.left_open(); }
            Expr hi;
            bool hi_open;
            if (*chi < 0) { hi = ia.hi(); hi_open = ia.right_open(); }
            else if (*chi > 0) { hi = ib.hi(); hi_open = ib.right_open(); }
            else { hi = ia.hi(); hi_open = ia.right_open() || ib.right_open(); }
            if (auto c = endpoint_cmp(lo, hi)) {
                if (*c > 0) return empty_set();
                if (*c == 0) {
                    return (lo_open || hi_open) ? empty_set()
                                               : finite_set({lo});
                }
                return interval(lo, hi, lo_open, hi_open);
            }
        }
    }
    return std::make_shared<Intersection>(std::move(a), std::move(b));
}
SetPtr set_complement(SetPtr universal, SetPtr removed) {
    if (removed->kind() == SetKind::Empty) return universal;
    // ℝ \ [a, b] = (−∞, a) ∪ (b, ∞), with each boundary's open/closed flipped
    // (a point removed from ℝ is excluded from the complement, and vice versa).
    // A ±∞ endpoint drops the ray on that side.
    if (universal->kind() == SetKind::Reals
        && removed->kind() == SetKind::Interval) {
        const auto& iv = static_cast<const Interval&>(*removed);
        const bool a_inf = (iv.lo() == S::NegativeInfinity());
        const bool b_inf = (iv.hi() == S::Infinity());
        SetPtr left = a_inf ? nullptr
                            : interval(S::NegativeInfinity(), iv.lo(), true,
                                       !iv.left_open());
        SetPtr right = b_inf ? nullptr
                             : interval(iv.hi(), S::Infinity(),
                                        !iv.right_open(), true);
        if (!left && !right) return empty_set();  // removed was all of ℝ
        if (!left) return right;
        if (!right) return left;
        return set_union(std::move(left), std::move(right));  // disjoint → Union
    }
    return std::make_shared<Complement>(std::move(universal), std::move(removed));
}

SetPtr condition_set(const Expr& var, const Expr& condition, SetPtr base_set) {
    return std::make_shared<ConditionSet>(var, condition, std::move(base_set));
}

SetPtr image_set(const Expr& var, const Expr& expr, SetPtr base_set) {
    return std::make_shared<ImageSet>(var, expr, std::move(base_set));
}

}  // namespace sympp
