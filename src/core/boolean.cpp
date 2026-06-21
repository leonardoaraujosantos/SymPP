#include <sympp/core/boolean.hpp>

#include <algorithm>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include <gmpxx.h>
#include <mpfr.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/float.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/number.hpp>
#include <sympp/core/number_arith.hpp>
#include <sympp/core/queries.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/core/type_id.hpp>

namespace sympp {

namespace {

constexpr std::size_t kBoolTrueHash = 0xB001'B001'C0DE'1ULL;
constexpr std::size_t kBoolFalseHash = 0xB001'B001'C0DE'0ULL;
constexpr std::size_t kRelHashSeed = 0xBEEF'CAFE'1234'5678ULL;

[[nodiscard]] std::size_t hash_combine(std::size_t a, std::size_t b) noexcept {
    return a ^ (b + 0x9E37'79B9'7F4A'7C15ULL + (a << 6) + (a >> 2));
}

// Returns -1, 0, +1, or std::nullopt when not comparable.
[[nodiscard]] std::optional<int> compare_known(const Expr& a, const Expr& b) {
    if (is_number(a) && is_number(b)) {
        // Promote pair to mpq if both exact, mpfr otherwise.
        const auto& na = static_cast<const Number&>(*a);
        const auto& nb = static_cast<const Number&>(*b);
        if (na.type_id() == TypeId::Float || nb.type_id() == TypeId::Float) {
            int dps = kDefaultDps;
            if (na.type_id() == TypeId::Float) {
                dps = std::max(dps, static_cast<const Float&>(na).precision_dps());
            }
            if (nb.type_id() == TypeId::Float) {
                dps = std::max(dps, static_cast<const Float&>(nb).precision_dps());
            }
            mpfr_t ax, bx;
            mpfr_init2(ax, dps_to_prec(dps));
            mpfr_init2(bx, dps_to_prec(dps));
            if (na.type_id() == TypeId::Float)
                mpfr_set(ax, static_cast<const Float&>(na).value(), MPFR_RNDN);
            else if (na.type_id() == TypeId::Integer)
                mpfr_set_z(ax, static_cast<const Integer&>(na).value().get_mpz_t(), MPFR_RNDN);
            else
                mpfr_set_q(ax, static_cast<const Rational&>(na).value().get_mpq_t(), MPFR_RNDN);
            if (nb.type_id() == TypeId::Float)
                mpfr_set(bx, static_cast<const Float&>(nb).value(), MPFR_RNDN);
            else if (nb.type_id() == TypeId::Integer)
                mpfr_set_z(bx, static_cast<const Integer&>(nb).value().get_mpz_t(), MPFR_RNDN);
            else
                mpfr_set_q(bx, static_cast<const Rational&>(nb).value().get_mpq_t(), MPFR_RNDN);
            int c = mpfr_cmp(ax, bx);
            mpfr_clear(ax);
            mpfr_clear(bx);
            return c < 0 ? -1 : (c > 0 ? 1 : 0);
        }
        mpq_class qa, qb;
        if (na.type_id() == TypeId::Integer) qa = static_cast<const Integer&>(na).value();
        else qa = static_cast<const Rational&>(na).value();
        if (nb.type_id() == TypeId::Integer) qb = static_cast<const Integer&>(nb).value();
        else qb = static_cast<const Rational&>(nb).value();
        int c = cmp(qa, qb);
        return c < 0 ? -1 : (c > 0 ? 1 : 0);
    }
    return std::nullopt;
}

}  // namespace

// ----- BooleanTrue / BooleanFalse -------------------------------------------

std::size_t BooleanTrue::hash() const noexcept { return kBoolTrueHash; }
bool BooleanTrue::equals(const Basic& other) const noexcept {
    return dynamic_cast<const BooleanTrue*>(&other) != nullptr;
}

std::size_t BooleanFalse::hash() const noexcept { return kBoolFalseHash; }
bool BooleanFalse::equals(const Basic& other) const noexcept {
    return dynamic_cast<const BooleanFalse*>(&other) != nullptr;
}

// ----- Relational ------------------------------------------------------------

Relational::Relational(RelKind kind, Expr lhs, Expr rhs)
    : args_{std::move(lhs), std::move(rhs)}, kind_(kind) {
    hash_ = hash_combine(kRelHashSeed,
                        hash_combine(static_cast<std::size_t>(kind),
                                    hash_combine(args_[0]->hash(), args_[1]->hash())));
}

bool Relational::equals(const Basic& other) const noexcept {
    if (other.type_id() != TypeId::Relational) return false;
    const auto& o = static_cast<const Relational&>(other);
    if (kind_ != o.kind_) return false;
    return args_[0] == o.args_[0] && args_[1] == o.args_[1];
}

std::string Relational::str() const {
    const char* op = "?";
    switch (kind_) {
        case RelKind::Eq: op = "=="; break;
        case RelKind::Ne: op = "!="; break;
        case RelKind::Lt: op = "<"; break;
        case RelKind::Le: op = "<="; break;
        case RelKind::Gt: op = ">"; break;
        case RelKind::Ge: op = ">="; break;
    }
    return args_[0]->str() + " " + op + " " + args_[1]->str();
}

// ----- Factories -------------------------------------------------------------

bool is_boolean_true(const Expr& e) noexcept {
    return e && e->type_id() == TypeId::Boolean && e->str() == "True";
}
bool is_boolean_false(const Expr& e) noexcept {
    return e && e->type_id() == TypeId::Boolean && e->str() == "False";
}

namespace {

[[nodiscard]] Expr build_rel(RelKind k, const Expr& a, const Expr& b) {
    // Reflexive / antireflexive shortcut.
    if (a == b) {
        switch (k) {
            case RelKind::Eq:
            case RelKind::Le:
            case RelKind::Ge:
                return S::True();
            case RelKind::Ne:
            case RelKind::Lt:
            case RelKind::Gt:
                return S::False();
        }
    }
    // Numeric comparison short-circuit.
    if (auto c = compare_known(a, b); c.has_value()) {
        bool result = false;
        switch (k) {
            case RelKind::Eq: result = (*c == 0); break;
            case RelKind::Ne: result = (*c != 0); break;
            case RelKind::Lt: result = (*c < 0); break;
            case RelKind::Le: result = (*c <= 0); break;
            case RelKind::Gt: result = (*c > 0); break;
            case RelKind::Ge: result = (*c >= 0); break;
        }
        return result ? S::True() : S::False();
    }
    return make<Relational>(k, a, b);
}

}  // namespace

// ----- Logical connectives ---------------------------------------------------

bool is_bool_not(const Expr& e) noexcept {
    return e && dynamic_cast<const BoolNot*>(e.get()) != nullptr;
}
bool is_bool_and(const Expr& e) noexcept {
    return e && dynamic_cast<const BoolAnd*>(e.get()) != nullptr;
}
bool is_bool_or(const Expr& e) noexcept {
    return e && dynamic_cast<const BoolOr*>(e.get()) != nullptr;
}

namespace {
[[nodiscard]] std::string paren_if_compound(const Expr& e) {
    if (is_bool_and(e) || is_bool_or(e)) return "(" + e->str() + ")";
    return e->str();
}
[[nodiscard]] std::string join_ops(const std::vector<Expr>& xs, const char* sep) {
    std::string s;
    for (std::size_t i = 0; i < xs.size(); ++i) {
        if (i) s += sep;
        s += paren_if_compound(xs[i]);
    }
    return s;
}
}  // namespace

BoolNot::BoolNot(Expr arg) : args_{std::move(arg)} {
    hash_ = hash_combine(0xC0DE'0741ULL, args_[0]->hash());
}
bool BoolNot::equals(const Basic& other) const noexcept {
    const auto* p = dynamic_cast<const BoolNot*>(&other);
    return p && args_[0]->equals(*p->args_[0]);
}
std::string BoolNot::str() const { return "~" + paren_if_compound(args_[0]); }

BoolAnd::BoolAnd(std::vector<Expr> args) : args_(std::move(args)) {
    std::size_t h = 0xA0DD'0011ULL;
    for (const auto& x : args_) h = hash_combine(h, x->hash());
    hash_ = h;
}
bool BoolAnd::equals(const Basic& other) const noexcept {
    const auto* p = dynamic_cast<const BoolAnd*>(&other);
    if (!p || p->args_.size() != args_.size()) return false;
    for (std::size_t i = 0; i < args_.size(); ++i) {
        if (!args_[i]->equals(*p->args_[i])) return false;
    }
    return true;
}
std::string BoolAnd::str() const { return join_ops(args_, " & "); }

BoolOr::BoolOr(std::vector<Expr> args) : args_(std::move(args)) {
    std::size_t h = 0x0011'0022ULL;
    for (const auto& x : args_) h = hash_combine(h, x->hash());
    hash_ = h;
}
bool BoolOr::equals(const Basic& other) const noexcept {
    const auto* p = dynamic_cast<const BoolOr*>(&other);
    if (!p || p->args_.size() != args_.size()) return false;
    for (std::size_t i = 0; i < args_.size(); ++i) {
        if (!args_[i]->equals(*p->args_[i])) return false;
    }
    return true;
}
std::string BoolOr::str() const { return join_ops(args_, " | "); }

namespace {
// Sort + de-duplicate operands by their printed form (canonical, stable).
void canonicalize(std::vector<Expr>& xs) {
    std::sort(xs.begin(), xs.end(),
              [](const Expr& a, const Expr& b) { return a->str() < b->str(); });
    xs.erase(std::unique(xs.begin(), xs.end(),
                         [](const Expr& a, const Expr& b) { return a->equals(*b); }),
             xs.end());
}
[[nodiscard]] bool contains_str(const std::vector<Expr>& xs, const std::string& s) {
    for (const auto& x : xs) if (x->str() == s) return true;
    return false;
}
}  // namespace

Expr bool_not(const Expr& a) {
    if (is_boolean_true(a)) return S::False();
    if (is_boolean_false(a)) return S::True();
    if (is_bool_not(a)) return static_cast<const BoolNot&>(*a).operand();
    return make<BoolNot>(a);
}

Expr bool_and(std::vector<Expr> args) {
    std::vector<Expr> flat;
    for (auto& a : args) {
        if (is_boolean_true(a)) continue;            // identity
        if (is_boolean_false(a)) return S::False();  // annihilator
        if (is_bool_and(a)) {
            for (const auto& sub : a->args()) flat.push_back(sub);
        } else {
            flat.push_back(a);
        }
    }
    canonicalize(flat);
    for (const auto& x : flat) {                      // x ∧ ¬x = False
        if (contains_str(flat, bool_not(x)->str())) return S::False();
    }
    if (flat.empty()) return S::True();
    if (flat.size() == 1) return flat[0];
    return make<BoolAnd>(std::move(flat));
}
Expr bool_and(const Expr& a, const Expr& b) { return bool_and(std::vector<Expr>{a, b}); }

Expr bool_or(std::vector<Expr> args) {
    std::vector<Expr> flat;
    for (auto& a : args) {
        if (is_boolean_false(a)) continue;          // identity
        if (is_boolean_true(a)) return S::True();   // annihilator
        if (is_bool_or(a)) {
            for (const auto& sub : a->args()) flat.push_back(sub);
        } else {
            flat.push_back(a);
        }
    }
    canonicalize(flat);
    for (const auto& x : flat) {                     // x ∨ ¬x = True
        if (contains_str(flat, bool_not(x)->str())) return S::True();
    }
    if (flat.empty()) return S::False();
    if (flat.size() == 1) return flat[0];
    return make<BoolOr>(std::move(flat));
}
Expr bool_or(const Expr& a, const Expr& b) { return bool_or(std::vector<Expr>{a, b}); }

Expr bool_xor(const Expr& a, const Expr& b) {
    return bool_or(bool_and(a, bool_not(b)), bool_and(bool_not(a), b));
}
Expr implies(const Expr& a, const Expr& b) { return bool_or(bool_not(a), b); }
Expr equivalent(const Expr& a, const Expr& b) {
    return bool_or(bool_and(a, b), bool_and(bool_not(a), bool_not(b)));
}

Expr eq(const Expr& a, const Expr& b) { return build_rel(RelKind::Eq, a, b); }
Expr ne(const Expr& a, const Expr& b) { return build_rel(RelKind::Ne, a, b); }
Expr lt(const Expr& a, const Expr& b) { return build_rel(RelKind::Lt, a, b); }
Expr le(const Expr& a, const Expr& b) { return build_rel(RelKind::Le, a, b); }
Expr gt(const Expr& a, const Expr& b) { return build_rel(RelKind::Gt, a, b); }
Expr ge(const Expr& a, const Expr& b) { return build_rel(RelKind::Ge, a, b); }

}  // namespace sympp
