#include <sympp/core/satask.hpp>

#include <array>
#include <optional>
#include <utility>
#include <vector>

#include <sympp/core/assumption_key.hpp>
#include <sympp/core/assumption_mask.hpp>
#include <sympp/core/basic.hpp>
#include <sympp/core/queries.hpp>

namespace sympp {

Proposition::Proposition(Expr e, AssumptionKey k, bool expected)
    : op_(Op::Lit), expr_(std::move(e)), key_(k), expected_(expected) {}

Proposition Proposition::make(Op op, std::vector<Proposition> args) {
    Proposition p;
    p.op_ = op;
    p.args_ = std::move(args);
    return p;
}

Proposition Q(const Expr& e, AssumptionKey k) { return Proposition(e, k, true); }
Proposition Qn(const Expr& e, AssumptionKey k) { return Proposition(e, k, false); }

Proposition operator!(const Proposition& p) {
    return Proposition::make(Proposition::Op::Not, {p});
}
Proposition operator&&(const Proposition& a, const Proposition& b) {
    return Proposition::make(Proposition::Op::And, {a, b});
}
Proposition operator||(const Proposition& a, const Proposition& b) {
    return Proposition::make(Proposition::Op::Or, {a, b});
}

namespace {

using Op = Proposition::Op;

// ---- Three-valued (Kleene) evaluation ---------------------------------------

[[nodiscard]] std::optional<bool> eval_kleene(const Proposition& p) noexcept {
    switch (p.op()) {
        case Op::Lit: {
            auto v = ask(p.expr(), p.key());
            if (!v.has_value()) return std::nullopt;
            return *v == p.expected();
        }
        case Op::Not: {
            auto v = eval_kleene(p.args().front());
            if (!v.has_value()) return std::nullopt;
            return !*v;
        }
        case Op::And: {
            bool unknown = false;
            for (const auto& a : p.args()) {
                auto v = eval_kleene(a);
                if (v == false) return false;
                if (!v.has_value()) unknown = true;
            }
            return unknown ? std::nullopt : std::optional<bool>(true);
        }
        case Op::Or: {
            bool unknown = false;
            for (const auto& a : p.args()) {
                auto v = eval_kleene(a);
                if (v == true) return true;
                if (!v.has_value()) unknown = true;
            }
            return unknown ? std::nullopt : std::optional<bool>(false);
        }
    }
    return std::nullopt;
}

// ---- Refutation-based entailment for single-expression propositions ---------

struct Literal {
    AssumptionKey key;
    bool expected;
};
using Clause = std::vector<Literal>;  // a conjunction of literals (a DNF term)
using Dnf = std::vector<Clause>;      // a disjunction of conjunctions

// Negation-normal-form negation: pushes ¬ down to the literals so the result
// contains only Lit / And / Or.
[[nodiscard]] Proposition negate_nnf(const Proposition& p) {
    switch (p.op()) {
        case Op::Lit:
            return Proposition(p.expr(), p.key(), !p.expected());
        case Op::Not:
            return p.args().front();  // ¬¬a = a
        case Op::And: {
            std::vector<Proposition> args;
            for (const auto& a : p.args()) args.push_back(negate_nnf(a));
            return Proposition::make(Op::Or, std::move(args));
        }
        case Op::Or: {
            std::vector<Proposition> args;
            for (const auto& a : p.args()) args.push_back(negate_nnf(a));
            return Proposition::make(Op::And, std::move(args));
        }
    }
    return p;
}

// Disjunctive normal form. `Not` nodes are eliminated by recursing through
// negate_nnf, so the input may contain them.
[[nodiscard]] Dnf to_dnf(const Proposition& p) {
    switch (p.op()) {
        case Op::Lit:
            return {{Literal{p.key(), p.expected()}}};
        case Op::Not:
            return to_dnf(negate_nnf(p.args().front()));
        case Op::Or: {
            Dnf out;
            for (const auto& a : p.args()) {
                Dnf d = to_dnf(a);
                out.insert(out.end(), d.begin(), d.end());
            }
            return out;
        }
        case Op::And: {
            Dnf acc{{}};  // identity: one empty clause
            for (const auto& a : p.args()) {
                Dnf d = to_dnf(a);
                Dnf next;
                for (const auto& left : acc) {
                    for (const auto& right : d) {
                        Clause merged = left;
                        merged.insert(merged.end(), right.begin(), right.end());
                        next.push_back(std::move(merged));
                    }
                }
                acc = std::move(next);
            }
            return acc;
        }
    }
    return {};
}

// All predicate keys, so the base mask can be seeded from `ask`.
constexpr std::array<AssumptionKey, 21> kAllKeys{
    AssumptionKey::Real,           AssumptionKey::Rational,
    AssumptionKey::Integer,        AssumptionKey::Positive,
    AssumptionKey::Negative,       AssumptionKey::Zero,
    AssumptionKey::Nonzero,        AssumptionKey::Nonnegative,
    AssumptionKey::Nonpositive,    AssumptionKey::Finite,
    AssumptionKey::Even,           AssumptionKey::Odd,
    AssumptionKey::Complex,        AssumptionKey::Imaginary,
    AssumptionKey::Prime,          AssumptionKey::Composite,
    AssumptionKey::Irrational,     AssumptionKey::Algebraic,
    AssumptionKey::Transcendental, AssumptionKey::ExtendedReal,
    AssumptionKey::Infinite,
};

// Seed an AssumptionMask with everything currently known about `e`.
[[nodiscard]] AssumptionMask base_mask(const Expr& e) noexcept {
    AssumptionMask m;
    for (AssumptionKey k : kAllKeys) {
        if (auto v = ask(e, k); v.has_value()) m.set(k, *v);
    }
    return m;
}

// Is `base ∧ clause` satisfiable? A clause asserts each literal; a literal that
// contradicts a fact already in the mask makes the clause unsatisfiable, and so
// does an inconsistent closure.
[[nodiscard]] bool clause_satisfiable(const AssumptionMask& base, const Clause& clause) {
    AssumptionMask m = base;
    for (const Literal& lit : clause) {
        if (auto cur = m.get(lit.key); cur.has_value() && *cur != lit.expected) {
            return false;
        }
        m.set(lit.key, lit.expected);
    }
    return assumptions_consistent(m);
}

// Is `base ∧ formula` satisfiable? formula in DNF is satisfiable iff some term is.
[[nodiscard]] bool satisfiable(const AssumptionMask& base, const Proposition& formula) {
    for (const Clause& clause : to_dnf(formula)) {
        if (clause_satisfiable(base, clause)) return true;
    }
    return false;
}

// Collect the (single) expression referenced by every literal; returns null if
// the proposition mixes distinct expressions (entailment is single-expr only).
[[nodiscard]] bool single_expr(const Proposition& p, Expr& out) noexcept {
    switch (p.op()) {
        case Op::Lit:
            if (!out) {
                out = p.expr();
                return true;
            }
            return out->equals(*p.expr());
        default:
            for (const auto& a : p.args()) {
                if (!single_expr(a, out)) return false;
            }
            return true;
    }
}

}  // namespace

std::optional<bool> ask(const Proposition& p) noexcept {
    if (auto v = eval_kleene(p); v.has_value()) return v;

    // Upgrade: refute over the single expression the literals share, if any.
    Expr e;
    if (!single_expr(p, e) || !e) return std::nullopt;

    const AssumptionMask base = base_mask(e);
    if (!satisfiable(base, negate_nnf(p))) return true;   // ¬p unsatisfiable ⇒ p
    if (!satisfiable(base, p)) return false;              // p unsatisfiable ⇒ ¬p
    return std::nullopt;
}

}  // namespace sympp
