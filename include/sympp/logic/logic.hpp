#pragma once

// Propositional-logic algorithms over Boolean expressions built from the
// connectives in <core/boolean.hpp> (bool_and / bool_or / bool_not / bool_xor /
// implies / equivalent) with Symbols, Relationals, and the True/False atoms as
// propositional variables.
//
//   auto a = symbol("a"), b = symbol("b");
//   satisfiable(bool_and(a, bool_not(a)));   // → nullopt (UNSAT)
//   simplify_logic(bool_or(bool_and(a, b), bool_and(a, bool_not(b))));  // → a
//
// Reference: sympy/logic/inference.py (satisfiable),
//            sympy/logic/boolalg.py (to_cnf, to_dnf, simplify_logic).

#include <map>
#include <optional>
#include <string>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

// A satisfying assignment (variable name → truth value), or std::nullopt when
// the formula is unsatisfiable. An always-true formula returns an empty map.
// Decided by exhaustive truth-table search; throws if the formula has more than
// 22 distinct variables.
[[nodiscard]] SYMPP_EXPORT std::optional<std::map<std::string, bool>> satisfiable(
    const Expr& e);

// Conjunctive normal form (AND of ORs) and disjunctive normal form (OR of
// ANDs), as canonical (truth-table-derived) forms equivalent to `e`.
[[nodiscard]] SYMPP_EXPORT Expr to_cnf(const Expr& e);
[[nodiscard]] SYMPP_EXPORT Expr to_dnf(const Expr& e);

// Minimal sum-of-products equivalent to `e`, via Quine–McCluskey (prime
// implicants + essential/greedy cover). Throws above 16 variables.
[[nodiscard]] SYMPP_EXPORT Expr simplify_logic(const Expr& e);

}  // namespace sympp
