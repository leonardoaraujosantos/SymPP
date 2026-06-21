#include <sympp/logic/logic.hpp>

#include <algorithm>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <sympp/core/basic.hpp>
#include <sympp/core/boolean.hpp>
#include <sympp/core/singletons.hpp>

namespace sympp {

namespace {

using Atoms = std::vector<Expr>;  // distinct propositional variables, by str order
using Index = std::map<std::string, std::size_t>;

void collect(const Expr& e, std::map<std::string, Expr>& seen) {
    if (!e) return;
    if (is_boolean_true(e) || is_boolean_false(e)) return;
    if (is_bool_not(e) || is_bool_and(e) || is_bool_or(e)) {
        for (const auto& a : e->args()) collect(a, seen);
        return;
    }
    seen.emplace(e->str(), e);  // leaf ⇒ propositional variable
}

[[nodiscard]] Atoms atoms_of(const Expr& e, Index& idx) {
    std::map<std::string, Expr> seen;
    collect(e, seen);
    Atoms out;
    for (auto& [s, expr] : seen) {
        idx.emplace(s, out.size());
        out.push_back(expr);
    }
    return out;
}

// Evaluate under an assignment indexed by atom position.
[[nodiscard]] bool eval(const Expr& e, const std::vector<bool>& asg, const Index& idx) {
    if (is_boolean_true(e)) return true;
    if (is_boolean_false(e)) return false;
    if (is_bool_not(e)) return !eval(e->args()[0], asg, idx);
    if (is_bool_and(e)) {
        for (const auto& a : e->args()) if (!eval(a, asg, idx)) return false;
        return true;
    }
    if (is_bool_or(e)) {
        for (const auto& a : e->args()) if (eval(a, asg, idx)) return true;
        return false;
    }
    return asg[idx.at(e->str())];
}

// Literal for atom i with the given polarity.
[[nodiscard]] Expr literal(const Atoms& atoms, std::size_t i, bool positive) {
    return positive ? atoms[i] : bool_not(atoms[i]);
}

// ---- Quine–McCluskey ---------------------------------------------------------

struct Implicant {
    std::uint32_t bits = 0;  // fixed-bit values (meaningful where mask bit is 0)
    std::uint32_t mask = 0;  // 1 ⇒ "don't care" position
    bool operator<(const Implicant& o) const {
        return bits != o.bits ? bits < o.bits : mask < o.mask;
    }
    bool operator==(const Implicant& o) const { return bits == o.bits && mask == o.mask; }
};

[[nodiscard]] bool covers(const Implicant& imp, std::uint32_t minterm) {
    return (minterm & ~imp.mask) == (imp.bits & ~imp.mask);
}

[[nodiscard]] std::vector<Implicant> prime_implicants(std::vector<std::uint32_t> minterms) {
    std::vector<Implicant> current;
    current.reserve(minterms.size());
    for (auto m : minterms) current.push_back({m, 0});
    std::set<Implicant> primes;
    while (!current.empty()) {
        std::vector<bool> used(current.size(), false);
        std::set<Implicant> next;
        for (std::size_t i = 0; i < current.size(); ++i) {
            for (std::size_t j = i + 1; j < current.size(); ++j) {
                if (current[i].mask != current[j].mask) continue;
                std::uint32_t diff = current[i].bits ^ current[j].bits;
                if (diff == 0 || (diff & (diff - 1)) != 0) continue;  // not one bit
                used[i] = used[j] = true;
                next.insert({current[i].bits & ~diff, current[i].mask | diff});
            }
        }
        for (std::size_t i = 0; i < current.size(); ++i) {
            if (!used[i]) primes.insert(current[i]);
        }
        current.assign(next.begin(), next.end());
    }
    return {primes.begin(), primes.end()};
}

// Essential-prime-implicant + greedy minimal cover.
[[nodiscard]] std::vector<Implicant> minimal_cover(const std::vector<Implicant>& primes,
                                                   const std::vector<std::uint32_t>& minterms) {
    std::set<std::uint32_t> uncovered(minterms.begin(), minterms.end());
    std::vector<Implicant> chosen;
    std::vector<bool> picked(primes.size(), false);

    // Essential primes: a minterm covered by exactly one prime.
    for (auto m : minterms) {
        int count = 0;
        std::size_t last = 0;
        for (std::size_t p = 0; p < primes.size(); ++p) {
            if (covers(primes[p], m)) { ++count; last = p; }
        }
        if (count == 1 && !picked[last]) {
            picked[last] = true;
            chosen.push_back(primes[last]);
            for (auto it = uncovered.begin(); it != uncovered.end();) {
                it = covers(primes[last], *it) ? uncovered.erase(it) : std::next(it);
            }
        }
    }
    // Greedy for the rest.
    while (!uncovered.empty()) {
        std::size_t best = 0;
        std::size_t best_cov = 0;
        for (std::size_t p = 0; p < primes.size(); ++p) {
            if (picked[p]) continue;
            std::size_t c = 0;
            for (auto m : uncovered) if (covers(primes[p], m)) ++c;
            if (c > best_cov) { best_cov = c; best = p; }
        }
        if (best_cov == 0) break;  // safety
        picked[best] = true;
        chosen.push_back(primes[best]);
        for (auto it = uncovered.begin(); it != uncovered.end();) {
            it = covers(primes[best], *it) ? uncovered.erase(it) : std::next(it);
        }
    }
    return chosen;
}

}  // namespace

std::optional<std::map<std::string, bool>> satisfiable(const Expr& e) {
    Index idx;
    Atoms atoms = atoms_of(e, idx);
    if (atoms.size() > 22) {
        throw std::runtime_error("satisfiable: too many variables (> 22)");
    }
    std::size_t n = atoms.size();
    std::uint64_t total = std::uint64_t{1} << n;
    for (std::uint64_t a = 0; a < total; ++a) {
        std::vector<bool> asg(n);
        for (std::size_t i = 0; i < n; ++i) asg[i] = (a >> i) & 1U;
        if (eval(e, asg, idx)) {
            std::map<std::string, bool> model;
            for (std::size_t i = 0; i < n; ++i) model[atoms[i]->str()] = asg[i];
            return model;
        }
    }
    return std::nullopt;
}

Expr to_dnf(const Expr& e) {
    Index idx;
    Atoms atoms = atoms_of(e, idx);
    std::size_t n = atoms.size();
    if (n > 16) throw std::runtime_error("to_dnf: too many variables (> 16)");
    std::vector<Expr> terms;
    for (std::uint64_t a = 0; a < (std::uint64_t{1} << n); ++a) {
        std::vector<bool> asg(n);
        for (std::size_t i = 0; i < n; ++i) asg[i] = (a >> i) & 1U;
        if (!eval(e, asg, idx)) continue;
        std::vector<Expr> lits;
        for (std::size_t i = 0; i < n; ++i) lits.push_back(literal(atoms, i, asg[i]));
        terms.push_back(bool_and(std::move(lits)));
    }
    if (terms.empty()) return S::False();
    return bool_or(std::move(terms));
}

Expr to_cnf(const Expr& e) {
    Index idx;
    Atoms atoms = atoms_of(e, idx);
    std::size_t n = atoms.size();
    if (n > 16) throw std::runtime_error("to_cnf: too many variables (> 16)");
    std::vector<Expr> clauses;
    for (std::uint64_t a = 0; a < (std::uint64_t{1} << n); ++a) {
        std::vector<bool> asg(n);
        for (std::size_t i = 0; i < n; ++i) asg[i] = (a >> i) & 1U;
        if (eval(e, asg, idx)) continue;  // maxterm where the formula is false
        std::vector<Expr> lits;
        for (std::size_t i = 0; i < n; ++i) lits.push_back(literal(atoms, i, !asg[i]));
        clauses.push_back(bool_or(std::move(lits)));
    }
    if (clauses.empty()) return S::True();
    return bool_and(std::move(clauses));
}

Expr simplify_logic(const Expr& e) {
    Index idx;
    Atoms atoms = atoms_of(e, idx);
    std::size_t n = atoms.size();
    if (n > 16) throw std::runtime_error("simplify_logic: too many variables (> 16)");

    std::vector<std::uint32_t> minterms;
    for (std::uint32_t a = 0; a < (std::uint32_t{1} << n); ++a) {
        std::vector<bool> asg(n);
        for (std::size_t i = 0; i < n; ++i) asg[i] = (a >> i) & 1U;
        if (eval(e, asg, idx)) minterms.push_back(a);
    }
    if (minterms.empty()) return S::False();
    if (minterms.size() == (std::size_t{1} << n)) return S::True();

    auto primes = prime_implicants(minterms);
    auto cover = minimal_cover(primes, minterms);

    std::vector<Expr> terms;
    for (const auto& imp : cover) {
        std::vector<Expr> lits;
        for (std::size_t i = 0; i < n; ++i) {
            if (imp.mask & (std::uint32_t{1} << i)) continue;  // don't care
            lits.push_back(literal(atoms, i, (imp.bits >> i) & 1U));
        }
        terms.push_back(bool_and(std::move(lits)));
    }
    return bool_or(std::move(terms));
}

}  // namespace sympp
