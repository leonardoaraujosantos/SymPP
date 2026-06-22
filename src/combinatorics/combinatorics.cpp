#include <sympp/combinatorics/combinatorics.hpp>

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

#include <gmpxx.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/integer.hpp>

namespace sympp::combinatorics {

// ----- Permutation -----------------------------------------------------------

Permutation::Permutation(std::vector<int> image) : image_(std::move(image)) {
    // Validate that image_ is a genuine bijection of {0, …, n−1}.
    std::vector<bool> seen(image_.size(), false);
    for (int v : image_) {
        if (v < 0 || v >= static_cast<int>(image_.size()) || seen[static_cast<std::size_t>(v)]) {
            throw std::invalid_argument("Permutation: image must be a permutation of 0..n-1");
        }
        seen[static_cast<std::size_t>(v)] = true;
    }
}

Permutation Permutation::identity(int n) {
    std::vector<int> id(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) id[static_cast<std::size_t>(i)] = i;
    return Permutation(std::move(id));
}

Permutation Permutation::compose(const Permutation& other) const {
    if (size() != other.size()) {
        throw std::invalid_argument("Permutation::compose: size mismatch");
    }
    std::vector<int> r(image_.size());
    for (int i = 0; i < size(); ++i) r[static_cast<std::size_t>(i)] = apply(other.apply(i));
    return Permutation(std::move(r));
}

Permutation Permutation::inverse() const {
    std::vector<int> r(image_.size());
    for (int i = 0; i < size(); ++i) r[static_cast<std::size_t>(apply(i))] = i;
    return Permutation(std::move(r));
}

std::vector<std::vector<int>> Permutation::cyclic_form() const {
    std::vector<std::vector<int>> cycles;
    std::vector<bool> visited(image_.size(), false);
    for (int i = 0; i < size(); ++i) {
        if (visited[static_cast<std::size_t>(i)]) continue;
        std::vector<int> cycle;
        int j = i;
        while (!visited[static_cast<std::size_t>(j)]) {
            visited[static_cast<std::size_t>(j)] = true;
            cycle.push_back(j);
            j = apply(j);
        }
        if (cycle.size() > 1) cycles.push_back(std::move(cycle));  // drop fixed points
    }
    return cycles;
}

int Permutation::sign() const {
    // Parity = (−1)^(number of transpositions) = (−1)^(Σ (len−1) over cycles).
    long transpositions = 0;
    for (const auto& c : cyclic_form()) transpositions += static_cast<long>(c.size()) - 1;
    return (transpositions % 2 == 0) ? 1 : -1;
}

long Permutation::order() const {
    // lcm of the cycle lengths.
    long ord = 1;
    auto gcd = [](long a, long b) {
        while (b) { long t = a % b; a = b; b = t; }
        return a;
    };
    for (const auto& c : cyclic_form()) {
        long len = static_cast<long>(c.size());
        ord = ord / gcd(ord, len) * len;
    }
    return ord;
}

bool Permutation::is_identity() const {
    for (int i = 0; i < size(); ++i)
        if (apply(i) != i) return false;
    return true;
}

// ----- PermutationGroup ------------------------------------------------------

PermutationGroup::PermutationGroup(std::vector<Permutation> generators)
    : gens_(std::move(generators)), degree_(0) {
    if (gens_.empty()) {
        throw std::invalid_argument("PermutationGroup: at least one generator required");
    }
    degree_ = gens_.front().size();
    for (const auto& g : gens_) {
        if (g.size() != degree_) {
            throw std::invalid_argument("PermutationGroup: generators must share a degree");
        }
    }
}

std::vector<Permutation> PermutationGroup::elements() const {
    // BFS closure over the generators, keyed by array form for set membership.
    std::set<std::vector<int>> seen;
    std::vector<Permutation> result;
    std::vector<Permutation> frontier;

    Permutation id = Permutation::identity(degree_);
    seen.insert(id.array_form());
    result.push_back(id);
    frontier.push_back(id);

    while (!frontier.empty()) {
        std::vector<Permutation> next;
        for (const auto& p : frontier) {
            for (const auto& g : gens_) {
                Permutation q = g.compose(p);
                if (seen.insert(q.array_form()).second) {
                    result.push_back(q);
                    next.push_back(q);
                }
            }
        }
        frontier = std::move(next);
    }
    return result;
}

long PermutationGroup::order() const { return static_cast<long>(elements().size()); }

bool PermutationGroup::contains(const Permutation& p) const {
    if (p.size() != degree_) return false;
    for (const auto& e : elements())
        if (e == p) return true;
    return false;
}

bool PermutationGroup::is_abelian() const {
    // Abelian iff every pair of generators commutes.
    for (std::size_t i = 0; i < gens_.size(); ++i) {
        for (std::size_t j = i + 1; j < gens_.size(); ++j) {
            if (gens_[i].compose(gens_[j]) != gens_[j].compose(gens_[i])) return false;
        }
    }
    return true;
}

// ----- Standard groups -------------------------------------------------------

PermutationGroup symmetric_group(int n) {
    if (n < 1) throw std::invalid_argument("symmetric_group: n >= 1 required");
    if (n == 1) return PermutationGroup({Permutation::identity(1)});
    // Generated by the transposition (0 1) and the n-cycle (0 1 … n−1).
    std::vector<int> transp(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) transp[static_cast<std::size_t>(i)] = i;
    std::swap(transp[0], transp[1]);
    std::vector<int> cyc(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) cyc[static_cast<std::size_t>(i)] = (i + 1) % n;
    return PermutationGroup({Permutation(std::move(transp)), Permutation(std::move(cyc))});
}

PermutationGroup cyclic_group(int n) {
    if (n < 1) throw std::invalid_argument("cyclic_group: n >= 1 required");
    std::vector<int> cyc(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) cyc[static_cast<std::size_t>(i)] = (i + 1) % n;
    return PermutationGroup({Permutation(std::move(cyc))});
}

PermutationGroup dihedral_group(int n) {
    if (n < 1) throw std::invalid_argument("dihedral_group: n >= 1 required");
    // Rotation r = (0 1 … n−1) and reflection s: i ↦ (n−i) mod n.
    std::vector<int> rot(static_cast<std::size_t>(n));
    std::vector<int> refl(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        rot[static_cast<std::size_t>(i)] = (i + 1) % n;
        refl[static_cast<std::size_t>(i)] = (n - i) % n;
    }
    return PermutationGroup({Permutation(std::move(rot)), Permutation(std::move(refl))});
}

PermutationGroup alternating_group(int n) {
    if (n < 1) throw std::invalid_argument("alternating_group: n >= 1 required");
    if (n < 3) return PermutationGroup({Permutation::identity(std::max(n, 1))});
    // Aₙ is generated by the 3-cycles (0 1 2), (0 1 3), …, (0 1 n−1).
    std::vector<Permutation> gens;
    for (int k = 2; k < n; ++k) {
        std::vector<int> img(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) img[static_cast<std::size_t>(i)] = i;
        // 3-cycle (0 1 k): 0→1, 1→k, k→0.
        img[0] = 1;
        img[1] = k;
        img[static_cast<std::size_t>(k)] = 0;
        gens.push_back(Permutation(std::move(img)));
    }
    return PermutationGroup(std::move(gens));
}

// ----- Integer partitions ----------------------------------------------------

Expr partition_count(int n) {
    if (n < 0) throw std::invalid_argument("partition_count: n >= 0 required");
    // Euler's pentagonal-number recurrence:
    //   p(n) = Σ_k (−1)^{k−1} [ p(n − g_k) + p(n − g'_k) ],  g_k = k(3k−1)/2.
    std::vector<mpz_class> p(static_cast<std::size_t>(n) + 1, mpz_class(0));
    p[0] = 1;
    for (int m = 1; m <= n; ++m) {
        mpz_class total = 0;
        for (int k = 1;; ++k) {
            int g1 = k * (3 * k - 1) / 2;
            int g2 = k * (3 * k + 1) / 2;
            if (g1 > m && g2 > m) break;
            int sign = (k % 2 == 1) ? 1 : -1;
            if (g1 <= m) total += sign * p[static_cast<std::size_t>(m - g1)];
            if (g2 <= m) total += sign * p[static_cast<std::size_t>(m - g2)];
        }
        p[static_cast<std::size_t>(m)] = total;
    }
    return make<Integer>(p[static_cast<std::size_t>(n)]);
}

std::vector<std::vector<int>> integer_partitions(int n) {
    if (n < 0) throw std::invalid_argument("integer_partitions: n >= 0 required");
    std::vector<std::vector<int>> out;
    std::vector<int> cur;
    // Generate non-increasing partitions: at each step pick the next part ≤ the
    // previous one (and ≤ remaining).
    auto rec = [&](auto&& self, int remaining, int max_part) -> void {
        if (remaining == 0) {
            out.push_back(cur);
            return;
        }
        for (int part = std::min(max_part, remaining); part >= 1; --part) {
            cur.push_back(part);
            self(self, remaining - part, part);
            cur.pop_back();
        }
    };
    rec(rec, n, n);
    return out;
}

}  // namespace sympp::combinatorics
