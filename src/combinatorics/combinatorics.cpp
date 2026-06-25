#include <sympp/combinatorics/combinatorics.hpp>

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>

#include <gmpxx.h>

#include <sympp/core/basic.hpp>
#include <sympp/core/integer.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/pow.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/symbol.hpp>

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

namespace {

// One level of a Schreier–Sims stabilizer chain: a base point and the
// transversal (coset representatives) of its orbit under the level's generators.
struct ChainLevel {
    int base;
    std::map<int, Permutation> transversal;  // x ↦ rep r with r(base) = x
};

// Orbit of `base` under `gens` with a Schreier transversal (left-action,
// r(base) = x), built breadth-first.
[[nodiscard]] std::map<int, Permutation> orbit_transversal(int base,
                                                           const std::vector<Permutation>& gens,
                                                           int deg) {
    std::map<int, Permutation> t;
    t.emplace(base, Permutation::identity(deg));
    std::vector<int> queue{base};
    for (std::size_t qi = 0; qi < queue.size(); ++qi) {
        int x = queue[qi];
        for (const auto& g : gens) {
            int y = g.apply(x);
            if (t.find(y) == t.end()) {
                t.emplace(y, g.compose(t.at(x)));  // g ∘ rep_x sends base → y
                queue.push_back(y);
            }
        }
    }
    return t;
}

// Schreier generators of the stabilizer of `base`: T[xˢ]⁻¹ ∘ s ∘ T[x] for every
// orbit point x and generator s (each fixes `base`); deduplicated, identities
// dropped. By Schreier's lemma these generate the full point stabilizer.
[[nodiscard]] std::vector<Permutation> schreier_generators(
    int base, const std::vector<Permutation>& gens, const std::map<int, Permutation>& t) {
    std::set<std::vector<int>> seen;
    std::vector<Permutation> out;
    for (const auto& [x, ux] : t) {
        for (const auto& s : gens) {
            Permutation g = s.compose(ux);
            Permutation sg = t.at(g.apply(base)).inverse().compose(g);
            if (sg.is_identity()) continue;
            if (seen.insert(sg.array_form()).second) out.push_back(sg);
        }
    }
    return out;
}

// Build the full stabilizer chain (Schreier–Sims). Each level stabilizes all
// previous base points; the product of the orbit sizes is |G|.
[[nodiscard]] std::vector<ChainLevel> build_chain(std::vector<Permutation> gens, int deg) {
    std::vector<ChainLevel> chain;
    while (true) {
        int base = -1;  // least point moved by some current generator
        for (const auto& g : gens) {
            for (int i = 0; i < deg; ++i) {
                if (g.apply(i) != i) {
                    base = (base < 0) ? i : std::min(base, i);
                    break;
                }
            }
        }
        if (base < 0) break;  // trivial — chain complete
        auto t = orbit_transversal(base, gens, deg);
        gens = schreier_generators(base, gens, t);
        chain.push_back({base, std::move(t)});
    }
    return chain;
}

}  // namespace

long PermutationGroup::order() const {
    long ord = 1;  // |G| = ∏ orbit sizes along the stabilizer chain
    for (const auto& level : build_chain(gens_, degree_)) {
        ord *= static_cast<long>(level.transversal.size());
    }
    return ord;
}

bool PermutationGroup::contains(const Permutation& p) const {
    if (p.size() != degree_) return false;
    // Sift p through the chain: strip a coset rep at each level; in the group iff
    // the residue is the identity.
    Permutation residue = p;
    for (const auto& level : build_chain(gens_, degree_)) {
        int x = residue.apply(level.base);
        auto it = level.transversal.find(x);
        if (it == level.transversal.end()) return false;  // base point not in orbit
        residue = it->second.inverse().compose(residue);   // T[x]⁻¹ ∘ residue fixes base
    }
    return residue.is_identity();
}

bool PermutationGroup::is_transitive() const {
    return static_cast<int>(orbit(*this, 0).size()) == degree_;
}

bool PermutationGroup::is_subgroup(const PermutationGroup& other) const {
    if (degree_ != other.degree_) return false;
    for (const auto& g : gens_)
        if (!other.contains(g)) return false;
    return true;
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

// ----- Group actions: orbits & Pólya/Burnside enumeration --------------------

std::vector<int> orbit(const PermutationGroup& g, int point) {
    std::set<int> seen{point};
    std::vector<int> frontier{point};
    while (!frontier.empty()) {
        std::vector<int> next;
        for (int p : frontier) {
            for (const auto& gen : g.generators()) {
                int q = gen.apply(p);
                if (seen.insert(q).second) next.push_back(q);
            }
        }
        frontier = std::move(next);
    }
    return std::vector<int>(seen.begin(), seen.end());  // std::set is ordered
}

std::vector<std::vector<int>> orbits(const PermutationGroup& g) {
    std::vector<std::vector<int>> result;
    std::vector<bool> placed(static_cast<std::size_t>(g.degree()), false);
    for (int p = 0; p < g.degree(); ++p) {
        if (placed[static_cast<std::size_t>(p)]) continue;
        auto o = orbit(g, p);
        for (int q : o) placed[static_cast<std::size_t>(q)] = true;
        result.push_back(std::move(o));
    }
    return result;
}

namespace {
// Number of disjoint cycles of `p`, counting fixed points as 1-cycles.
[[nodiscard]] int total_cycles(const Permutation& p) {
    int cycles = 0;
    std::vector<bool> visited(static_cast<std::size_t>(p.size()), false);
    for (int i = 0; i < p.size(); ++i) {
        if (visited[static_cast<std::size_t>(i)]) continue;
        ++cycles;
        int j = i;
        while (!visited[static_cast<std::size_t>(j)]) {
            visited[static_cast<std::size_t>(j)] = true;
            j = p.apply(j);
        }
    }
    return cycles;
}
}  // namespace

Expr colorings_count(const PermutationGroup& g, int k) {
    if (k < 0) throw std::invalid_argument("colorings_count: k must be non-negative");
    // Σ_{g∈G} k^{cycles(g)}, then divide by |G| (Burnside guarantees divisibility).
    mpz_class total = 0;
    auto elems = g.elements();
    for (const auto& e : elems) {
        mpz_class term;
        mpz_pow_ui(term.get_mpz_t(), mpz_class(k).get_mpz_t(),
                   static_cast<unsigned long>(total_cycles(e)));
        total += term;
    }
    return make<Integer>(total / static_cast<long>(elems.size()));
}

namespace {
// Cycle type of `p`: cycle length → number of cycles of that length (fixed points
// counted as length-1 cycles).
[[nodiscard]] std::map<int, int> cycle_type(const Permutation& p) {
    std::map<int, int> type;
    std::vector<bool> visited(static_cast<std::size_t>(p.size()), false);
    for (int i = 0; i < p.size(); ++i) {
        if (visited[static_cast<std::size_t>(i)]) continue;
        int len = 0;
        int j = i;
        while (!visited[static_cast<std::size_t>(j)]) {
            visited[static_cast<std::size_t>(j)] = true;
            ++len;
            j = p.apply(j);
        }
        ++type[len];
    }
    return type;
}
}  // namespace

Expr cycle_index(const PermutationGroup& g) {
    // Z(G) = (1/|G|)·Σ_{g∈G} ∏_len a_len^{cycles_len(g)}.
    auto elems = g.elements();
    Expr sum = make<Integer>(mpz_class(0));
    for (const auto& e : elems) {
        Expr monomial = make<Integer>(mpz_class(1));
        for (const auto& [len, count] : cycle_type(e)) {
            Expr a = symbol("a" + std::to_string(len));
            monomial = monomial * pow(a, make<Integer>(mpz_class(count)));
        }
        sum = sum + monomial;
    }
    return rational(1, static_cast<long>(elems.size())) * sum;
}

Expr necklaces(int n, int k) {
    if (n < 1) throw std::invalid_argument("necklaces: n >= 1 required");
    if (k < 0) throw std::invalid_argument("necklaces: k must be non-negative");
    // (1/n)·Σ_{d|n} φ(d)·k^{n/d}, divisibility guaranteed; computed exactly.
    auto totient = [](long m) {
        long result = m;
        long t = m;
        for (long p = 2; p * p <= t; ++p) {
            if (t % p == 0) {
                while (t % p == 0) t /= p;
                result -= result / p;
            }
        }
        if (t > 1) result -= result / t;
        return result;
    };
    mpz_class total = 0;
    for (long d = 1; d <= n; ++d) {
        if (n % d != 0) continue;
        mpz_class term;
        mpz_pow_ui(term.get_mpz_t(), mpz_class(k).get_mpz_t(),
                   static_cast<unsigned long>(n / d));
        total += totient(d) * term;
    }
    return make<Integer>(total / static_cast<long>(n));
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

// ----- Sylow theory ----------------------------------------------------------

namespace {
// Largest power p^a dividing m (a ≥ 0); returns p^a, or 1 when p ∤ m.
[[nodiscard]] long largest_prime_power_factor(long m, int p) {
    long power = 1;
    while (m % p == 0) {
        m /= p;
        power *= p;
    }
    return power;
}

// Whether n > 0 is a power of p (1 = p^0 counts).
[[nodiscard]] bool is_power_of(long n, int p) {
    if (n < 1) return false;
    while (n % p == 0) n /= p;
    return n == 1;
}
}  // namespace

long sylow_order(const PermutationGroup& g, int p) {
    if (p < 2) throw std::invalid_argument("sylow_order: p must be a prime >= 2");
    return largest_prime_power_factor(g.order(), p);
}

PermutationGroup sylow_subgroup(const PermutationGroup& g, int p) {
    if (p < 2) throw std::invalid_argument("sylow_subgroup: p must be a prime >= 2");
    const long target = sylow_order(g, p);

    // p^0 = 1: the trivial subgroup (generated by the identity).
    if (target == 1) return PermutationGroup({Permutation::identity(g.degree())});

    // Collect the p-elements of G: non-identity elements whose order is a power of
    // p. A Sylow p-subgroup consists entirely of such elements.
    std::vector<Permutation> p_elements;
    for (const auto& e : g.elements()) {
        if (!e.is_identity() && is_power_of(e.order(), p)) p_elements.push_back(e);
    }

    // Order of the subgroup generated by `gens` (via BFS closure over the seeds).
    auto closure_order = [&](const std::vector<Permutation>& gens) -> long {
        return PermutationGroup(gens).order();
    };

    // Greedily grow a p-subgroup: keep adding a p-element that strictly enlarges
    // the closure while staying a p-group, until the order reaches p^a.
    std::vector<Permutation> gens;
    long current = 1;
    bool progress = true;
    while (current < target && progress) {
        progress = false;
        for (const auto& cand : p_elements) {
            std::vector<Permutation> trial = gens;
            trial.push_back(cand);
            long ord = closure_order(trial);
            if (ord > current && is_power_of(ord, p) && ord <= target) {
                gens = std::move(trial);
                current = ord;
                progress = true;
                if (current == target) break;
            }
        }
    }

    if (gens.empty()) {
        throw std::logic_error("sylow_subgroup: no p-elements available to build subgroup");
    }

    PermutationGroup result(gens);
    // Guarantees: the result is a subgroup of G of order exactly p^a.
    if (result.order() != target) {
        throw std::logic_error("sylow_subgroup: failed to reach Sylow order p^a");
    }
    if (!result.is_subgroup(g)) {
        throw std::logic_error("sylow_subgroup: constructed group is not a subgroup of G");
    }
    return result;
}

// ----- Conjugacy, center & derived subgroup ----------------------------------

std::vector<std::vector<Permutation>> conjugacy_classes(const PermutationGroup& g) {
    auto elems = g.elements();
    std::set<std::vector<int>> classified;
    std::vector<std::vector<Permutation>> classes;
    for (const auto& a : elems) {
        if (classified.count(a.array_form())) continue;
        // Build the orbit of `a` under conjugation by every group element.
        std::set<std::vector<int>> members;
        std::vector<Permutation> cls;
        for (const auto& x : elems) {
            Permutation conj = x.compose(a).compose(x.inverse());  // x·a·x⁻¹
            if (members.insert(conj.array_form()).second) {
                cls.push_back(conj);
                classified.insert(conj.array_form());
            }
        }
        classes.push_back(std::move(cls));
    }
    return classes;
}

long num_conjugacy_classes(const PermutationGroup& g) {
    return static_cast<long>(conjugacy_classes(g).size());
}

PermutationGroup group_center(const PermutationGroup& g) {
    // Z(G) = elements commuting with every generator (hence with all of G).
    std::vector<Permutation> central;
    for (const auto& e : g.elements()) {
        bool commutes = true;
        for (const auto& gen : g.generators()) {
            if (e.compose(gen) != gen.compose(e)) {
                commutes = false;
                break;
            }
        }
        if (commutes) central.push_back(e);
    }
    if (central.empty()) central.push_back(Permutation::identity(g.degree()));
    return PermutationGroup(std::move(central));
}

PermutationGroup derived_subgroup(const PermutationGroup& g) {
    // G' is generated by all commutators a·b·a⁻¹·b⁻¹ of group elements.
    auto elems = g.elements();
    std::vector<Permutation> commutators;
    std::set<std::vector<int>> seen;
    for (const auto& a : elems) {
        for (const auto& b : elems) {
            Permutation c = a.compose(b).compose(a.inverse()).compose(b.inverse());
            if (!c.is_identity() && seen.insert(c.array_form()).second) {
                commutators.push_back(c);
            }
        }
    }
    if (commutators.empty()) {
        // Abelian group: G' is trivial.
        return PermutationGroup({Permutation::identity(g.degree())});
    }
    return PermutationGroup(std::move(commutators));  // BFS closure on use
}

bool is_solvable(const PermutationGroup& g) {
    // Iterate the derived series G ▷ G' ▷ G'' ▷ … Each step's order is
    // non-increasing; stop when it stabilizes. Solvable iff it reaches 1.
    long prev = g.order();
    PermutationGroup cur = derived_subgroup(g);
    while (true) {
        long ord = cur.order();
        if (ord == 1) return true;       // reached the trivial group
        if (ord == prev) return false;   // series stabilized above 1 (e.g. A₅)
        prev = ord;
        cur = derived_subgroup(cur);
    }
}

namespace {

// One step of the lower central series: [G, H] — the subgroup generated by all
// commutators [a, b] = a·b·a⁻¹·b⁻¹ with a ∈ G and b ∈ H.
[[nodiscard]] PermutationGroup lower_central_step(const PermutationGroup& g,
                                                  const PermutationGroup& h) {
    auto gelems = g.elements();
    auto helems = h.elements();
    std::vector<Permutation> commutators;
    std::set<std::vector<int>> seen;
    for (const auto& a : gelems) {
        for (const auto& b : helems) {
            Permutation c = a.compose(b).compose(a.inverse()).compose(b.inverse());
            if (!c.is_identity() && seen.insert(c.array_form()).second) {
                commutators.push_back(c);
            }
        }
    }
    if (commutators.empty()) {
        return PermutationGroup({Permutation::identity(g.degree())});
    }
    return PermutationGroup(std::move(commutators));  // BFS closure on use
}

}  // namespace

bool is_nilpotent(const PermutationGroup& g) {
    // Lower central series γ₁ = G, γ_{i+1} = [G, γ_i]. Each order is
    // non-increasing; nilpotent iff it reaches the trivial group.
    long prev = g.order();
    PermutationGroup cur = lower_central_step(g, g);  // γ₂ = [G, G]
    while (true) {
        long ord = cur.order();
        if (ord == 1) return true;       // reached the trivial group
        if (ord == prev) return false;   // series stabilized above 1
        prev = ord;
        cur = lower_central_step(g, cur);
    }
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
