#pragma once

// Combinatorics & group theory — permutations, permutation groups, and integer
// partitions. Self-contained (GMP for big counts); mirrors the core of
// sympy.combinatorics.
//
//   Permutation p({1, 2, 0});          // array form: 0→1, 1→2, 2→0  (a 3-cycle)
//   p.order();                          // 3
//   p.sign();                           // +1 (even)
//   auto S3 = symmetric_group(3);
//   S3.order();                         // 6
//   S3.is_abelian();                    // false
//
// Reference: sympy/combinatorics/{permutations,perm_groups}.py, sympy.ntheory
// partition function.

#include <cstddef>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp::combinatorics {

// A permutation of {0, …, n−1} in array form: array_form()[i] is the image of i.
// Composition is function composition: (a∘b)(i) = a(b(i)).
class SYMPP_EXPORT Permutation {
public:
    explicit Permutation(std::vector<int> image);
    [[nodiscard]] static Permutation identity(int n);

    [[nodiscard]] int size() const noexcept { return static_cast<int>(image_.size()); }
    [[nodiscard]] int apply(int i) const { return image_[static_cast<std::size_t>(i)]; }
    [[nodiscard]] const std::vector<int>& array_form() const noexcept { return image_; }

    // a.compose(b) is a∘b. inverse() is the group inverse.
    [[nodiscard]] Permutation compose(const Permutation& other) const;
    [[nodiscard]] Permutation inverse() const;

    [[nodiscard]] int sign() const;                    // +1 if even, −1 if odd
    [[nodiscard]] long order() const;                  // least k>0 with pᵏ = id
    [[nodiscard]] bool is_identity() const;
    // Disjoint-cycle decomposition, fixed points omitted (SymPy's cyclic_form).
    [[nodiscard]] std::vector<std::vector<int>> cyclic_form() const;

    [[nodiscard]] bool operator==(const Permutation& o) const { return image_ == o.image_; }
    [[nodiscard]] bool operator!=(const Permutation& o) const { return image_ != o.image_; }

private:
    std::vector<int> image_;
};

// A finite permutation group given by generators. Element enumeration uses
// breadth-first closure (suitable for the small/medium groups used in teaching
// and tests, not industrial-scale groups).
class SYMPP_EXPORT PermutationGroup {
public:
    explicit PermutationGroup(std::vector<Permutation> generators);

    [[nodiscard]] int degree() const noexcept { return degree_; }
    [[nodiscard]] const std::vector<Permutation>& generators() const noexcept { return gens_; }

    [[nodiscard]] std::vector<Permutation> elements() const;  // full closure
    [[nodiscard]] long order() const;                         // |G|
    [[nodiscard]] bool contains(const Permutation& p) const;
    [[nodiscard]] bool is_abelian() const;
    // Transitive: the action on {0, …, degree−1} has a single orbit.
    [[nodiscard]] bool is_transitive() const;
    // Whether this group is a subgroup of `other` (every generator lies in it).
    [[nodiscard]] bool is_subgroup(const PermutationGroup& other) const;

private:
    std::vector<Permutation> gens_;
    int degree_;
};

// ----- Group actions: orbits & Pólya/Burnside enumeration --------------------
// Orbit of `point` under the group, ascending.
[[nodiscard]] SYMPP_EXPORT std::vector<int> orbit(const PermutationGroup& g, int point);
// Orbits partitioning {0, …, degree−1} (each ascending; by least element).
[[nodiscard]] SYMPP_EXPORT std::vector<std::vector<int>> orbits(const PermutationGroup& g);
// Number of distinct colorings of the `degree` points with `k` colors up to the
// group action — Burnside's lemma / Pólya: (1/|G|)·Σ_{g∈G} k^{cycles(g)}, where
// cycles(g) counts every disjoint cycle including fixed points.
[[nodiscard]] SYMPP_EXPORT Expr colorings_count(const PermutationGroup& g, int k);

// ----- Standard groups -------------------------------------------------------
[[nodiscard]] SYMPP_EXPORT PermutationGroup symmetric_group(int n);    // Sₙ, order n!
[[nodiscard]] SYMPP_EXPORT PermutationGroup cyclic_group(int n);       // Cₙ, order n
[[nodiscard]] SYMPP_EXPORT PermutationGroup dihedral_group(int n);     // Dₙ, order 2n
[[nodiscard]] SYMPP_EXPORT PermutationGroup alternating_group(int n);  // Aₙ, order n!/2

// ----- Integer partitions ----------------------------------------------------
// Partition function p(n): number of ways to write n as a sum of positive
// integers (order-insensitive). Arbitrary precision.
[[nodiscard]] SYMPP_EXPORT Expr partition_count(int n);
// All partitions of n, each as a non-increasing vector of parts.
[[nodiscard]] SYMPP_EXPORT std::vector<std::vector<int>> integer_partitions(int n);

}  // namespace sympp::combinatorics
