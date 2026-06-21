#pragma once

// Discrete transforms over sequences of expressions: the discrete Fourier
// transform (exact, via roots of unity), the number-theoretic transform, linear
// convolution, and the subset (Möbius) transform.
//
//   fft({1, 2, 3, 4});                 // → {10, -2 - 2i, -2, -2 + 2i}
//   convolution({1, 2, 3}, {4, 5, 6}); // → {4, 13, 28, 27, 18}
//
// fft/ifft/ntt/intt/mobius pad the input to the next power of two (as SymPy
// does). ntt/intt require a prime p with the transform length dividing p − 1.
//
// Reference: sympy/discrete/{transforms,convolutions}.py.

#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp::discrete {

[[nodiscard]] SYMPP_EXPORT std::vector<Expr> fft(const std::vector<Expr>& seq);
[[nodiscard]] SYMPP_EXPORT std::vector<Expr> ifft(const std::vector<Expr>& seq);
[[nodiscard]] SYMPP_EXPORT std::vector<Expr> ntt(const std::vector<Expr>& seq,
                                                 const Expr& prime);
[[nodiscard]] SYMPP_EXPORT std::vector<Expr> intt(const std::vector<Expr>& seq,
                                                  const Expr& prime);
[[nodiscard]] SYMPP_EXPORT std::vector<Expr> convolution(const std::vector<Expr>& a,
                                                         const std::vector<Expr>& b);
[[nodiscard]] SYMPP_EXPORT std::vector<Expr> mobius_transform(const std::vector<Expr>& seq);
[[nodiscard]] SYMPP_EXPORT std::vector<Expr> inverse_mobius_transform(
    const std::vector<Expr>& seq);

}  // namespace sympp::discrete
