// Engineering Algebra — 04: Complex numbers.
//
// Phasors, impedances and stability all live in the complex plane. Covers
// complex arithmetic, real/imaginary parts and conjugates, Euler's formula,
// and the roots of unity (DFT twiddle factors / equally-spaced poles).

#include <sympp/sympp.hpp>

#include <iostream>
#include <sympp/core/rewrite.hpp>

int main() {
    using namespace sympp;
    auto I = S::I();
    auto x = symbol("x");

    std::cout << "=== Complex arithmetic ===\n";
    auto z1 = integer(2) + integer(3) * I;
    auto z2 = integer(1) - integer(2) * I;
    std::cout << "z1 = 2+3i, z2 = 1-2i\n";
    std::cout << "z1 * z2     = " << simplify(mul(z1, z2))->str() << "\n";
    std::cout << "z1 / z2     = " << simplify(z1 / z2)->str() << "\n";

    std::cout << "\n=== Parts, conjugate, magnitude ===\n";
    std::cout << "Re(z1)      = " << simplify(re(z1))->str() << "\n";
    std::cout << "Im(z1)      = " << simplify(im(z1))->str() << "\n";
    std::cout << "conj(z1)    = " << simplify(conjugate(z1))->str() << "\n";
    std::cout << "|z1|        = " << simplify(abs(z1))->str()
              << "  (= sqrt(13))\n";
    // |z|^2 = z * conj(z).
    std::cout << "z1*conj(z1) = " << simplify(mul(z1, conjugate(z1)))->str()
              << "\n";

    std::cout << "\n=== Euler's formula ===\n";
    // e^{i pi} + 1 = 0.
    std::cout << "e^(i*pi) + 1 = "
              << simplify(exp(mul(I, S::Pi())) + integer(1))->str()
              << "   (Euler's identity)\n";
    // cos(x) as the real part of e^{ix}.
    std::cout << "cos(x) via exp: "
              << rewrite(cos(x), RewriteTarget::Exp)->str() << "\n";

    std::cout << "\n=== Roots of unity (DFT / pole placement) ===\n";
    // The three cube roots of unity.
    std::cout << "x^3 = 1  ->";
    for (const auto& r : solve(pow(x, integer(3)) - integer(1), x))
        std::cout << "  " << simplify(r)->str();
    std::cout << "\n";
    // Fourth roots: 1, i, -1, -i.
    std::cout << "x^4 = 1  ->";
    for (const auto& r : solve(pow(x, integer(4)) - integer(1), x))
        std::cout << "  " << simplify(r)->str();
    std::cout << "\n";
    return 0;
}
