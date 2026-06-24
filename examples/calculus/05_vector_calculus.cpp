// Engineering Calculus — 05: Vector calculus.
//
// The language of fields: electromagnetics, fluid flow, heat transfer.
// Gradient, divergence, curl and the Laplacian, plus the two fundamental
// vanishing identities (curl grad = 0, div curl = 0).

#include <sympp/sympp.hpp>
#include <sympp/vector/vector_calculus.hpp>

#include <iostream>

int main() {
    using namespace sympp;
    namespace vc = sympp::vector;

    auto x = symbol("x");
    auto y = symbol("y");
    auto z = symbol("z");
    std::vector<Expr> coords = {x, y, z};

    auto show = [](const std::string& name, const vc::VectorField& v) {
        std::cout << name << " = {" << v[0]->str() << ",  " << v[1]->str()
                  << ",  " << v[2]->str() << "}\n";
    };

    std::cout << "=== Gradient (steepest ascent of a scalar field) ===\n";
    auto phi = x * y * z;  // a scalar potential
    std::cout << "phi = x y z\n";
    show("grad phi", vc::gradient(phi, coords));

    std::cout << "\n=== Divergence (source density of a vector field) ===\n";
    vc::VectorField F = {x, y, z};  // radial field
    std::cout << "F = {x, y, z}\n";
    std::cout << "div F = " << vc::divergence(F, coords)->str()
              << "  (uniform source)\n";

    std::cout << "\n=== Curl (rotation of a vector field) ===\n";
    vc::VectorField G = {mul(integer(-1), y), x, integer(0)};  // rigid rotation
    std::cout << "G = {-y, x, 0}\n";
    show("curl G", vc::curl(G, coords));

    std::cout << "\n=== Laplacian (diffusion / steady heat) ===\n";
    auto u = pow(x, integer(2)) + pow(y, integer(2)) + pow(z, integer(2));
    std::cout << "u = x^2 + y^2 + z^2\n";
    std::cout << "laplacian u = " << vc::laplacian(u, coords)->str() << "\n";
    // A harmonic function has zero Laplacian (Laplace's equation).
    auto harm = pow(x, integer(2)) - pow(y, integer(2));
    std::cout << "laplacian (x^2 - y^2) = "
              << vc::laplacian(harm, coords)->str() << "  (harmonic)\n";

    std::cout << "\n=== Fundamental identities ===\n";
    // curl(grad phi) = 0 for any scalar field.
    auto cg = vc::curl(vc::gradient(phi, coords), coords);
    std::cout << "curl(grad phi) = {" << simplify(cg[0])->str() << ", "
              << simplify(cg[1])->str() << ", " << simplify(cg[2])->str()
              << "}\n";
    // div(curl G) = 0 for any vector field.
    std::cout << "div(curl G)    = "
              << simplify(vc::divergence(vc::curl(G, coords), coords))->str()
              << "\n";
    return 0;
}
