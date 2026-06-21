#include <sympp/vector/vector_calculus.hpp>

#include <stdexcept>

#include <sympp/calculus/diff.hpp>
#include <sympp/core/operators.hpp>
#include <sympp/core/rational.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp::vector {

VectorField gradient(const Expr& f, const std::vector<Expr>& coords) {
    VectorField g;
    g.reserve(coords.size());
    for (const auto& c : coords) g.push_back(diff(f, c));
    return g;
}

Expr divergence(const VectorField& field, const std::vector<Expr>& coords) {
    if (field.size() != coords.size()) {
        throw std::invalid_argument("divergence: field and coords sizes differ");
    }
    Expr s = S::Zero();
    for (std::size_t i = 0; i < coords.size(); ++i) s = s + diff(field[i], coords[i]);
    return simplify(s);
}

VectorField curl(const VectorField& f, const std::vector<Expr>& x) {
    if (f.size() != 3 || x.size() != 3) {
        throw std::invalid_argument("curl: requires three components and coordinates");
    }
    return {simplify(diff(f[2], x[1]) - diff(f[1], x[2])),
            simplify(diff(f[0], x[2]) - diff(f[2], x[0])),
            simplify(diff(f[1], x[0]) - diff(f[0], x[1]))};
}

Expr laplacian(const Expr& f, const std::vector<Expr>& coords) {
    Expr s = S::Zero();
    for (const auto& c : coords) s = s + diff(f, c, 2);
    return simplify(s);
}

ChristoffelSymbols christoffel(const Matrix& g, const std::vector<Expr>& x) {
    std::size_t n = x.size();
    Matrix gi = g.inverse();
    ChristoffelSymbols gamma;
    gamma.reserve(n);
    for (std::size_t k = 0; k < n; ++k) {
        Matrix gk = Matrix::zeros(n, n);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                Expr s = S::Zero();
                for (std::size_t l = 0; l < n; ++l) {
                    // ∂ᵢ g_lj + ∂ⱼ g_li − ∂_l g_ij
                    Expr term = diff(g.at(l, j), x[i]) + diff(g.at(l, i), x[j]) -
                                diff(g.at(i, j), x[l]);
                    s = s + gi.at(k, l) * term;
                }
                gk.set(i, j, simplify(mul(rational(1, 2), s)));
            }
        }
        gamma.push_back(gk);
    }
    return gamma;
}

Matrix ricci_tensor(const Matrix& g, const std::vector<Expr>& x) {
    std::size_t n = x.size();
    ChristoffelSymbols G = christoffel(g, x);
    Matrix R = Matrix::zeros(n, n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            Expr r = S::Zero();
            for (std::size_t k = 0; k < n; ++k) {
                r = r + diff(G[k].at(i, j), x[k]) - diff(G[k].at(i, k), x[j]);
                for (std::size_t l = 0; l < n; ++l) {
                    r = r + G[k].at(k, l) * G[l].at(i, j) -
                            G[k].at(j, l) * G[l].at(i, k);
                }
            }
            R.set(i, j, simplify(r));
        }
    }
    return R;
}

Expr ricci_scalar(const Matrix& g, const std::vector<Expr>& x) {
    std::size_t n = x.size();
    Matrix gi = g.inverse();
    Matrix R = ricci_tensor(g, x);
    Expr s = S::Zero();
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) s = s + gi.at(i, j) * R.at(i, j);
    }
    return simplify(s);
}

}  // namespace sympp::vector
