#include <sympp/tensor/tensor.hpp>

#include <numeric>
#include <stdexcept>

#include <sympp/core/operators.hpp>
#include <sympp/core/singletons.hpp>
#include <sympp/simplify/simplify.hpp>

namespace sympp::tensor {

namespace {
[[nodiscard]] std::size_t product(const std::vector<std::size_t>& shape) {
    std::size_t n = 1;
    for (auto d : shape) n *= d;
    return n;
}
// Row-major flat offset for a multi-index.
[[nodiscard]] std::size_t flat(const std::vector<std::size_t>& shape,
                               const std::vector<std::size_t>& idx) {
    std::size_t off = 0;
    for (std::size_t k = 0; k < shape.size(); ++k) off = off * shape[k] + idx[k];
    return off;
}
}  // namespace

Tensor::Tensor(std::vector<std::size_t> shape, std::vector<Expr> data)
    : shape_(std::move(shape)), data_(std::move(data)) {
    if (data_.size() != product(shape_)) {
        throw std::invalid_argument("Tensor: data size does not match shape");
    }
}

const Expr& Tensor::at(const std::vector<std::size_t>& idx) const {
    return data_[flat(shape_, idx)];
}
void Tensor::set(const std::vector<std::size_t>& idx, Expr value) {
    data_[flat(shape_, idx)] = std::move(value);
}

bool Tensor::operator==(const Tensor& o) const {
    if (shape_ != o.shape_) return false;
    for (std::size_t i = 0; i < data_.size(); ++i) {
        if (!(simplify(data_[i] - o.data_[i]) == S::Zero())) return false;
    }
    return true;
}

Tensor tensor_product(const Tensor& a, const Tensor& b) {
    std::vector<std::size_t> shape = a.shape();
    shape.insert(shape.end(), b.shape().begin(), b.shape().end());
    std::vector<Expr> data;
    data.reserve(a.size() * b.size());
    for (const auto& x : a.data()) {
        for (const auto& y : b.data()) data.push_back(x * y);
    }
    return {std::move(shape), std::move(data)};
}

Tensor contract(const Tensor& t, std::size_t i, std::size_t j) {
    if (i == j || i >= t.rank() || j >= t.rank()) {
        throw std::invalid_argument("contract: invalid axes");
    }
    if (i > j) std::swap(i, j);
    if (t.shape()[i] != t.shape()[j]) {
        throw std::invalid_argument("contract: axes have different dimensions");
    }
    std::size_t dim = t.shape()[i];
    std::vector<std::size_t> out_shape;
    for (std::size_t k = 0; k < t.rank(); ++k) {
        if (k != i && k != j) out_shape.push_back(t.shape()[k]);
    }
    std::size_t out_size = product(out_shape);
    std::vector<Expr> out(out_size, S::Zero());

    // Iterate over every result multi-index, summing the diagonal of i,j.
    std::vector<std::size_t> oidx(out_shape.size(), 0);
    for (std::size_t flat_out = 0; flat_out < out_size; ++flat_out) {
        Expr s = S::Zero();
        for (std::size_t d = 0; d < dim; ++d) {
            std::vector<std::size_t> full(t.rank());
            std::size_t p = 0;
            for (std::size_t k = 0; k < t.rank(); ++k) {
                if (k == i || k == j) full[k] = d;
                else full[k] = oidx[p++];
            }
            s = s + t.at(full);
        }
        out[flat_out] = simplify(s);
        // increment oidx (row-major)
        for (std::size_t k = out_shape.size(); k-- > 0;) {
            if (++oidx[k] < out_shape[k]) break;
            oidx[k] = 0;
        }
    }
    if (out_shape.empty()) return {{}, {out[0]}};  // rank-0 scalar
    return {std::move(out_shape), std::move(out)};
}

Tensor add(const Tensor& a, const Tensor& b) {
    if (a.shape() != b.shape()) throw std::invalid_argument("add: shape mismatch");
    std::vector<Expr> data;
    data.reserve(a.size());
    for (std::size_t k = 0; k < a.size(); ++k) data.push_back(simplify(a.data()[k] + b.data()[k]));
    return {a.shape(), std::move(data)};
}

Tensor scalar_mul(const Tensor& t, const Expr& s) {
    std::vector<Expr> data;
    data.reserve(t.size());
    for (const auto& x : t.data()) data.push_back(simplify(mul(s, x)));
    return {t.shape(), std::move(data)};
}

namespace {
// Generic index transform: T'…ₐ… = Σ_b M_{ab} T…ᵇ…  (M a rank-2 tensor).
[[nodiscard]] Tensor transform_index(const Tensor& t, std::size_t axis, const Tensor& m) {
    if (m.rank() != 2 || m.shape()[0] != m.shape()[1]) {
        throw std::invalid_argument("index transform: metric must be square rank-2");
    }
    if (axis >= t.rank() || t.shape()[axis] != m.shape()[0]) {
        throw std::invalid_argument("index transform: axis/metric dimension mismatch");
    }
    std::size_t dim = m.shape()[0];
    std::vector<Expr> out(t.size(), S::Zero());
    std::vector<std::size_t> idx(t.rank(), 0);
    for (std::size_t flat_i = 0; flat_i < t.size(); ++flat_i) {
        Expr s = S::Zero();
        std::vector<std::size_t> src = idx;
        for (std::size_t b = 0; b < dim; ++b) {
            src[axis] = b;
            s = s + m.at({idx[axis], b}) * t.at(src);
        }
        out[flat_i] = simplify(s);
        for (std::size_t k = t.rank(); k-- > 0;) {
            if (++idx[k] < t.shape()[k]) break;
            idx[k] = 0;
        }
    }
    return {t.shape(), std::move(out)};
}
}  // namespace

Tensor lower_index(const Tensor& t, std::size_t axis, const Tensor& metric) {
    return transform_index(t, axis, metric);
}
Tensor raise_index(const Tensor& t, std::size_t axis, const Tensor& inverse_metric) {
    return transform_index(t, axis, inverse_metric);
}

}  // namespace sympp::tensor
