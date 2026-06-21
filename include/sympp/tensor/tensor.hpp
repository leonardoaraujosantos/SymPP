#pragma once

// Dense N-dimensional tensors with symbolic components — tensor (outer)
// product, index contraction, addition/scaling, and metric-based index
// raising/lowering.
//
//   Tensor m({2, 2}, {integer(1), integer(2), integer(3), integer(4)});
//   contract(m, 0, 1);                 // → trace 5 (rank-0 tensor)
//   tensor_product(u, v);              // → outer product
//
// Reference: sympy/tensor/array/* (Array, tensorproduct, tensorcontraction).

#include <cstddef>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp::tensor {

class SYMPP_EXPORT Tensor {
public:
    Tensor(std::vector<std::size_t> shape, std::vector<Expr> data);

    [[nodiscard]] const std::vector<std::size_t>& shape() const noexcept { return shape_; }
    [[nodiscard]] std::size_t rank() const noexcept { return shape_.size(); }
    [[nodiscard]] std::size_t size() const noexcept { return data_.size(); }
    [[nodiscard]] const std::vector<Expr>& data() const noexcept { return data_; }

    [[nodiscard]] const Expr& at(const std::vector<std::size_t>& idx) const;
    void set(const std::vector<std::size_t>& idx, Expr value);

    [[nodiscard]] bool operator==(const Tensor& o) const;

private:
    std::vector<std::size_t> shape_;
    std::vector<Expr> data_;  // row-major
};

// Outer product: shape is a.shape ++ b.shape.
[[nodiscard]] SYMPP_EXPORT Tensor tensor_product(const Tensor& a, const Tensor& b);
// Contract axes i and j (which must have equal dimension) by summing the
// diagonal; the result drops both axes.
[[nodiscard]] SYMPP_EXPORT Tensor contract(const Tensor& t, std::size_t i, std::size_t j);
[[nodiscard]] SYMPP_EXPORT Tensor add(const Tensor& a, const Tensor& b);
[[nodiscard]] SYMPP_EXPORT Tensor scalar_mul(const Tensor& t, const Expr& s);

// Lower / raise `axis` of `t` with a rank-2 metric (or its inverse):
//   T'…ₐ… = Σ_b g_{ab} T…ᵇ…
[[nodiscard]] SYMPP_EXPORT Tensor lower_index(const Tensor& t, std::size_t axis,
                                              const Tensor& metric);
[[nodiscard]] SYMPP_EXPORT Tensor raise_index(const Tensor& t, std::size_t axis,
                                              const Tensor& inverse_metric);

}  // namespace sympp::tensor
