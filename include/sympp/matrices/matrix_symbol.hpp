#pragma once

// MatrixSymbol — opaque symbolic matrix of a known shape, plus a small
// expression tree (MatMul, MatAdd, MatPow, MTranspose, MTrace,
// MDeterminant, MInverse, BlockMatrix) that records operations without
// materializing the underlying entries.
//
// Modeled on SymPy's matrices/expressions/* layer but pared down: each
// node carries the shape it produces, plus its operands, and prints
// in the same form as SymPy's matrix expressions for cross-validation.
//
// The nodes are *not* sympp::Basic-derived because they operate on
// matrix shapes rather than scalar expressions; they live alongside
// the dense Matrix class. Callers can mix the two (e.g. by
// substituting a concrete Matrix for a MatrixSymbol via .as_dense()).
//
// Reference: sympy/matrices/expressions/{matexpr,matmul,matadd,
// matpow,transpose,trace,determinant,inverse,blockmatrix}.py

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <sympp/core/api.hpp>
#include <sympp/fwd.hpp>

namespace sympp {

class MatrixExpr;
using MatrixExprPtr = std::shared_ptr<const MatrixExpr>;

enum class MatExprKind : std::uint8_t {
    Symbol,
    Add,
    Mul,
    Pow,
    Transpose,
    Trace,
    Determinant,
    Inverse,
    Block,
};

class SYMPP_EXPORT MatrixExpr {
public:
    virtual ~MatrixExpr() = default;
    [[nodiscard]] virtual MatExprKind kind() const noexcept = 0;
    [[nodiscard]] virtual std::size_t rows() const noexcept = 0;
    [[nodiscard]] virtual std::size_t cols() const noexcept = 0;
    [[nodiscard]] virtual std::string str() const = 0;
};

// Concrete: opaque symbol "A" with shape m×n.
class SYMPP_EXPORT MatrixSymbol final : public MatrixExpr {
public:
    MatrixSymbol(std::string name, std::size_t rows, std::size_t cols);
    [[nodiscard]] MatExprKind kind() const noexcept override {
        return MatExprKind::Symbol;
    }
    [[nodiscard]] std::size_t rows() const noexcept override { return rows_; }
    [[nodiscard]] std::size_t cols() const noexcept override { return cols_; }
    [[nodiscard]] std::string str() const override { return name_; }
    [[nodiscard]] const std::string& name() const noexcept { return name_; }

private:
    std::string name_;
    std::size_t rows_, cols_;
};

class SYMPP_EXPORT MatAdd final : public MatrixExpr {
public:
    MatAdd(MatrixExprPtr a, MatrixExprPtr b);
    [[nodiscard]] MatExprKind kind() const noexcept override {
        return MatExprKind::Add;
    }
    [[nodiscard]] std::size_t rows() const noexcept override { return a_->rows(); }
    [[nodiscard]] std::size_t cols() const noexcept override { return a_->cols(); }
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] const MatrixExprPtr& lhs() const { return a_; }
    [[nodiscard]] const MatrixExprPtr& rhs() const { return b_; }

private:
    MatrixExprPtr a_, b_;
};

class SYMPP_EXPORT MatMul final : public MatrixExpr {
public:
    MatMul(MatrixExprPtr a, MatrixExprPtr b);
    [[nodiscard]] MatExprKind kind() const noexcept override {
        return MatExprKind::Mul;
    }
    [[nodiscard]] std::size_t rows() const noexcept override { return a_->rows(); }
    [[nodiscard]] std::size_t cols() const noexcept override { return b_->cols(); }
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] const MatrixExprPtr& lhs() const { return a_; }
    [[nodiscard]] const MatrixExprPtr& rhs() const { return b_; }

private:
    MatrixExprPtr a_, b_;
};

class SYMPP_EXPORT MatPow final : public MatrixExpr {
public:
    MatPow(MatrixExprPtr a, int n);
    [[nodiscard]] MatExprKind kind() const noexcept override {
        return MatExprKind::Pow;
    }
    [[nodiscard]] std::size_t rows() const noexcept override { return a_->rows(); }
    [[nodiscard]] std::size_t cols() const noexcept override { return a_->cols(); }
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] int exponent() const noexcept { return n_; }
    [[nodiscard]] const MatrixExprPtr& base() const { return a_; }

private:
    MatrixExprPtr a_;
    int n_;
};

class SYMPP_EXPORT MTranspose final : public MatrixExpr {
public:
    explicit MTranspose(MatrixExprPtr a);
    [[nodiscard]] MatExprKind kind() const noexcept override {
        return MatExprKind::Transpose;
    }
    [[nodiscard]] std::size_t rows() const noexcept override { return a_->cols(); }
    [[nodiscard]] std::size_t cols() const noexcept override { return a_->rows(); }
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] const MatrixExprPtr& operand() const { return a_; }

private:
    MatrixExprPtr a_;
};

class SYMPP_EXPORT MInverse final : public MatrixExpr {
public:
    explicit MInverse(MatrixExprPtr a);
    [[nodiscard]] MatExprKind kind() const noexcept override {
        return MatExprKind::Inverse;
    }
    [[nodiscard]] std::size_t rows() const noexcept override { return a_->rows(); }
    [[nodiscard]] std::size_t cols() const noexcept override { return a_->cols(); }
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] const MatrixExprPtr& operand() const { return a_; }

private:
    MatrixExprPtr a_;
};

// Trace and Determinant return a *scalar* shape — we model them as
// MatrixExpr returning a 1×1 matrix and offer a convenience helper to
// extract the underlying scalar Expr when one becomes available.
class SYMPP_EXPORT MTrace final : public MatrixExpr {
public:
    explicit MTrace(MatrixExprPtr a);
    [[nodiscard]] MatExprKind kind() const noexcept override {
        return MatExprKind::Trace;
    }
    [[nodiscard]] std::size_t rows() const noexcept override { return 1; }
    [[nodiscard]] std::size_t cols() const noexcept override { return 1; }
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] const MatrixExprPtr& operand() const { return a_; }

private:
    MatrixExprPtr a_;
};

class SYMPP_EXPORT MDeterminant final : public MatrixExpr {
public:
    explicit MDeterminant(MatrixExprPtr a);
    [[nodiscard]] MatExprKind kind() const noexcept override {
        return MatExprKind::Determinant;
    }
    [[nodiscard]] std::size_t rows() const noexcept override { return 1; }
    [[nodiscard]] std::size_t cols() const noexcept override { return 1; }
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] const MatrixExprPtr& operand() const { return a_; }

private:
    MatrixExprPtr a_;
};

// BlockMatrix: 2D grid of MatrixExpr blocks with compatible shapes.
class SYMPP_EXPORT BlockMatrix final : public MatrixExpr {
public:
    using BlockGrid = std::vector<std::vector<MatrixExprPtr>>;
    explicit BlockMatrix(BlockGrid blocks);
    [[nodiscard]] MatExprKind kind() const noexcept override {
        return MatExprKind::Block;
    }
    [[nodiscard]] std::size_t rows() const noexcept override { return total_rows_; }
    [[nodiscard]] std::size_t cols() const noexcept override { return total_cols_; }
    [[nodiscard]] std::string str() const override;
    [[nodiscard]] const BlockGrid& blocks() const { return blocks_; }

private:
    BlockGrid blocks_;
    std::size_t total_rows_, total_cols_;
};

// ---- Factories -------------------------------------------------------------
[[nodiscard]] SYMPP_EXPORT MatrixExprPtr matrix_symbol(
    std::string_view name, std::size_t rows, std::size_t cols);

[[nodiscard]] SYMPP_EXPORT MatrixExprPtr matadd(MatrixExprPtr a, MatrixExprPtr b);
[[nodiscard]] SYMPP_EXPORT MatrixExprPtr matmul(MatrixExprPtr a, MatrixExprPtr b);
[[nodiscard]] SYMPP_EXPORT MatrixExprPtr matpow(MatrixExprPtr a, int n);
[[nodiscard]] SYMPP_EXPORT MatrixExprPtr mtranspose(MatrixExprPtr a);
[[nodiscard]] SYMPP_EXPORT MatrixExprPtr minverse(MatrixExprPtr a);
[[nodiscard]] SYMPP_EXPORT MatrixExprPtr mtrace(MatrixExprPtr a);
[[nodiscard]] SYMPP_EXPORT MatrixExprPtr mdeterminant(MatrixExprPtr a);
[[nodiscard]] SYMPP_EXPORT MatrixExprPtr block_matrix(
    BlockMatrix::BlockGrid blocks);

}  // namespace sympp
