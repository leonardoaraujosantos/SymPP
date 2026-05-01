#include <sympp/matrices/matrix_symbol.hpp>

#include <stdexcept>
#include <string>

namespace sympp {

namespace {

// Return a parenthesized form of an operand when its outer printing
// would associate ambiguously inside a binary operator.
[[nodiscard]] std::string maybe_paren(const MatrixExprPtr& e) {
    switch (e->kind()) {
        case MatExprKind::Add:
        case MatExprKind::Mul:
            return "(" + e->str() + ")";
        default:
            return e->str();
    }
}

}  // namespace

// ----- MatrixSymbol --------------------------------------------------------

MatrixSymbol::MatrixSymbol(std::string name, std::size_t rows, std::size_t cols)
    : name_(std::move(name)), rows_(rows), cols_(cols) {}

// ----- MatAdd --------------------------------------------------------------

MatAdd::MatAdd(MatrixExprPtr a, MatrixExprPtr b)
    : a_(std::move(a)), b_(std::move(b)) {
    if (a_->rows() != b_->rows() || a_->cols() != b_->cols()) {
        throw std::invalid_argument("MatAdd: shape mismatch");
    }
}

std::string MatAdd::str() const {
    return a_->str() + " + " + b_->str();
}

// ----- MatMul --------------------------------------------------------------

MatMul::MatMul(MatrixExprPtr a, MatrixExprPtr b)
    : a_(std::move(a)), b_(std::move(b)) {
    if (a_->cols() != b_->rows()) {
        throw std::invalid_argument("MatMul: inner dimension mismatch");
    }
}

std::string MatMul::str() const {
    return maybe_paren(a_) + "*" + maybe_paren(b_);
}

// ----- MatPow --------------------------------------------------------------

MatPow::MatPow(MatrixExprPtr a, int n)
    : a_(std::move(a)), n_(n) {
    if (a_->rows() != a_->cols()) {
        throw std::invalid_argument("MatPow: square operand required");
    }
}

std::string MatPow::str() const {
    return maybe_paren(a_) + "**" + std::to_string(n_);
}

// ----- MTranspose ----------------------------------------------------------

MTranspose::MTranspose(MatrixExprPtr a) : a_(std::move(a)) {}
std::string MTranspose::str() const { return maybe_paren(a_) + ".T"; }

// ----- MInverse ------------------------------------------------------------

MInverse::MInverse(MatrixExprPtr a) : a_(std::move(a)) {
    if (a_->rows() != a_->cols()) {
        throw std::invalid_argument("MInverse: square operand required");
    }
}
std::string MInverse::str() const { return maybe_paren(a_) + "**(-1)"; }

// ----- MTrace --------------------------------------------------------------

MTrace::MTrace(MatrixExprPtr a) : a_(std::move(a)) {
    if (a_->rows() != a_->cols()) {
        throw std::invalid_argument("MTrace: square operand required");
    }
}
std::string MTrace::str() const { return "Trace(" + a_->str() + ")"; }

// ----- MDeterminant --------------------------------------------------------

MDeterminant::MDeterminant(MatrixExprPtr a) : a_(std::move(a)) {
    if (a_->rows() != a_->cols()) {
        throw std::invalid_argument("MDeterminant: square operand required");
    }
}
std::string MDeterminant::str() const {
    return "Determinant(" + a_->str() + ")";
}

// ----- BlockMatrix ---------------------------------------------------------

BlockMatrix::BlockMatrix(BlockGrid blocks) : blocks_(std::move(blocks)) {
    if (blocks_.empty() || blocks_[0].empty()) {
        throw std::invalid_argument("BlockMatrix: empty grid");
    }
    // Validate row-shapes: each block in row r must share the same row count.
    // Validate col-shapes: each block in col c must share the same col count.
    const std::size_t R = blocks_.size();
    const std::size_t C = blocks_[0].size();
    for (std::size_t r = 0; r < R; ++r) {
        if (blocks_[r].size() != C) {
            throw std::invalid_argument("BlockMatrix: ragged grid");
        }
    }
    std::vector<std::size_t> row_h(R), col_w(C);
    for (std::size_t r = 0; r < R; ++r) row_h[r] = blocks_[r][0]->rows();
    for (std::size_t c = 0; c < C; ++c) col_w[c] = blocks_[0][c]->cols();
    for (std::size_t r = 0; r < R; ++r) {
        for (std::size_t c = 0; c < C; ++c) {
            if (blocks_[r][c]->rows() != row_h[r]
                || blocks_[r][c]->cols() != col_w[c]) {
                throw std::invalid_argument(
                    "BlockMatrix: block shape inconsistent with row/col header");
            }
        }
    }
    total_rows_ = 0;
    for (auto h : row_h) total_rows_ += h;
    total_cols_ = 0;
    for (auto w : col_w) total_cols_ += w;
}

std::string BlockMatrix::str() const {
    std::string s = "BlockMatrix([";
    for (std::size_t r = 0; r < blocks_.size(); ++r) {
        if (r > 0) s += ", ";
        s += "[";
        for (std::size_t c = 0; c < blocks_[r].size(); ++c) {
            if (c > 0) s += ", ";
            s += blocks_[r][c]->str();
        }
        s += "]";
    }
    s += "])";
    return s;
}

// ---- Factories ------------------------------------------------------------

MatrixExprPtr matrix_symbol(std::string_view name, std::size_t rows, std::size_t cols) {
    return std::make_shared<MatrixSymbol>(std::string(name), rows, cols);
}
MatrixExprPtr matadd(MatrixExprPtr a, MatrixExprPtr b) {
    return std::make_shared<MatAdd>(std::move(a), std::move(b));
}
MatrixExprPtr matmul(MatrixExprPtr a, MatrixExprPtr b) {
    return std::make_shared<MatMul>(std::move(a), std::move(b));
}
MatrixExprPtr matpow(MatrixExprPtr a, int n) {
    return std::make_shared<MatPow>(std::move(a), n);
}
MatrixExprPtr mtranspose(MatrixExprPtr a) {
    return std::make_shared<MTranspose>(std::move(a));
}
MatrixExprPtr minverse(MatrixExprPtr a) {
    return std::make_shared<MInverse>(std::move(a));
}
MatrixExprPtr mtrace(MatrixExprPtr a) {
    return std::make_shared<MTrace>(std::move(a));
}
MatrixExprPtr mdeterminant(MatrixExprPtr a) {
    return std::make_shared<MDeterminant>(std::move(a));
}
MatrixExprPtr block_matrix(BlockMatrix::BlockGrid blocks) {
    return std::make_shared<BlockMatrix>(std::move(blocks));
}

}  // namespace sympp
