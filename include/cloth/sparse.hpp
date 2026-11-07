#pragma once

#include <cstddef>
#include <vector>

namespace cloth {

struct Triplet {
    std::size_t row;
    std::size_t column;
    double value;
};

// Minimal immutable CSR storage for assembly and matrix-vector products.
// It is not a linear solver. Explicit zero entries are retained so a topology's
// sparsity pattern does not depend on a particular deformation.
class SparseMatrix {
public:
    SparseMatrix() = default; // Valid 0 x 0 matrix.

    static SparseMatrix from_triplets(std::size_t rows, std::size_t cols,
                                      std::vector<Triplet> entries);

    std::size_t rows() const noexcept { return rows_; }
    std::size_t cols() const noexcept { return cols_; }
    // Number of stored entries, including explicit zeros.
    std::size_t non_zeros() const noexcept { return values_.size(); }
    const std::vector<std::size_t>& row_offsets() const noexcept { return row_offsets_; }
    const std::vector<std::size_t>& column_indices() const noexcept { return columns_; }
    const std::vector<double>& values() const noexcept { return values_; }

    double coeff(std::size_t row, std::size_t col) const;
    std::vector<double> multiply(const std::vector<double>& x) const;
    // Small-reference/debug convenience; checks dense allocation size first.
    std::vector<double> dense_row_major() const;

private:
    std::size_t rows_ = 0;
    std::size_t cols_ = 0;
    std::vector<std::size_t> row_offsets_{0};
    std::vector<std::size_t> columns_;
    std::vector<double> values_;
};

} // namespace cloth
