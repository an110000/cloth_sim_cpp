#include "cloth/sparse.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cloth {

SparseMatrix SparseMatrix::from_triplets(std::size_t rows, std::size_t cols,
                                         std::vector<Triplet> entries) {
    SparseMatrix result;
    if (rows == std::numeric_limits<std::size_t>::max() ||
        rows + 1 > result.row_offsets_.max_size()) {
        throw std::length_error("CSR row offsets exceed storage capacity");
    }
    if (entries.size() > result.columns_.max_size() ||
        entries.size() > result.values_.max_size()) {
        throw std::length_error("CSR entries exceed storage capacity");
    }
    for (const auto& entry : entries) {
        if (entry.row >= rows || entry.column >= cols) {
            throw std::out_of_range("Sparse triplet index is outside matrix dimensions");
        }
        if (!std::isfinite(entry.value)) {
            throw std::invalid_argument("Sparse triplet value must be finite");
        }
    }
    // Keep original order among equal indices for deterministic duplicate sums.
    std::stable_sort(entries.begin(), entries.end(), [](const Triplet& a, const Triplet& b) {
        return a.row < b.row || (a.row == b.row && a.column < b.column);
    });
    result.rows_ = rows;
    result.cols_ = cols;
    result.row_offsets_.assign(rows + 1, 0);
    result.columns_.reserve(entries.size());
    result.values_.reserve(entries.size());
    std::size_t index = 0;
    while (index < entries.size()) {
        const auto row = entries[index].row;
        const auto column = entries[index].column;
        double value = entries[index].value;
        ++index;
        while (index < entries.size() && entries[index].row == row &&
               entries[index].column == column) {
            value += entries[index].value;
            if (!std::isfinite(value)) {
                throw std::overflow_error("Duplicate sparse entries overflow double precision");
            }
            ++index;
        }
        result.columns_.push_back(column);
        result.values_.push_back(value); // Keep even an exactly cancelled zero.
        ++result.row_offsets_[row + 1];
    }
    for (std::size_t row = 0; row < rows; ++row) {
        result.row_offsets_[row + 1] += result.row_offsets_[row];
    }
    return result;
}

double SparseMatrix::coeff(std::size_t row, std::size_t col) const {
    if (row >= rows_ || col >= cols_) {
        throw std::out_of_range("Sparse coefficient index is outside matrix dimensions");
    }
    const auto first = columns_.begin() + static_cast<std::ptrdiff_t>(row_offsets_[row]);
    const auto last = columns_.begin() + static_cast<std::ptrdiff_t>(row_offsets_[row + 1]);
    const auto found = std::lower_bound(first, last, col);
    if (found == last || *found != col) return 0.0;
    return values_[static_cast<std::size_t>(found - columns_.begin())];
}

std::vector<double> SparseMatrix::multiply(const std::vector<double>& x) const {
    if (x.size() != cols_) {
        throw std::invalid_argument("Sparse matrix-vector dimension mismatch");
    }
    for (const double value : x) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("Matrix-vector input must be finite");
        }
    }
    std::vector<double> result(rows_, 0.0);
    for (std::size_t row = 0; row < rows_; ++row) {
        double sum = 0.0;
        for (std::size_t k = row_offsets_[row]; k < row_offsets_[row + 1]; ++k) {
            const double term = values_[k] * x[columns_[k]];
            if (!std::isfinite(term)) {
                throw std::overflow_error("Matrix-vector product term overflows double precision");
            }
            sum += term;
            if (!std::isfinite(sum)) {
                throw std::overflow_error("Matrix-vector row sum overflows double precision");
            }
        }
        result[row] = sum;
    }
    return result;
}

std::vector<double> SparseMatrix::dense_row_major() const {
    if (cols_ != 0 && rows_ > std::numeric_limits<std::size_t>::max() / cols_) {
        throw std::length_error("Dense matrix dimensions overflow size_t");
    }
    const std::size_t count = rows_ * cols_;
    std::vector<double> result;
    if (count > result.max_size()) {
        throw std::length_error("Dense matrix exceeds vector capacity");
    }
    result.assign(count, 0.0);
    for (std::size_t row = 0; row < rows_; ++row) {
        for (std::size_t k = row_offsets_[row]; k < row_offsets_[row + 1]; ++k) {
            result[row * cols_ + columns_[k]] = values_[k];
        }
    }
    return result;
}

} // namespace cloth
