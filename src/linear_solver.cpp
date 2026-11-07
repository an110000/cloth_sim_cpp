#include "cloth/linear_solver.hpp"

#include <Eigen/SparseLU>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace cloth {

LinearSolveResult solve_sparse_lu(const SparseMatrix& matrix,
                                  const std::vector<double>& rhs,
                                  double relative_tolerance) {
    if (matrix.rows() != matrix.cols() || rhs.size() != matrix.rows() ||
        !std::isfinite(relative_tolerance) || relative_tolerance <= 0.0) {
        throw std::invalid_argument("Sparse LU requires a square matrix, matching RHS and positive tolerance");
    }
    for (const auto value : rhs) {
        if (!std::isfinite(value)) throw std::invalid_argument("Sparse LU RHS must be finite");
    }
    if (rhs.empty()) return {true, {}, 0.0, "empty system"};
    using StorageIndex = int;
    if (matrix.rows() > static_cast<std::size_t>(std::numeric_limits<StorageIndex>::max()) ||
        matrix.non_zeros() > static_cast<std::size_t>(std::numeric_limits<StorageIndex>::max())) {
        throw std::length_error("Sparse LU exceeds Eigen storage index range");
    }
    using EigenMatrix = Eigen::SparseMatrix<double, Eigen::ColMajor, StorageIndex>;
    std::vector<Eigen::Triplet<double, StorageIndex>> entries;
    entries.reserve(matrix.non_zeros());
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        for (auto k = matrix.row_offsets()[row]; k < matrix.row_offsets()[row + 1]; ++k) {
            entries.emplace_back(static_cast<StorageIndex>(row),
                                 static_cast<StorageIndex>(matrix.column_indices()[k]), matrix.values()[k]);
        }
    }
    const auto n = static_cast<StorageIndex>(matrix.rows());
    EigenMatrix system(n, n);
    // CSR input is ordered by row; Eigen's compressed column-major constructor
    // requires outer (column), then inner (row) order for this direct path.
    std::sort(entries.begin(), entries.end(), [](const auto& a, const auto& b) {
        return a.col() != b.col() ? a.col() < b.col() : a.row() < b.row();
    });
    system.setFromSortedTriplets(entries.begin(), entries.end());
    system.makeCompressed();
    Eigen::SparseLU<EigenMatrix, Eigen::COLAMDOrdering<StorageIndex>> lu;
    lu.compute(system);
    if (lu.info() != Eigen::Success) {
        return {false, {}, std::numeric_limits<double>::infinity(), "Sparse LU factorization failed"};
    }
    const Eigen::Map<const Eigen::VectorXd> b(rhs.data(), n);
    const Eigen::VectorXd x = lu.solve(b);
    if (lu.info() != Eigen::Success || !x.allFinite()) {
        return {false, {}, std::numeric_limits<double>::infinity(), "Sparse LU solution is not finite"};
    }
    std::vector<double> solution(x.data(), x.data() + x.size());
    long double matrix_norm = 0.0L, x_norm = 0.0L, b_norm = 0.0L, error_norm = 0.0L;
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        long double ax = 0.0L, row_sum = 0.0L;
        for (auto k = matrix.row_offsets()[row]; k < matrix.row_offsets()[row + 1]; ++k) {
            ax += static_cast<long double>(matrix.values()[k]) * solution[matrix.column_indices()[k]];
            row_sum += std::abs(static_cast<long double>(matrix.values()[k]));
        }
        const long double difference = std::abs(ax - rhs[row]);
        if (!std::isfinite(ax) || !std::isfinite(row_sum) || !std::isfinite(difference)) {
            return {false, {}, std::numeric_limits<double>::infinity(), "Sparse LU residual arithmetic overflow"};
        }
        matrix_norm = std::max(matrix_norm, row_sum);
        x_norm = std::max(x_norm, std::abs(static_cast<long double>(solution[row])));
        b_norm = std::max(b_norm, std::abs(static_cast<long double>(rhs[row])));
        error_norm = std::max(error_norm, difference);
    }
    const long double denominator = matrix_norm * x_norm + b_norm;
    if (!std::isfinite(denominator)) {
        return {false, {}, std::numeric_limits<double>::infinity(), "Sparse LU residual scale overflow"};
    }
    const double error = static_cast<double>(denominator == 0.0L ? error_norm : error_norm / denominator);
    if (!std::isfinite(error) || error > relative_tolerance) {
        return {false, {}, error, "Sparse LU backward error exceeds tolerance"};
    }
    return {true, std::move(solution), error, "success"};
}

} // namespace cloth
