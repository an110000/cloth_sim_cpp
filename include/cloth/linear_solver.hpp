#pragma once

#include "cloth/sparse.hpp"

#include <limits>
#include <string>
#include <vector>

namespace cloth {

struct LinearSolveResult {
    bool success = false;
    std::vector<double> solution;
    // Infinity-norm backward error: ||Ax-b|| / (||A||*||x||+||b||).
    double relative_residual = std::numeric_limits<double>::infinity();
    std::string message;
};

// General square sparse LU, including indefinite systems. No dense conversion.
// Invalid dimensions/non-finite RHS/tolerance throw invalid_argument;
// numerical factorization or residual failure is returned explicitly.
LinearSolveResult solve_sparse_lu(const SparseMatrix& matrix,
                                  const std::vector<double>& rhs,
                                  double relative_tolerance = 1.0e-10);

} // namespace cloth
