#include "cloth/linear_solver.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

int main() {
    try {
        const auto a = cloth::SparseMatrix::from_triplets(3, 3,
            {{0,1,2}, {1,0,2}, {1,1,-1}, {2,2,3}});
        const auto result = cloth::solve_sparse_lu(a, {-4,4,1.5});
        if (!result.success || std::abs(result.solution[0]-1) > 1e-12 ||
            std::abs(result.solution[1]+2) > 1e-12 ||
            std::abs(result.solution[2]-.5) > 1e-12 || result.relative_residual > 1e-12) {
            throw std::runtime_error("Indefinite system requiring pivoting did not solve");
        }
        const auto singular = cloth::SparseMatrix::from_triplets(2, 2, {{0,0,1}});
        if (cloth::solve_sparse_lu(singular, {1,1}).success) throw std::runtime_error("Singular accepted");
        const auto empty = cloth::solve_sparse_lu({}, {});
        if (!empty.success || !empty.solution.empty()) throw std::runtime_error("Empty system failed");
        bool rejected = false;
        try { (void)cloth::solve_sparse_lu(a, {1,2}); } catch (const std::invalid_argument&) { rejected = true; }
        if (!rejected) throw std::runtime_error("Wrong RHS size accepted");
        rejected = false;
        try { (void)cloth::solve_sparse_lu(a, {1,2,std::numeric_limits<double>::quiet_NaN()}); }
        catch (const std::invalid_argument&) { rejected = true; }
        if (!rejected) throw std::runtime_error("Non-finite RHS accepted");
        std::cout << "Sparse LU pivoting, singular, empty and input checks passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
