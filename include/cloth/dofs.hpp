#pragma once

#include "cloth/sparse.hpp"

#include <cstddef>
#include <limits>
#include <vector>

namespace cloth {

class DofMap {
public:
    static constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();

    // Fixing a vertex removes all three coordinates. Duplicate/out-of-range
    // fixed vertices are rejected. Free DOFs retain ascending global order.
    DofMap(std::size_t vertex_count, const std::vector<std::size_t>& fixed_vertices);

    std::size_t full_size() const { return full_to_free_.size(); }
    std::size_t free_size() const { return free_to_full_.size(); }
    const std::vector<std::size_t>& free_to_full() const { return free_to_full_; }
    const std::vector<std::size_t>& full_to_free() const { return full_to_free_; }

private:
    std::vector<std::size_t> free_to_full_;
    std::vector<std::size_t> full_to_free_;
};

std::vector<double> restrict_to_free(const std::vector<double>& full, const DofMap& map);

// Lift a Newton/search increment, filling fixed coordinates with zero.
// This is not a position reconstruction for nonzero prescribed positions.
std::vector<double> lift_from_free(const std::vector<double>& free, const DofMap& map);

// Extract A_FF. For nonzero prescribed unknowns, a caller must separately form
// the corrected right-hand side b_F - A_FC*x_C. No such correction is hidden.
SparseMatrix restrict_to_free(const SparseMatrix& full, const DofMap& map);

} // namespace cloth
