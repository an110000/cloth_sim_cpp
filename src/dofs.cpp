#include "cloth/dofs.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace cloth {
namespace {

void require_finite(const std::vector<double>& values) {
    for (const auto value : values) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("DOF vector must be finite");
        }
    }
}

} // namespace

DofMap::DofMap(std::size_t vertex_count, const std::vector<std::size_t>& fixed_vertices) {
    if (vertex_count > std::numeric_limits<std::size_t>::max() / 3) {
        throw std::length_error("DOF count overflows size_t");
    }
    const auto ndof = 3 * vertex_count;
    if (ndof > full_to_free_.max_size() || ndof > free_to_full_.max_size()) {
        throw std::length_error("DOF count exceeds vector capacity");
    }
    std::vector<bool> fixed(vertex_count, false);
    for (const auto vertex : fixed_vertices) {
        if (vertex >= vertex_count) {
            throw std::invalid_argument("Fixed vertex is out of range");
        }
        if (fixed[vertex]) {
            throw std::invalid_argument("Duplicate fixed vertex");
        }
        fixed[vertex] = true;
    }
    full_to_free_.assign(ndof, npos);
    free_to_full_.reserve(3 * (vertex_count - fixed_vertices.size()));
    for (std::size_t vertex = 0; vertex < vertex_count; ++vertex) {
        if (fixed[vertex]) continue;
        for (std::size_t axis = 0; axis < 3; ++axis) {
            const auto dof = 3 * vertex + axis;
            full_to_free_[dof] = free_to_full_.size();
            free_to_full_.push_back(dof);
        }
    }
}

std::vector<double> restrict_to_free(const std::vector<double>& full, const DofMap& map) {
    if (full.size() != map.full_size()) {
        throw std::invalid_argument("Full vector size does not match DOF map");
    }
    require_finite(full);
    std::vector<double> free(map.free_size());
    for (std::size_t i = 0; i < free.size(); ++i) {
        free[i] = full[map.free_to_full()[i]];
    }
    return free;
}

std::vector<double> lift_from_free(const std::vector<double>& free, const DofMap& map) {
    if (free.size() != map.free_size()) {
        throw std::invalid_argument("Free vector size does not match DOF map");
    }
    require_finite(free);
    std::vector<double> full(map.full_size(), 0.0);
    for (std::size_t i = 0; i < free.size(); ++i) {
        full[map.free_to_full()[i]] = free[i];
    }
    return full;
}

SparseMatrix restrict_to_free(const SparseMatrix& full, const DofMap& map) {
    if (full.rows() != map.full_size() || full.cols() != map.full_size()) {
        throw std::invalid_argument("Matrix shape does not match DOF map");
    }
    if (map.free_size() == 0) {
        return SparseMatrix::from_triplets(0, 0, {});
    }
    std::vector<Triplet> entries;
    entries.reserve(full.non_zeros());
    for (std::size_t row = 0; row < full.rows(); ++row) {
        const auto free_row = map.full_to_free()[row];
        if (free_row == DofMap::npos) continue;
        for (auto i = full.row_offsets()[row]; i < full.row_offsets()[row + 1]; ++i) {
            const auto free_column = map.full_to_free()[full.column_indices()[i]];
            if (free_column != DofMap::npos) {
                entries.push_back({free_row, free_column, full.values()[i]});
            }
        }
    }
    return SparseMatrix::from_triplets(map.free_size(), map.free_size(), std::move(entries));
}

} // namespace cloth
