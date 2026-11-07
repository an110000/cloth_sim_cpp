#include "cloth/mass.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace cloth {

SparseMatrix assemble_mass_matrix(const Mesh& mesh, double density,
                                  double thickness, MassMode mode) {
    if (mode != MassMode::lumped && mode != MassMode::consistent) {
        throw std::invalid_argument("Unsupported mass mode");
    }
    // Reuse the existing validation of parameters, reference geometry, indices,
    // positive masses, matching vertex arrays, and isolated vertices.
    const auto vertex_masses = lumped_masses(mesh, density, thickness);
    if (vertex_masses.size() > std::numeric_limits<std::size_t>::max() / 3) {
        throw std::length_error("Mass degrees of freedom overflow size_t");
    }
    const std::size_t ndof = 3 * vertex_masses.size();
    std::vector<Triplet> entries;
    if (mode == MassMode::lumped) {
        if (ndof > entries.max_size()) {
            throw std::length_error("Lumped mass entries exceed vector capacity");
        }
        entries.reserve(ndof);
        for (std::size_t vertex = 0; vertex < vertex_masses.size(); ++vertex) {
            for (std::size_t axis = 0; axis < 3; ++axis) {
                const auto dof = 3 * vertex + axis;
                entries.push_back({dof, dof, vertex_masses[vertex]});
            }
        }
        return SparseMatrix::from_triplets(ndof, ndof, std::move(entries));
    }
    // Nine scalar vertex pairs, each contributing three matching-axis entries.
    if (mesh.triangles.size() > entries.max_size() / 27) {
        throw std::length_error("Consistent mass entries exceed vector capacity");
    }
    entries.reserve(27 * mesh.triangles.size());
    for (const auto& triangle : mesh.triangles) {
        const auto& a = mesh.reference_positions[triangle[0]];
        const auto& b = mesh.reference_positions[triangle[1]];
        const auto& c = mesh.reference_positions[triangle[2]];
        const long double bx = static_cast<long double>(b[0]) - a[0];
        const long double by = static_cast<long double>(b[1]) - a[1];
        const long double cx = static_cast<long double>(c[0]) - a[0];
        const long double cy = static_cast<long double>(c[1]) - a[1];
        const long double twice_area = std::abs(bx * cy - by * cx);
        const double base = static_cast<double>((twice_area / 24.0L) * density * thickness);
        if (!std::isfinite(base) || base <= 0.0 || !std::isfinite(2.0 * base)) {
            throw std::overflow_error("Consistent mass entry is outside double precision range");
        }
        for (std::size_t local_a = 0; local_a < 3; ++local_a) {
            for (std::size_t local_b = 0; local_b < 3; ++local_b) {
                const double value = (local_a == local_b ? 2.0 : 1.0) * base;
                for (std::size_t axis = 0; axis < 3; ++axis) {
                    entries.push_back({3 * triangle[local_a] + axis,
                                       3 * triangle[local_b] + axis, value});
                }
            }
        }
    }
    return SparseMatrix::from_triplets(ndof, ndof, std::move(entries));
}

} // namespace cloth
