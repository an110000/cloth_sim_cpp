#include "cloth/assembly.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace cloth {

ElasticAssembly assemble_elastic(const Mesh& mesh, const Material& material) {
    const std::size_t vertices = mesh.reference_positions.size();
    if (vertices == 0 || mesh.positions.size() != vertices || mesh.triangles.empty()) {
        throw std::invalid_argument("Elastic assembly requires matching vertices and triangles");
    }
    if (vertices > std::numeric_limits<std::size_t>::max() / 3) {
        throw std::length_error("Elastic assembly DOF count overflows size_t");
    }
    const std::size_t ndof = 3 * vertices;
    if (ndof > std::vector<double>{}.max_size()) {
        throw std::length_error("Elastic assembly gradient exceeds vector capacity");
    }
    if (mesh.triangles.size() > std::vector<Triplet>{}.max_size() / 81) {
        throw std::length_error("Elastic assembly triplet count exceeds vector capacity");
    }
    for (const auto& point : mesh.reference_positions) {
        for (const auto value : point) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument("Reference positions must be finite");
            }
        }
    }
    for (const auto& point : mesh.positions) {
        for (const auto value : point) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument("Current positions must be finite");
            }
        }
    }
    std::vector<bool> referenced(vertices, false);
    for (const auto& triangle : mesh.triangles) {
        for (const auto vertex : triangle) {
            if (vertex >= vertices) {
                throw std::invalid_argument("Elastic assembly triangle index is out of range");
            }
            referenced[vertex] = true;
        }
    }
    if (std::find(referenced.begin(), referenced.end(), false) != referenced.end()) {
        throw std::invalid_argument("Elastic assembly mesh contains an isolated vertex");
    }

    double energy = 0.0;
    std::vector<double> gradient(ndof, 0.0);
    std::vector<Triplet> triplets;
    triplets.reserve(mesh.triangles.size() * 81);
    for (const auto& ids : mesh.triangles) {
        std::array<Vec2, 3> reference{};
        std::array<Vec3, 3> positions{};
        std::array<std::size_t, 9> dofs{};
        for (std::size_t vertex = 0; vertex < 3; ++vertex) {
            reference[vertex] = mesh.reference_positions[ids[vertex]];
            positions[vertex] = mesh.positions[ids[vertex]];
            for (std::size_t axis = 0; axis < 3; ++axis) {
                dofs[3 * vertex + axis] = 3 * ids[vertex] + axis;
            }
        }
        const auto local = evaluate_triangle_with_hessian(reference, positions, material);
        energy += local.value.energy;
        if (!std::isfinite(energy)) {
            throw std::overflow_error("Assembled elastic energy is not finite");
        }
        for (std::size_t row = 0; row < 9; ++row) {
            gradient[dofs[row]] += local.value.gradient[row / 3][row % 3];
            if (!std::isfinite(gradient[dofs[row]])) {
                throw std::overflow_error("Assembled elastic gradient is not finite");
            }
            for (std::size_t column = 0; column < 9; ++column) {
                triplets.push_back({dofs[row], dofs[column], local.hessian[row][column]});
            }
        }
    }
    return {energy, std::move(gradient),
            SparseMatrix::from_triplets(ndof, ndof, std::move(triplets))};
}

} // namespace cloth
