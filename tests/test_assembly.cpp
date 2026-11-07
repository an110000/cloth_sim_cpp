#include "cloth/assembly.hpp"
#include "cloth/dofs.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const std::string& message) {
    check(std::isfinite(actual) && std::isfinite(expected) &&
              std::abs(actual - expected) <= tolerance, message);
}

template<class F> void must_throw(F action, const std::string& message) {
    bool caught = false;
    try { action(); } catch (const std::exception&) { caught = true; }
    check(caught, message);
}

cloth::Mesh deformed_mesh() {
    auto mesh = cloth::make_rectangular_mesh(2, 1, 1.3, 0.7);
    // Nonuniform, skewed reference triangles exercise inverse transpose and
    // shared-vertex indexing; a regular axis-aligned triangle can hide errors.
    mesh.reference_positions[1] = {0.47, -0.08};
    mesh.reference_positions[4] = {0.79, 0.76};
    for (std::size_t i = 0; i < mesh.positions.size(); ++i) {
        const auto& x = mesh.reference_positions[i];
        mesh.positions[i] = {1.08 * x[0] + 0.14 * x[1],
                             0.05 * x[0] * x[1] + 0.01 * static_cast<double>(i),
                             0.93 * x[1] - 0.06 * x[0]};
    }
    return mesh;
}

// Independent scalar energy sum, without the global scatter implementation.
double direct_energy(const cloth::Mesh& mesh, const cloth::Material& material) {
    double sum = 0.0;
    for (const auto& triangle : mesh.triangles) {
        const std::array<cloth::Vec2, 3> reference{
            mesh.reference_positions[triangle[0]], mesh.reference_positions[triangle[1]],
            mesh.reference_positions[triangle[2]]};
        const std::array<cloth::Vec3, 3> positions{
            mesh.positions[triangle[0]], mesh.positions[triangle[1]], mesh.positions[triangle[2]]};
        sum += cloth::evaluate_triangle(reference, positions, material).energy;
    }
    return sum;
}

void test_global_derivatives() {
    const auto mesh = deformed_mesh();
    const cloth::Material material{};
    const auto assembled = cloth::assemble_elastic(mesh, material);
    const auto n = 3 * mesh.positions.size();
    check(assembled.gradient.size() == n, "Global gradient size");
    check(assembled.hessian.rows() == n && assembled.hessian.cols() == n, "Global Hessian shape");
    near(assembled.energy, direct_energy(mesh, material), 1.0e-11, "Global energy sum");
    const double step = 1.0e-6;
    double error_g2 = 0.0, norm_g2 = 0.0;
    double error_h2 = 0.0, norm_h2 = 0.0;
    double max_h = 0.0, max_asymmetry = 0.0;
    for (std::size_t column = 0; column < n; ++column) {
        auto plus = mesh;
        auto minus = mesh;
        plus.positions[column / 3][column % 3] += step;
        minus.positions[column / 3][column % 3] -= step;
        const double numeric_g = (direct_energy(plus, material) - direct_energy(minus, material)) /
                                 (2.0 * step);
        const double dg = numeric_g - assembled.gradient[column];
        error_g2 += dg * dg;
        norm_g2 += numeric_g * numeric_g;
        const auto gp = cloth::assemble_elastic(plus, material).gradient;
        const auto gm = cloth::assemble_elastic(minus, material).gradient;
        for (std::size_t row = 0; row < n; ++row) {
            const double numeric_h = (gp[row] - gm[row]) / (2.0 * step);
            const double actual_h = assembled.hessian.coeff(row, column);
            const double dh = numeric_h - actual_h;
            error_h2 += dh * dh;
            norm_h2 += numeric_h * numeric_h;
            max_h = std::max(max_h, std::abs(actual_h));
            max_asymmetry = std::max(max_asymmetry,
                std::abs(actual_h - assembled.hessian.coeff(column, row)));
        }
    }
    check(std::sqrt(error_g2) / std::max(1.0, std::sqrt(norm_g2)) < 2.0e-7,
          "Assembled gradient finite difference");
    check(std::sqrt(error_h2) / std::max(1.0, std::sqrt(norm_h2)) < 2.0e-7,
          "Assembled Hessian finite difference");
    check(max_asymmetry < 1.0e-11 * std::max(1.0, max_h), "Raw assembled symmetry");

    for (std::size_t axis = 0; axis < 3; ++axis) {
        std::vector<double> translation(n, 0.0);
        double net_gradient = 0.0;
        for (std::size_t vertex = 0; vertex < mesh.positions.size(); ++vertex) {
            translation[3 * vertex + axis] = 1.0;
            net_gradient += assembled.gradient[3 * vertex + axis];
        }
        near(net_gradient, 0.0, 1.0e-10, "Global internal force balance");
        const auto action = assembled.hessian.multiply(translation);
        for (const auto value : action)
            near(value, 0.0, 1.0e-10 * std::max(1.0, max_h), "Hessian translation null mode");
    }

    auto reordered = mesh;
    std::reverse(reordered.triangles.begin(), reordered.triangles.end());
    const auto other = cloth::assemble_elastic(reordered, material);
    near(other.energy, assembled.energy, 1.0e-10, "Element order energy invariance");
    for (std::size_t row = 0; row < n; ++row) {
        near(other.gradient[row], assembled.gradient[row], 1.0e-10, "Element order gradient");
        for (std::size_t column = 0; column < n; ++column)
            near(other.hessian.coeff(row, column), assembled.hessian.coeff(row, column),
                  1.0e-10 * std::max(1.0, max_h), "Element order Hessian");
    }
    auto rest = mesh;
    for (std::size_t vertex = 0; vertex < rest.positions.size(); ++vertex) {
        rest.positions[vertex] = {rest.reference_positions[vertex][0], 0.0,
                                  rest.reference_positions[vertex][1]};
    }
    const auto rest_hessian = cloth::assemble_elastic(rest, material).hessian;
    check(rest_hessian.row_offsets() == assembled.hessian.row_offsets() &&
          rest_hessian.column_indices() == assembled.hessian.column_indices(),
          "CSR topology pattern must be independent of the current deformation");
    // Assembly is a read-only operation, not a time step.
    check(mesh.positions == deformed_mesh().positions, "Assembly mutated input positions");
}

void test_rest_patch_and_shared_scatter() {
    const auto mesh = cloth::make_rectangular_mesh(2, 2, 1.0, 1.0);
    const auto assembly = cloth::assemble_elastic(mesh, {});
    near(assembly.energy, 0.0, 1.0e-12, "Rest patch energy");
    for (const auto value : assembly.gradient) near(value, 0.0, 1.0e-10, "Rest patch gradient");

    auto double_patch = mesh;
    double_patch.triangles.insert(double_patch.triangles.end(), mesh.triangles.begin(), mesh.triangles.end());
    const auto doubled = cloth::assemble_elastic(double_patch, {});
    for (std::size_t row = 0; row < assembly.hessian.rows(); ++row)
        for (std::size_t column = 0; column < assembly.hessian.cols(); ++column)
            near(doubled.hessian.coeff(row, column), 2.0 * assembly.hessian.coeff(row, column),
                  1.0e-9, "Repeated triangles must accumulate, not overwrite");
}

void test_dof_mapping() {
    const cloth::DofMap map(4, {2, 0});
    check(map.full_size() == 12 && map.free_size() == 6, "DOF sizes");
    const std::vector<std::size_t> expected{3, 4, 5, 9, 10, 11};
    check(map.free_to_full() == expected, "Ascending free DOF order");
    check(map.full_to_free()[0] == cloth::DofMap::npos &&
          map.full_to_free()[6] == cloth::DofMap::npos, "Fixed DOF sentinel");
    check(map.full_to_free()[9] == 3, "Inverse DOF map");
    std::vector<double> full(12);
    for (std::size_t i = 0; i < full.size(); ++i) full[i] = static_cast<double>(i + 1);
    const auto reduced = cloth::restrict_to_free(full, map);
    check(reduced == std::vector<double>{4, 5, 6, 10, 11, 12}, "Restrict vector");
    const auto lifted = cloth::lift_from_free(reduced, map);
    for (std::size_t i = 0; i < full.size(); ++i)
        near(lifted[i], map.full_to_free()[i] == cloth::DofMap::npos ? 0.0 : full[i],
              0.0, "Lift increment fixes zero");

    std::vector<cloth::Triplet> entries;
    for (std::size_t row = 0; row < 12; ++row)
        for (std::size_t column = 0; column < 12; ++column)
            entries.push_back({row, column, static_cast<double>(100 * row + column)});
    const auto matrix = cloth::SparseMatrix::from_triplets(12, 12, entries);
    const auto block = cloth::restrict_to_free(matrix, map);
    for (std::size_t row = 0; row < expected.size(); ++row)
        for (std::size_t column = 0; column < expected.size(); ++column)
            near(block.coeff(row, column), static_cast<double>(100 * expected[row] + expected[column]),
                  0.0, "A_FF analytic indices");

    const cloth::DofMap all_fixed(4, {0, 1, 2, 3});
    check(cloth::restrict_to_free(full, all_fixed).empty(), "All-fixed vector");
    check(cloth::lift_from_free({}, all_fixed) == std::vector<double>(12, 0.0), "All-fixed lift");
    const auto empty = cloth::restrict_to_free(matrix, all_fixed);
    check(empty.rows() == 0 && empty.cols() == 0 && empty.non_zeros() == 0, "All-fixed matrix");
    const cloth::DofMap no_fixed(4, {});
    check(cloth::restrict_to_free(full, no_fixed) == full, "Identity map");
    const cloth::DofMap no_vertices(0, {});
    check(no_vertices.full_size() == 0 && no_vertices.free_size() == 0, "Empty DOF map");

    must_throw([] { cloth::DofMap invalid(4, {4}); }, "Out-of-range fixed vertex accepted");
    must_throw([] { cloth::DofMap invalid(4, {1, 1}); }, "Duplicate fixed vertex accepted");
    must_throw([] { cloth::DofMap invalid(std::numeric_limits<std::size_t>::max(), {}); },
               "Overflow DOF count accepted");
    must_throw([&] { cloth::restrict_to_free(std::vector<double>{1}, map); }, "Wrong full size accepted");
    must_throw([&] { cloth::lift_from_free(std::vector<double>{1}, map); }, "Wrong free size accepted");
    auto nan_vector = full;
    nan_vector[0] = std::numeric_limits<double>::quiet_NaN();
    must_throw([&] { cloth::restrict_to_free(nan_vector, map); }, "Nonfinite full vector accepted");
    must_throw([&] { cloth::restrict_to_free(cloth::SparseMatrix::from_triplets(2, 2, {}), map); },
               "Wrong matrix shape accepted");
}

void test_invalid_mesh() {
    must_throw([] { cloth::assemble_elastic({}, {}); }, "Empty mesh accepted");
    auto mesh = deformed_mesh();
    mesh.positions.pop_back();
    must_throw([&] { cloth::assemble_elastic(mesh, {}); }, "Position mismatch accepted");
    mesh = deformed_mesh();
    mesh.triangles[0][0] = mesh.positions.size();
    must_throw([&] { cloth::assemble_elastic(mesh, {}); }, "Bad triangle index accepted");
    mesh = deformed_mesh();
    mesh.reference_positions.push_back({2.0, 2.0});
    mesh.positions.push_back({2.0, 0.0, 2.0});
    must_throw([&] { cloth::assemble_elastic(mesh, {}); }, "Isolated vertex accepted");
    mesh = deformed_mesh();
    mesh.reference_positions[0][0] = std::numeric_limits<double>::infinity();
    must_throw([&] { cloth::assemble_elastic(mesh, {}); }, "Infinite reference accepted");
    mesh = deformed_mesh();
    mesh.positions[0][1] = std::numeric_limits<double>::quiet_NaN();
    must_throw([&] { cloth::assemble_elastic(mesh, {}); }, "NaN position accepted");
    mesh = deformed_mesh();
    mesh.triangles[0][1] = mesh.triangles[0][0];
    must_throw([&] { cloth::assemble_elastic(mesh, {}); }, "Degenerate triangle accepted");
    mesh = deformed_mesh();
    for (auto& point : mesh.positions) point = {0.0, 0.0, 0.0};
    must_throw([&] { cloth::assemble_elastic(mesh, {}); }, "Collapsed current mesh accepted");
}

} // namespace

int main() {
    try {
        test_global_derivatives();
        test_rest_patch_and_shared_scatter();
        test_dof_mapping();
        test_invalid_mesh();
        std::cout << "Global assembly derivatives, scatter, DOF restriction and invalid inputs passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
