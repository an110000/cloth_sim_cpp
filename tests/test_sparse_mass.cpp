#include "cloth/mass.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

namespace {

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, const std::string& message) {
    check(std::isfinite(actual) && std::abs(actual - expected) <=
              2.0e-12 * std::max(1.0, std::abs(expected)), message);
}

template<class Exception, class Function>
void must_throw(Function function, const std::string& message) {
    try { function(); } catch (const Exception&) { return; }
    throw std::runtime_error(message);
}

void sparse_storage_and_product() {
    const auto matrix = cloth::SparseMatrix::from_triplets(4, 3, {
        {2, 2, 4.0}, {0, 2, 3.0}, {0, 0, 1.0}, {0, 2, -3.0},
        {2, 0, -2.0}, {3, 1, 5.0}});
    check(matrix.rows() == 4 && matrix.cols() == 3, "CSR dimensions");
    check(matrix.non_zeros() == 5, "Duplicate must merge, explicit zero retained");
    check(matrix.row_offsets() == std::vector<std::size_t>{0, 2, 2, 4, 5}, "CSR row offsets");
    check(matrix.column_indices() == std::vector<std::size_t>{0, 2, 0, 2, 1}, "Sorted columns");
    near(matrix.coeff(0, 2), 0.0, "Cancelled zero coefficient");
    near(matrix.coeff(1, 1), 0.0, "Missing entry in empty row");
    const auto product = matrix.multiply({2.0, 3.0, 4.0});
    check(product == std::vector<double>{2.0, 0.0, 12.0, 15.0}, "CSR multiplication");
    const auto dense = matrix.dense_row_major();
    check(dense.size() == 12, "Dense matrix shape");
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        double expected = 0.0;
        for (std::size_t col = 0; col < matrix.cols(); ++col) {
            near(dense[row * 3 + col], matrix.coeff(row, col), "Dense coefficient agreement");
            expected += dense[row * 3 + col] * std::vector<double>{2, 3, 4}[col];
        }
        near(product[row], expected, "Sparse versus dense product");
    }
    const auto empty = cloth::SparseMatrix::from_triplets(0, 0, {});
    check(empty.rows() == 0 && empty.row_offsets() == std::vector<std::size_t>{0}, "0x0 CSR");
    check(empty.multiply({}).empty() && empty.dense_row_major().empty(), "0x0 operations");
    const auto no_columns = cloth::SparseMatrix::from_triplets(2, 0, {});
    check(no_columns.multiply({}) == std::vector<double>{0, 0}, "Matrix with no columns");
    const auto no_rows = cloth::SparseMatrix::from_triplets(0, 2, {});
    check(no_rows.multiply({1, 2}).empty(), "Matrix with no rows");
}

void sparse_rejects_errors() {
    const auto nan = std::numeric_limits<double>::quiet_NaN();
    const auto infinity = std::numeric_limits<double>::infinity();
    const auto largest = std::numeric_limits<std::size_t>::max();
    must_throw<std::out_of_range>([] { cloth::SparseMatrix::from_triplets(2, 2, {{2, 0, 1}}); }, "Out-of-range row accepted");
    must_throw<std::out_of_range>([] { cloth::SparseMatrix::from_triplets(2, 2, {{0, 2, 1}}); }, "Out-of-range column accepted");
    must_throw<std::out_of_range>([] { cloth::SparseMatrix::from_triplets(0, 0, {{0, 0, 0}}); }, "Entry in 0x0 accepted");
    for (const double bad : {nan, infinity})
        must_throw<std::invalid_argument>([bad] { cloth::SparseMatrix::from_triplets(1, 1, {{0, 0, bad}}); }, "Non-finite entry accepted");
    must_throw<std::overflow_error>([] {
        const double huge = std::numeric_limits<double>::max();
        cloth::SparseMatrix::from_triplets(1, 1, {{0, 0, huge}, {0, 0, huge}});
    }, "Duplicate sum overflow accepted");
    must_throw<std::length_error>([largest] { cloth::SparseMatrix::from_triplets(largest, 0, {}); }, "Row-offset size overflow accepted");
    must_throw<std::length_error>([largest] { cloth::SparseMatrix::from_triplets(2, largest, {}).dense_row_major(); }, "Dense product overflow accepted");
    must_throw<std::length_error>([largest] { cloth::SparseMatrix::from_triplets(1, largest, {}).dense_row_major(); }, "Dense capacity overflow accepted");
    const auto identity = cloth::SparseMatrix::from_triplets(2, 2, {{0, 0, 1}, {1, 1, 1}});
    must_throw<std::out_of_range>([&] { identity.coeff(2, 0); }, "Invalid coefficient lookup accepted");
    must_throw<std::invalid_argument>([&] { identity.multiply({1}); }, "Dimension mismatch accepted");
    must_throw<std::invalid_argument>([&] { identity.multiply({nan, 1}); }, "Non-finite vector accepted");
    must_throw<std::overflow_error>([] { cloth::SparseMatrix::from_triplets(1, 1, {{0, 0, 2}}).multiply({std::numeric_limits<double>::max()}); }, "Product overflow accepted");
}

void check_mass_properties(const cloth::Mesh& mesh, double density,
                           double thickness, double expected_mass) {
    const auto lumped = cloth::assemble_mass_matrix(mesh, density, thickness, cloth::MassMode::lumped);
    const auto consistent = cloth::assemble_mass_matrix(mesh, density, thickness, cloth::MassMode::consistent);
    const auto vertex_masses = cloth::lumped_masses(mesh, density, thickness);
    const std::size_t n = mesh.positions.size() * 3;
    check(lumped.rows() == n && consistent.rows() == n, "Mass matrix dimensions");
    const auto row_sums = consistent.multiply(std::vector<double>(n, 1.0));
    near(std::accumulate(consistent.values().begin(), consistent.values().end(), 0.0) / 3.0,
         expected_mass, "Consistent total mass counts xyz only once");
    near(std::accumulate(lumped.values().begin(), lumped.values().end(), 0.0) / 3.0,
         expected_mass, "Lumped total mass counts xyz only once");
    for (std::size_t i = 0; i < n; ++i) {
        near(lumped.coeff(i, i), vertex_masses[i / 3], "Scalar mass repeated per coordinate");
        near(row_sums[i], lumped.coeff(i, i), "Consistent row-sum lumping");
        for (std::size_t j = 0; j < n; ++j) {
            near(consistent.coeff(i, j), consistent.coeff(j, i), "Consistent mass symmetry");
            if (i % 3 != j % 3) near(consistent.coeff(i, j), 0.0, "Mass does not mix axes");
            if (i != j) near(lumped.coeff(i, j), 0.0, "Lumped mass diagonal");
        }
    }
    std::vector<double> velocity(n);
    for (std::size_t i = 0; i < n; ++i)
        velocity[i] = (static_cast<double>(i % 7) - 3.0) / 2.0;
    const auto momentum = consistent.multiply(velocity);
    const double kinetic_twice = std::inner_product(velocity.begin(), velocity.end(), momentum.begin(), 0.0);
    check(kinetic_twice > 0.0 && std::isfinite(kinetic_twice), "Consistent kinetic energy positive");
    for (const auto mode : {cloth::MassMode::lumped, cloth::MassMode::consistent}) {
        const auto original = cloth::assemble_mass_matrix(mesh, density, thickness, mode);
        const auto twice_density = cloth::assemble_mass_matrix(mesh, density * 2, thickness, mode);
        const auto twice_thickness = cloth::assemble_mass_matrix(mesh, density, thickness * 2, mode);
        for (std::size_t row = 0; row < n; ++row)
            for (std::size_t col = 0; col < n; ++col) {
                near(twice_density.coeff(row, col), 2.0 * original.coeff(row, col), "Density scaling");
                near(twice_thickness.coeff(row, col), 2.0 * original.coeff(row, col), "Thickness scaling");
            }
    }
}

void mass_analytic_cases() {
    const auto unit = cloth::make_rectangular_mesh(1, 1, 1, 1);
    check_mass_properties(unit, 1, 1, 1);
    const auto mass = cloth::assemble_mass_matrix(unit, 1, 1, cloth::MassMode::consistent);
    near(mass.coeff(0, 0), 1.0 / 6.0, "Shared corner diagonal");
    near(mass.coeff(3, 3), 1.0 / 12.0, "Single-triangle corner diagonal");
    near(mass.coeff(0, 9), 1.0 / 12.0, "Shared diagonal off-diagonal");
    near(mass.coeff(0, 3), 1.0 / 24.0, "Single-edge off-diagonal");
    near(mass.coeff(3, 6), 0.0, "Unconnected vertices have zero coefficient");
    auto irregular = unit;
    irregular.reference_positions = {{0, 0}, {2, 0}, {0, 1}, {1, 1}};
    check_mass_properties(irregular, 1, 1, 1.5);
    const auto irregular_lumped = cloth::assemble_mass_matrix(irregular, 1, 1, cloth::MassMode::lumped);
    near(irregular_lumped.coeff(0, 0), 0.5, "Unequal triangle contributions accumulate");
    near(irregular_lumped.coeff(3, 3), 1.0 / 3.0, "Large triangle corner");
    near(irregular_lumped.coeff(6, 6), 1.0 / 6.0, "Small triangle corner");
    check_mass_properties(cloth::make_rectangular_mesh(2, 3, 2, 3), 1000, 0.001, 6);
}

void mass_rejects_errors() {
    const auto mesh = cloth::make_rectangular_mesh(1, 1, 1, 1);
    for (const auto mode : {cloth::MassMode::lumped, cloth::MassMode::consistent}) {
        for (const double bad : {0.0, -1.0, std::numeric_limits<double>::infinity(),
                                 std::numeric_limits<double>::quiet_NaN()}) {
            must_throw<std::invalid_argument>([&] { cloth::assemble_mass_matrix(mesh, bad, 1, mode); }, "Invalid density accepted");
            must_throw<std::invalid_argument>([&] { cloth::assemble_mass_matrix(mesh, 1, bad, mode); }, "Invalid thickness accepted");
        }
        auto malformed = mesh;
        malformed.triangles[0][0] = 99;
        must_throw<std::invalid_argument>([&] { cloth::assemble_mass_matrix(malformed, 1, 1, mode); }, "Invalid index accepted");
        malformed = mesh;
        malformed.reference_positions[1] = malformed.reference_positions[0];
        must_throw<std::invalid_argument>([&] { cloth::assemble_mass_matrix(malformed, 1, 1, mode); }, "Degenerate element accepted");
        malformed = mesh;
        malformed.positions.push_back({2, 0, 2});
        malformed.reference_positions.push_back({2, 2});
        must_throw<std::invalid_argument>([&] { cloth::assemble_mass_matrix(malformed, 1, 1, mode); }, "Isolated vertex accepted");
        must_throw<std::overflow_error>([&] { cloth::assemble_mass_matrix(mesh, std::numeric_limits<double>::max(), 100, mode); }, "Mass overflow accepted");
    }
    must_throw<std::invalid_argument>([&] { cloth::assemble_mass_matrix(mesh, 1, 1, static_cast<cloth::MassMode>(99)); }, "Invalid mass mode accepted");
}

} // namespace

int main() {
    try {
        sparse_storage_and_product();
        sparse_rejects_errors();
        mass_analytic_cases();
        mass_rejects_errors();
        std::cout << "CSR storage/product and lumped/consistent mass tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Sparse/mass test failed: " << error.what() << '\n';
        return 1;
    }
}
