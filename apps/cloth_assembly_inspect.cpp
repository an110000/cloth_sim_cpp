#include "cloth/assembly.hpp"
#include "cloth/dofs.hpp"
#include "cloth/mass.hpp"

#include <iomanip>
#include <iostream>
#include <numeric>

int main() {
    const auto mesh = cloth::make_rectangular_mesh(1, 1, 1.0, 1.0);
    const cloth::Material material{};
    const auto elastic = cloth::assemble_elastic(mesh, material);
    const auto mass = cloth::assemble_mass_matrix(mesh, 1000.0, material.thickness, cloth::MassMode::consistent);
    const auto row_sums = mass.multiply(std::vector<double>(mass.cols(), 1.0));
    const cloth::DofMap map(mesh.positions.size(), {2, 3});
    const auto reduced = cloth::restrict_to_free(elastic.hessian, map);
    std::cout << std::setprecision(17)
              << "vertices=" << mesh.positions.size() << '\n'
              << "triangles=" << mesh.triangles.size() << '\n'
              << "full_dofs=" << map.full_size() << '\n'
              << "free_dofs=" << map.free_size() << '\n'
              << "elastic_energy_J=" << elastic.energy << '\n'
              << "hessian_shape=" << elastic.hessian.rows() << 'x' << elastic.hessian.cols() << '\n'
              << "hessian_stored_entries=" << elastic.hessian.non_zeros() << '\n'
              << "reduced_shape=" << reduced.rows() << 'x' << reduced.cols() << '\n'
              << "mass_stored_entries=" << mass.non_zeros() << '\n'
              << "mass_00_kg=" << mass.coeff(0, 0) << '\n'
              << "mass_03_kg=" << mass.coeff(0, 3) << '\n'
              << "total_mass_kg=" << std::accumulate(row_sums.begin(), row_sums.end(), 0.0) / 3.0 << '\n'
              << "vertex0_mass_kg=" << row_sums[0] << '\n'
              << "vertex1_mass_kg=" << row_sums[3] << '\n';
}
