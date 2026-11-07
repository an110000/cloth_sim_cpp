#include "cloth/material.hpp"

#include <array>
#include <iomanip>
#include <iostream>
#include <numeric>

int main() {
    const cloth::Material material{};
    const auto mesh = cloth::make_rectangular_mesh(2, 2, 1.0, 1.0);
    const auto mass = cloth::lumped_masses(mesh, 1000.0, material.thickness);
    const auto& ids = mesh.triangles.front();
    std::array<cloth::Vec2, 3> reference{};
    std::array<cloth::Vec3, 3> positions{};
    for (std::size_t i = 0; i < 3; ++i) {
        reference[i] = mesh.reference_positions[ids[i]];
        positions[i] = mesh.positions[ids[i]];
    }
    const auto element = cloth::evaluate_triangle(reference, positions, material);
    std::cout << std::setprecision(17)
              << "vertices=" << mesh.positions.size() << '\n'
              << "triangles=" << mesh.triangles.size() << '\n'
              << "total_mass_kg=" << std::accumulate(mass.begin(), mass.end(), 0.0) << '\n'
              << "rest_element_energy_J=" << element.energy << '\n'
              << "rest_area_ratio=" << element.area_ratio << '\n';
}
