#include "cloth/mesh.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

namespace {

void check(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void near(double actual, double expected, const std::string& message) {
    check(std::isfinite(actual) &&
              std::abs(actual - expected) <= 1.0e-12 * std::max(1.0, std::abs(expected)),
          message);
}

template <typename Exception, typename Function>
void expect_throw(Function action, const std::string& message) {
    try {
        action();
    } catch (const Exception&) {
        return;
    }
    throw std::runtime_error(message);
}

void test_topology_and_embedding() {
    const auto mesh = cloth::make_rectangular_mesh(2, 3, 2.0, 3.0);
    check(mesh.reference_positions.size() == 12, "vertex count");
    check(mesh.triangles.size() == 12, "triangle count");
    check(mesh.triangles.front() == cloth::Triangle{0, 1, 4}, "first triangle topology");
    check(mesh.triangles[1] == cloth::Triangle{0, 4, 3}, "second triangle topology");
    check(mesh.triangles.back() == cloth::Triangle{7, 11, 10}, "last triangle topology");
    for (std::size_t i = 0; i < mesh.positions.size(); ++i) {
        near(mesh.positions[i][0], mesh.reference_positions[i][0], "world x = xi");
        near(mesh.positions[i][1], 0.0, "world y = 0");
        near(mesh.positions[i][2], mesh.reference_positions[i][1], "world z = eta");
    }
    check(mesh.reference_positions.back() == cloth::Vec2{2.0, 3.0}, "rectangle endpoint");
}

void test_mass_conservation_and_distribution() {
    const auto mesh = cloth::make_rectangular_mesh(2, 3, 2.0, 3.0);
    const auto masses = cloth::lumped_masses(mesh, 1000.0, 0.001);
    near(std::accumulate(masses.begin(), masses.end(), 0.0), 6.0, "total mass = rho*t*area");
    check(std::all_of(masses.begin(), masses.end(), [](double m) { return m > 0.0; }),
          "all vertex masses positive");

    const auto unit = cloth::make_rectangular_mesh(1, 1, 1.0, 1.0);
    const auto corner = cloth::lumped_masses(unit, 1.0, 1.0);
    near(corner[0], 1.0 / 3.0, "shared diagonal corner receives two triangles");
    near(corner[1], 1.0 / 6.0, "boundary corner receives one triangle");
    near(corner[2], 1.0 / 6.0, "other boundary corner");
    near(corner[3], 1.0 / 3.0, "other shared diagonal corner");
    const auto doubled_density = cloth::lumped_masses(unit, 2.0, 1.0);
    const auto doubled_thickness = cloth::lumped_masses(unit, 1.0, 2.0);
    for (std::size_t i = 0; i < corner.size(); ++i) {
        near(doubled_density[i], 2.0 * corner[i], "density scaling");
        near(doubled_thickness[i], 2.0 * corner[i], "thickness scaling");
    }
}

void test_invalid_input() {
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    expect_throw<std::invalid_argument>([] { cloth::make_rectangular_mesh(0, 1, 1, 1); }, "zero nx accepted");
    expect_throw<std::invalid_argument>([] { cloth::make_rectangular_mesh(1, 0, 1, 1); }, "zero ny accepted");
    for (const double bad : {0.0, -1.0, infinity, nan}) {
        expect_throw<std::invalid_argument>([bad] { cloth::make_rectangular_mesh(1, 1, bad, 1); }, "invalid width accepted");
        expect_throw<std::invalid_argument>([bad] { cloth::make_rectangular_mesh(1, 1, 1, bad); }, "invalid height accepted");
    }
    const auto largest = std::numeric_limits<std::size_t>::max();
    expect_throw<std::length_error>([largest] { cloth::make_rectangular_mesh(largest, 1, 1, 1); }, "dimension addition overflow accepted");
    expect_throw<std::length_error>([largest] { cloth::make_rectangular_mesh(largest / 2, 3, 1, 1); }, "dimension multiplication overflow accepted");

    const auto unit = cloth::make_rectangular_mesh(1, 1, 1.0, 1.0);
    for (const double bad : {0.0, -1.0, infinity, nan}) {
        expect_throw<std::invalid_argument>([&unit, bad] { cloth::lumped_masses(unit, bad, 1); }, "invalid density accepted");
        expect_throw<std::invalid_argument>([&unit, bad] { cloth::lumped_masses(unit, 1, bad); }, "invalid thickness accepted");
    }
    auto malformed = unit;
    malformed.triangles[0][0] = unit.reference_positions.size();
    expect_throw<std::invalid_argument>([&malformed] { cloth::lumped_masses(malformed, 1, 1); }, "invalid index accepted");
    malformed = unit;
    malformed.reference_positions[1] = malformed.reference_positions[0];
    expect_throw<std::invalid_argument>([&malformed] { cloth::lumped_masses(malformed, 1, 1); }, "degenerate triangle accepted");
    malformed = unit;
    malformed.reference_positions[0][0] = nan;
    expect_throw<std::invalid_argument>([&malformed] { cloth::lumped_masses(malformed, 1, 1); }, "nonfinite reference point accepted");
    malformed = unit;
    malformed.positions.pop_back();
    expect_throw<std::invalid_argument>([&malformed] { cloth::lumped_masses(malformed, 1, 1); }, "mismatched vertex arrays accepted");
    malformed = unit;
    malformed.reference_positions.push_back({2.0, 2.0});
    malformed.positions.push_back({2.0, 0.0, 2.0});
    expect_throw<std::invalid_argument>([&malformed] { cloth::lumped_masses(malformed, 1, 1); }, "isolated vertex accepted");
    expect_throw<std::overflow_error>([&unit] { cloth::lumped_masses(unit, std::numeric_limits<double>::max(), 100.0); }, "mass overflow accepted");
    expect_throw<std::invalid_argument>([] { cloth::lumped_masses(cloth::Mesh{}, 1, 1); }, "empty mesh accepted");
}

}  // namespace

int main() {
    try {
        test_topology_and_embedding();
        test_mass_conservation_and_distribution();
        test_invalid_input();
        std::cout << "mesh and lumped mass tests passed\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "mesh test failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
