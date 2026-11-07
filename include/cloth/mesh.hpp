#pragma once

#include <array>
#include <cstddef>
#include <vector>

namespace cloth {

using Vec2 = std::array<double, 2>;
using Vec3 = std::array<double, 3>;
using Triangle = std::array<std::size_t, 3>;

struct Mesh {
    std::vector<Vec2> reference_positions;
    std::vector<Vec3> positions;
    std::vector<Triangle> triangles;
};

// nx and ny count cells, so there are (nx + 1) * (ny + 1) vertices.
// Reference coordinates are (xi, eta); world coordinates are (xi, 0, eta).
Mesh make_rectangular_mesh(std::size_t nx, std::size_t ny,
                           double width, double height);

// Each triangle contributes density * thickness * reference_area / 3
// to each of its three vertices. Returned values are scalar masses in kg.
std::vector<double> lumped_masses(const Mesh& mesh, double density,
                                  double thickness);

}  // namespace cloth
