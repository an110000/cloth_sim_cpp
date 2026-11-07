#include "cloth/mesh.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace cloth {
namespace {

void require_positive_finite(double value, const char* message) {
    if (!std::isfinite(value) || value <= 0.0) {
        throw std::invalid_argument(message);
    }
}

std::size_t checked_product(std::size_t a, std::size_t b) {
    if (b != 0 && a > std::numeric_limits<std::size_t>::max() / b) {
        throw std::length_error("mesh dimensions overflow size_t");
    }
    return a * b;
}

}  // namespace

Mesh make_rectangular_mesh(std::size_t nx, std::size_t ny,
                           double width, double height) {
    if (nx == 0 || ny == 0) {
        throw std::invalid_argument("nx and ny must be positive");
    }
    require_positive_finite(width, "width must be finite and positive");
    require_positive_finite(height, "height must be finite and positive");
    const auto largest = std::numeric_limits<std::size_t>::max();
    if (nx == largest || ny == largest) {
        throw std::length_error("mesh dimensions overflow size_t");
    }
    const std::size_t stride = nx + 1;
    const std::size_t vertex_count = checked_product(stride, ny + 1);
    const std::size_t triangle_count = checked_product(checked_product(nx, ny), 2);
    Mesh mesh;
    if (vertex_count > mesh.reference_positions.max_size() ||
        vertex_count > mesh.positions.max_size() ||
        triangle_count > mesh.triangles.max_size()) {
        throw std::length_error("mesh dimensions exceed vector capacity");
    }
    if (width / static_cast<double>(nx) == 0.0 ||
        height / static_cast<double>(ny) == 0.0) {
        throw std::invalid_argument("mesh spacing underflows double precision");
    }
    mesh.reference_positions.reserve(vertex_count);
    mesh.positions.reserve(vertex_count);
    mesh.triangles.reserve(triangle_count);

    for (std::size_t j = 0; j <= ny; ++j) {
        for (std::size_t i = 0; i <= nx; ++i) {
            // Divide first so width * i cannot overflow for finite width.
            const double xi = width * (static_cast<double>(i) / static_cast<double>(nx));
            const double eta = height * (static_cast<double>(j) / static_cast<double>(ny));
            mesh.reference_positions.push_back({xi, eta});
            mesh.positions.push_back({xi, 0.0, eta});
        }
    }
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const std::size_t n0 = j * stride + i;
            const std::size_t n1 = n0 + 1;
            const std::size_t n2 = n0 + stride;
            const std::size_t n3 = n2 + 1;
            mesh.triangles.push_back({n0, n1, n3});
            mesh.triangles.push_back({n0, n3, n2});
        }
    }
    return mesh;
}

std::vector<double> lumped_masses(const Mesh& mesh, double density,
                                  double thickness) {
    require_positive_finite(density, "density must be finite and positive");
    require_positive_finite(thickness, "thickness must be finite and positive");
    if (mesh.reference_positions.empty() || mesh.triangles.empty() ||
        mesh.positions.size() != mesh.reference_positions.size()) {
        throw std::invalid_argument("mesh must contain matching vertices and triangles");
    }
    for (const auto& point : mesh.reference_positions) {
        if (!std::isfinite(point[0]) || !std::isfinite(point[1])) {
            throw std::invalid_argument("reference positions must be finite");
        }
    }
    std::vector<double> masses(mesh.reference_positions.size(), 0.0);
    for (const auto& triangle : mesh.triangles) {
        for (const auto vertex : triangle) {
            if (vertex >= masses.size()) {
                throw std::invalid_argument("triangle vertex index is out of range");
            }
        }
        const auto& a = mesh.reference_positions[triangle[0]];
        const auto& b = mesh.reference_positions[triangle[1]];
        const auto& c = mesh.reference_positions[triangle[2]];
        // Extended intermediates reduce overflow in area * density * thickness.
        const long double bx = static_cast<long double>(b[0]) - a[0];
        const long double by = static_cast<long double>(b[1]) - a[1];
        const long double cx = static_cast<long double>(c[0]) - a[0];
        const long double cy = static_cast<long double>(c[1]) - a[1];
        const long double twice_area = std::abs(bx * cy - by * cx);
        if (!std::isfinite(twice_area) || twice_area <= 0.0L) {
            throw std::invalid_argument("reference triangle must have finite positive area");
        }
        const long double share = (twice_area / 6.0L) * density * thickness;
        const double contribution = static_cast<double>(share);
        if (!std::isfinite(contribution) || contribution <= 0.0) {
            throw std::overflow_error("triangle mass is outside double precision range");
        }
        for (const auto vertex : triangle) {
            const double updated = masses[vertex] + contribution;
            if (!std::isfinite(updated)) {
                throw std::overflow_error("vertex mass overflows double precision");
            }
            masses[vertex] = updated;
        }
    }
    for (const double mass : masses) {
        if (mass <= 0.0) {
            throw std::invalid_argument("mesh contains an isolated vertex with zero mass");
        }
    }
    return masses;
}

}  // namespace cloth
