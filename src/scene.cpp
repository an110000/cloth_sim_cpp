#include "cloth/scene.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace cloth {
namespace {

bool finite_vec3(const Vec3& value) {
    return std::all_of(value.begin(), value.end(), [](double x) { return std::isfinite(x); });
}

bool finite_vec2(const Vec2& value) {
    return std::all_of(value.begin(), value.end(), [](double x) { return std::isfinite(x); });
}

void validate_mesh(const Mesh& mesh) {
    if (mesh.reference_positions.empty() || mesh.positions.size() != mesh.reference_positions.size() ||
        mesh.triangles.empty()) {
        throw std::invalid_argument("Custom scene mesh must contain vertices and triangles");
    }
    for (const auto& point : mesh.reference_positions) {
        if (!finite_vec2(point)) throw std::invalid_argument("Custom reference positions must be finite");
    }
    for (const auto& point : mesh.positions) {
        if (!finite_vec3(point)) throw std::invalid_argument("Custom rest positions must be finite");
    }
    for (const auto& triangle : mesh.triangles) {
        if (triangle[0] >= mesh.positions.size() || triangle[1] >= mesh.positions.size() ||
            triangle[2] >= mesh.positions.size() || triangle[0] == triangle[1] ||
            triangle[1] == triangle[2] || triangle[0] == triangle[2]) {
            throw std::invalid_argument("Custom mesh triangle index is invalid");
        }
        const auto& a = mesh.reference_positions[triangle[0]];
        const auto& b = mesh.reference_positions[triangle[1]];
        const auto& c = mesh.reference_positions[triangle[2]];
        const double area2 = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]);
        if (!std::isfinite(area2) || std::abs(area2) <= 1.0e-14)
            throw std::invalid_argument("Custom mesh contains a degenerate reference triangle");
    }
}

std::size_t vertex_count(const RectangularSceneSpec& rectangle) {
    const auto maximum = std::numeric_limits<std::size_t>::max();
    if (rectangle.nx == maximum || rectangle.ny == maximum ||
        rectangle.nx + 1 > maximum / (rectangle.ny + 1)) {
        throw std::length_error("Scene rectangle vertex count overflows size_t");
    }
    const auto count = (rectangle.nx + 1) * (rectangle.ny + 1);
    if (count == 0) throw std::length_error("Scene rectangle vertex count is zero");
    return count;
}

void validate_material(const Material& material) {
    if (!std::isfinite(material.young) || material.young <= 0.0 ||
        !std::isfinite(material.poisson) || material.poisson <= -1.0 || material.poisson >= 0.5 ||
        !std::isfinite(material.thickness) || material.thickness <= 0.0 ||
        !std::isfinite(material.j_min) || material.j_min < 0.0) {
        throw std::invalid_argument("Scene material is outside its finite physical domain");
    }
}

void validate_newton(const NewtonOptions& options) {
    if (!std::isfinite(options.residual_tolerance) || options.residual_tolerance <= 0.0 ||
        options.line_search_max_iterations == 0 ||
        !std::isfinite(options.backtrack_factor) || options.backtrack_factor <= 0.0 || options.backtrack_factor >= 1.0 ||
        !std::isfinite(options.armijo) || options.armijo <= 0.0 || options.armijo >= 1.0 ||
        !std::isfinite(options.linear_relative_tolerance) || options.linear_relative_tolerance <= 0.0 ||
        options.diagonal_shifts.empty()) {
        throw std::invalid_argument("Scene Newton options are invalid");
    }
    double previous = -1.0;
    for (const double value : options.diagonal_shifts) {
        if (!std::isfinite(value) || value < 0.0 || value <= previous) {
            throw std::invalid_argument("Scene Newton diagonal shifts must be increasing and finite");
        }
        previous = value;
    }
}

void validate_xpbd(const XpbdOptions& options) {
    if (options.iterations == 0 ||
        !std::isfinite(options.distance_tolerance) || options.distance_tolerance <= 0.0 ||
        !std::isfinite(options.angle_tolerance) || options.angle_tolerance <= 0.0 ||
        !std::isfinite(options.energy_tolerance) || options.energy_tolerance <= 0.0) {
        throw std::invalid_argument("Scene XPBD options are invalid");
    }
}

void validate_hash(const std::string& value, const char* name) {
    if (value.empty()) return;
    if (value.size() != 64 || !std::all_of(value.begin(), value.end(), [](unsigned char c) {
            return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
        })) {
        throw std::invalid_argument(std::string(name) + " must be lowercase SHA-256 hex or empty");
    }
}

} // namespace

bool valid_sha256(const std::string& value) noexcept {
    return value.empty() || (value.size() == 64 && std::all_of(value.begin(), value.end(), [](unsigned char c) {
        return (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f');
    }));
}

std::vector<std::size_t> canonical_fixed_vertices(
    std::vector<std::size_t> vertices, std::size_t count) {
    for (const auto vertex : vertices) {
        if (vertex >= count) throw std::invalid_argument("Scene fixed vertex is out of range");
    }
    std::sort(vertices.begin(), vertices.end());
    if (std::adjacent_find(vertices.begin(), vertices.end()) != vertices.end()) {
        throw std::invalid_argument("Scene fixed vertices contain a duplicate");
    }
    return vertices;
}

void validate_scene_config(const SceneConfigInput& input) {
    const std::size_t count = input.custom_mesh ? input.custom_mesh->positions.size() : vertex_count(input.rectangle);
    if (!input.custom_mesh) {
        if (input.rectangle.nx == 0 || input.rectangle.ny == 0 ||
            !std::isfinite(input.rectangle.width) || input.rectangle.width <= 0.0 ||
            !std::isfinite(input.rectangle.height) || input.rectangle.height <= 0.0) {
            throw std::invalid_argument("Scene rectangle dimensions must be positive and finite");
        }
        if (input.rectangle.width / static_cast<double>(input.rectangle.nx) == 0.0 ||
            input.rectangle.height / static_cast<double>(input.rectangle.ny) == 0.0) {
            throw std::invalid_argument("Scene rectangle spacing underflows double precision");
        }
    } else {
        validate_mesh(*input.custom_mesh);
    }
    if ((!input.initial_positions.empty() && input.initial_positions.size() != count) ||
        (!input.initial_velocities.empty() && input.initial_velocities.size() != count)) {
        throw std::invalid_argument("Scene initial arrays must match rectangle vertex count");
    }
    for (const auto& point : input.initial_positions) {
        if (!finite_vec3(point)) throw std::invalid_argument("Scene initial positions must be finite");
    }
    for (const auto& velocity : input.initial_velocities) {
        if (!finite_vec3(velocity)) throw std::invalid_argument("Scene initial velocities must be finite");
    }
    for (const auto value : input.gravity) {
        if (!std::isfinite(value)) throw std::invalid_argument("Scene gravity must be finite");
    }
    if (!std::isfinite(input.dt) || input.dt <= 0.0 || !std::isfinite(input.dt * input.dt) || input.dt * input.dt == 0.0) {
        throw std::invalid_argument("Scene dt must be positive and have a finite nonzero square");
    }
    if (static_cast<long double>(input.steps) * static_cast<long double>(input.dt) >
        static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::invalid_argument("Scene requested duration exceeds finite double range");
    }
    if (!std::isfinite(input.density) || input.density <= 0.0) {
        throw std::invalid_argument("Scene density must be positive and finite");
    }
    validate_material(input.material);
    validate_ground_contact_options(input.contact);
    validate_newton(input.newton_options);
    validate_xpbd(input.xpbd_options);
    if (input.mass_mode != MassMode::lumped && input.mass_mode != MassMode::consistent) {
        throw std::invalid_argument("Unknown scene mass mode");
    }
    if (input.solver != SceneSolver::implicit_euler &&
        input.solver != SceneSolver::xpbd_distance && input.solver != SceneSolver::xpbd_energy) {
        throw std::invalid_argument("Unknown scene solver selection");
    }
    if (input.solver == SceneSolver::implicit_euler && input.contact.enabled &&
        (input.contact.restitution != 0.0 || input.contact.stabilization_beta != 0.0)) {
        throw std::invalid_argument("Implicit scene contact supports penalty friction only; restitution/beta are XPBD options");
    }
    (void)canonical_fixed_vertices(input.fixed_vertices, count);
    validate_hash(input.asset_sha256, "asset_sha256");
    validate_hash(input.scene_sha256, "scene_sha256");
    validate_hash(input.state_sha256, "state_sha256");
}

SceneConfig::SceneConfig(SceneConfigInput input)
    : rectangle_(input.rectangle), gravity_(input.gravity), dt_(input.dt), steps_(input.steps),
      run_seed_(input.run_seed), material_(input.material), density_(input.density),
      mass_mode_(input.mass_mode), contact_(input.contact), solver_(input.solver),
      newton_options_(input.newton_options), xpbd_options_(input.xpbd_options),
      asset_sha256_(input.asset_sha256), scene_sha256_(input.scene_sha256),
      state_sha256_(input.state_sha256) {
    validate_scene_config(input);
    custom_mesh_ = input.custom_mesh.has_value();
    mesh_ = input.custom_mesh ? std::move(*input.custom_mesh)
                              : make_rectangular_mesh(rectangle_.nx, rectangle_.ny, rectangle_.width, rectangle_.height);
    fixed_vertices_ = canonical_fixed_vertices(std::move(input.fixed_vertices), mesh_.positions.size());
    initial_state_.positions = input.initial_positions.empty() ? mesh_.positions : std::move(input.initial_positions);
    initial_state_.velocities = input.initial_velocities.empty()
        ? std::vector<Vec3>(mesh_.positions.size(), Vec3{0.0, 0.0, 0.0})
        : std::move(input.initial_velocities);
    masses_ = lumped_masses(mesh_, density_, material_.thickness);
}

SceneConfigInput default_scene_config_input() {
    SceneConfigInput input;
    input.fixed_vertices = {3, 4, 5};
    input.solver = SceneSolver::xpbd_distance;
    const auto mesh = make_rectangular_mesh(input.rectangle.nx, input.rectangle.ny,
                                            input.rectangle.width, input.rectangle.height);
    input.initial_positions = mesh.positions;
    input.initial_velocities.assign(mesh.positions.size(), Vec3{0.0, 0.0, 0.0});
    for (std::size_t vertex = 0; vertex < 3; ++vertex) {
        input.initial_positions[vertex][0] = 1.06 * input.initial_positions[vertex][0] + 0.025;
        input.initial_positions[vertex][1] = -0.025 * (1.0 + mesh.reference_positions[vertex][0]);
        input.initial_positions[vertex][2] -= 0.025;
        input.initial_velocities[vertex] = {0.015, -0.02, 0.01};
    }
    input.material = {2.0e4, 0.3, 0.001, 0.05};
    input.density = 800.0;
    return input;
}

SceneConfig make_default_scene_config() {
    return SceneConfig(default_scene_config_input());
}

} // namespace cloth
