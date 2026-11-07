#include "cloth/scene.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const char* message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance)
        throw std::runtime_error(message);
}

template<class Exception, class Function>
void rejects(Function function, const char* message) {
    try { function(); } catch (const Exception&) { return; }
    throw std::runtime_error(message);
}

std::string hash(char digit) {
    return std::string(64, digit);
}

void default_scene_matches_teaching_cli() {
    const auto config = cloth::make_default_scene_config();
    const auto mesh = cloth::make_rectangular_mesh(2, 1, 1.0, .7);
    check(config.rectangle().nx == 2 && config.rectangle().ny == 1,
          "Default rectangle cell counts");
    check(config.mesh().triangles == mesh.triangles &&
          config.mesh().reference_positions == mesh.reference_positions,
          "Default mesh topology and reference coordinates");
    check(config.fixed_vertices() == std::vector<std::size_t>{3, 4, 5},
          "Default pinned vertices");
    check(config.solver() == cloth::SceneSolver::xpbd_distance,
          "Default solver is distance XPBD");
    check(config.mass_mode() == cloth::MassMode::lumped && config.density() == 800.0,
          "Default physical mass parameters");
    check(config.gravity() == cloth::Vec3{0.0, -9.81, 0.0} &&
          config.dt() == .01 && config.steps() == 8 && config.run_seed() == 0,
          "Default integration parameters");
    check(config.contact().enabled == false && config.asset_sha256().empty() &&
          config.scene_sha256().empty() && config.state_sha256().empty(),
          "Default contact and provenance fields");
    const auto expected_masses = cloth::lumped_masses(mesh, 800.0, .001);
    check(config.masses() == expected_masses, "Default masses use physical thickness and density");
    check(config.initial_state().positions.size() == 6 && config.initial_state().velocities.size() == 6,
          "Default state has six vertices");
    near(config.initial_state().positions[0][0], .025, 1e-15, "Default initial x0 matches CLI transform");
    near(config.initial_state().positions[0][1], -.025, 1e-15, "Default initial y0 matches CLI transform");
    near(config.initial_state().positions[0][2], -.025, 1e-15, "Default initial z0 matches CLI transform");
    near(config.initial_state().positions[1][0], .555, 1e-15, "Default second x matches CLI transform");
    near(config.initial_state().positions[1][1], -.0375, 1e-15, "Default second y matches CLI transform");
    near(config.initial_state().positions[1][2], -.025, 1e-15, "Default second z matches CLI transform");
    check(config.initial_state().velocities[0] == cloth::Vec3{.015, -.02, .01},
          "Default initial velocity matches CLI");
    check(config.initial_state().velocities[3] == cloth::Vec3{0, 0, 0},
          "Default fixed velocity starts at zero");
}

void canonical_and_immutable_value_semantics() {
    auto input = cloth::default_scene_config_input();
    input.fixed_vertices = {5, 3, 4};
    input.asset_sha256 = hash('a');
    input.scene_sha256 = hash('b');
    input.state_sha256 = hash('c');
    input.steps = 0;
    input.run_seed = std::numeric_limits<std::uint64_t>::max();
    input.initial_positions[0] = {10, 20, 30};
    const auto config = cloth::SceneConfig(input);
    check(config.fixed_vertices() == std::vector<std::size_t>{3, 4, 5},
          "Fixed vertices have deterministic ascending order");
    check(config.asset_sha256() == hash('a') && config.scene_sha256() == hash('b') &&
          config.state_sha256() == hash('c'), "Provenance hashes are retained");
    check(config.steps() == 0 && config.run_seed() == std::numeric_limits<std::uint64_t>::max(),
          "Zero-step construction and seed are retained");
    input.initial_positions[0] = {-100, -100, -100};
    input.fixed_vertices.clear();
    input.asset_sha256.clear();
    check(config.initial_state().positions[0] == cloth::Vec3{10, 20, 30} &&
          config.fixed_vertices() == std::vector<std::size_t>{3, 4, 5} &&
          config.asset_sha256() == hash('a'), "Scene owns immutable copies of input data");
    check(cloth::valid_sha256("") && cloth::valid_sha256(hash('0')) &&
          !cloth::valid_sha256(std::string(63, '0')) && !cloth::valid_sha256(std::string(64, 'A')),
          "SHA-256 format predicate");
    const auto sorted = cloth::canonical_fixed_vertices({8, 2, 5}, 9);
    check(sorted == std::vector<std::size_t>{2, 5, 8}, "Canonical helper sorts stable set");
}

void custom_state_and_solver_options() {
    cloth::SceneConfigInput input;
    input.rectangle = {1, 2, 2.0, 3.0};
    const auto count = (input.rectangle.nx + 1) * (input.rectangle.ny + 1);
    input.initial_positions.assign(count, cloth::Vec3{1, 2, 3});
    input.initial_velocities.assign(count, cloth::Vec3{.1, .2, .3});
    input.fixed_vertices = {0, 3};
    input.gravity = {1, 2, 3};
    input.dt = .02;
    input.steps = 17;
    input.solver = cloth::SceneSolver::xpbd_energy;
    input.mass_mode = cloth::MassMode::consistent;
    input.contact.enabled = true;
    input.contact.restitution = .4;
    input.contact.stabilization_beta = .2;
    input.xpbd_options.iterations = 7;
    const auto config = cloth::SceneConfig(input);
    check(config.mesh().positions.size() == count && config.initial_state().positions == input.initial_positions,
          "Custom positions are copied without rectangle overwrite");
    check(config.initial_state().velocities == input.initial_velocities && config.fixed_vertices() == std::vector<std::size_t>{0, 3},
          "Custom velocities and fixed set");
    check(config.solver() == cloth::SceneSolver::xpbd_energy && config.mass_mode() == cloth::MassMode::consistent &&
          config.xpbd_options().iterations == 7 && config.contact().restitution == .4,
          "Custom solver/contact options");
    const auto default_velocity = cloth::SceneConfigInput{};
    check(default_velocity.initial_positions.empty() && default_velocity.initial_velocities.empty(),
          "Input distinguishes omitted initial arrays from generated state");
}

void arbitrary_topology_is_owned_and_validated() {
    cloth::SceneConfigInput input;
    cloth::Mesh mesh;
    mesh.reference_positions = {{0,0},{1,0},{0,1},{1,1}};
    mesh.positions = {{0,0,0},{1,0,0},{0,0,1},{1,0,1}};
    mesh.triangles = {{0,1,2},{1,3,2}};
    input.custom_mesh = mesh;
    input.initial_positions = mesh.positions;
    input.initial_positions[3][1] = -.01;
    input.fixed_vertices = {0,1};
    input.steps = 1;
    const auto config = cloth::SceneConfig(input);
    check(config.has_custom_mesh() && config.mesh().triangles == mesh.triangles,
          "Arbitrary topology is preserved");
    check(config.initial_state().positions[3][1] == -.01 && config.masses().size() == 4,
          "Arbitrary topology owns initial state and masses");
    input.custom_mesh->positions[3][0] = 99;
    check(config.mesh().positions[3][0] == 1, "Custom mesh is copied immutably");
    auto bad = input;
    bad.custom_mesh->triangles = {{0,1,4}};
    rejects<std::invalid_argument>([&] { cloth::SceneConfig{bad}; }, "Out-of-range custom triangle accepted");
    bad = input;
    bad.custom_mesh->reference_positions[2] = {0,0};
    rejects<std::invalid_argument>([&] { cloth::SceneConfig{bad}; }, "Degenerate custom reference triangle accepted");
}

void invalid_configuration() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double infinity = std::numeric_limits<double>::infinity();
    rejects<std::invalid_argument>([] { cloth::SceneConfigInput i; i.rectangle.nx = 0; cloth::validate_scene_config(i); }, "Zero nx accepted");
    rejects<std::invalid_argument>([] { cloth::SceneConfigInput i; i.rectangle.height = -1; (void)cloth::SceneConfig{i}; }, "Negative rectangle size accepted");
    rejects<std::invalid_argument>([&] { cloth::SceneConfigInput i; i.rectangle.width = nan; (void)cloth::SceneConfig{i}; }, "Non-finite rectangle size accepted");
    rejects<std::invalid_argument>([] { cloth::SceneConfigInput i; i.fixed_vertices = {3, 3}; (void)cloth::SceneConfig{i}; }, "Duplicate fixed vertex accepted");
    rejects<std::invalid_argument>([] { cloth::SceneConfigInput i; i.fixed_vertices = {6}; (void)cloth::SceneConfig{i}; }, "Out-of-range fixed vertex accepted");
    rejects<std::invalid_argument>([] { cloth::SceneConfigInput i; i.initial_positions.resize(1); (void)cloth::SceneConfig{i}; }, "Wrong initial position size accepted");
    rejects<std::invalid_argument>([] { cloth::SceneConfigInput i; i.initial_velocities.resize(1); (void)cloth::SceneConfig{i}; }, "Wrong initial velocity size accepted");
    rejects<std::invalid_argument>([&] { cloth::SceneConfigInput i; i.gravity[1] = infinity; (void)cloth::SceneConfig{i}; }, "Non-finite gravity accepted");
    rejects<std::invalid_argument>([] { cloth::SceneConfigInput i; i.dt = 0; (void)cloth::SceneConfig{i}; }, "Zero dt accepted");
    rejects<std::invalid_argument>([] { cloth::SceneConfigInput i; i.dt = std::numeric_limits<double>::denorm_min(); (void)cloth::SceneConfig{i}; }, "Underflowing dt accepted");
    rejects<std::invalid_argument>([] { cloth::SceneConfigInput i; i.density = 0; (void)cloth::SceneConfig{i}; }, "Zero density accepted");
    rejects<std::invalid_argument>([] { cloth::SceneConfigInput i; i.material.poisson = .5; (void)cloth::SceneConfig{i}; }, "Inadmissible Poisson ratio accepted");
    rejects<std::invalid_argument>([] { cloth::SceneConfigInput i; i.material.thickness = -1; (void)cloth::SceneConfig{i}; }, "Negative material thickness accepted");
    rejects<std::invalid_argument>([] { cloth::SceneConfigInput i; i.xpbd_options.iterations = 0; (void)cloth::SceneConfig{i}; }, "Zero XPBD iterations accepted");
    rejects<std::invalid_argument>([] { cloth::SceneConfigInput i; i.newton_options.backtrack_factor = 1; (void)cloth::SceneConfig{i}; }, "Unit Newton backtrack accepted");
    rejects<std::invalid_argument>([] { cloth::SceneConfigInput i; i.newton_options.diagonal_shifts = {0, 0}; (void)cloth::SceneConfig{i}; }, "Duplicate Newton shift accepted");
    rejects<std::invalid_argument>([] { cloth::SceneConfigInput i; i.asset_sha256 = std::string(64, 'A'); (void)cloth::SceneConfig{i}; }, "Uppercase hash accepted");
    rejects<std::invalid_argument>([] { cloth::SceneConfigInput i; i.scene_sha256 = "not-a-hash"; (void)cloth::SceneConfig{i}; }, "Malformed scene hash accepted");
    rejects<std::invalid_argument>([] { cloth::SceneConfigInput i; i.solver = cloth::SceneSolver::implicit_euler; i.contact.enabled = true; i.contact.restitution = .1; (void)cloth::SceneConfig{i}; }, "Implicit restitution silently accepted");
    rejects<std::invalid_argument>([] { cloth::SceneConfigInput i; i.solver = static_cast<cloth::SceneSolver>(99); (void)cloth::SceneConfig{i}; }, "Unknown solver accepted");
    rejects<std::invalid_argument>([] { (void)cloth::canonical_fixed_vertices({1, 1}, 2); }, "Canonical helper hides duplicate");
    rejects<std::invalid_argument>([] { (void)cloth::canonical_fixed_vertices({2}, 2); }, "Canonical helper accepts out-of-range");
    const auto valid = cloth::SceneConfigInput{};
    cloth::validate_scene_config(valid);
}

} // namespace

int main() {
    try {
        default_scene_matches_teaching_cli();
        canonical_and_immutable_value_semantics();
        custom_state_and_solver_options();
        arbitrary_topology_is_owned_and_validated();
        invalid_configuration();
        std::cout << "scene configuration, provenance, canonical ordering and validation passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "scene test failed: " << error.what() << '\n';
        return 1;
    }
}
