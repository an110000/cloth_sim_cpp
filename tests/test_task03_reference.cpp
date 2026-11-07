#include "cloth/assembly.hpp"
#include "cloth/dofs.hpp"
#include "cloth/implicit.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

class Reader {
public:
    explicit Reader(const std::string& path) {
        std::ifstream file(path);
        if (!file) throw std::runtime_error("Cannot open task-03 reference data");
        std::string text, line;
        while (std::getline(file, line)) {
            if (!line.empty() && line.front() != '#') text += line + '\n';
        }
        input_.str(text);
    }

    template<class T> T get() {
        T value{};
        if (!(input_ >> value)) throw std::runtime_error("Malformed/truncated task-03 reference");
        return value;
    }

    void expect(const std::string& expected) {
        const auto actual = get<std::string>();
        if (actual != expected)
            throw std::runtime_error("Expected " + expected + ", got " + actual);
    }

    void finish() {
        std::string extra;
        if (input_ >> extra) throw std::runtime_error("Unexpected trailing task-03 reference data");
    }

private:
    std::istringstream input_;
};

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void compare(double actual, double expected, double absolute, const std::string& name) {
    const double tolerance = absolute + 1.e-8 * std::abs(expected);
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(actual - expected) > tolerance) {
        std::ostringstream detail;
        detail.precision(17);
        detail << name << ": actual=" << actual << " expected=" << expected
               << " error=" << std::abs(actual - expected) << " limit=" << tolerance;
        throw std::runtime_error(detail.str());
    }
}

std::vector<cloth::Vec3> read_vectors(Reader& reader, const std::string& tag,
                                    std::size_t count) {
    reader.expect(tag);
    std::vector<cloth::Vec3> result(count);
    for (auto& point : result) for (auto& value : point) value = reader.get<double>();
    return result;
}

void compare_vectors(const std::vector<cloth::Vec3>& actual,
                     const std::vector<cloth::Vec3>& expected,
                     double absolute, const std::string& name) {
    check(actual.size() == expected.size(), name + " size mismatch");
    for (std::size_t i = 0; i < actual.size(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            compare(actual[i][k], expected[i][k], absolute,
                    name + "[" + std::to_string(i) + "][" + std::to_string(k) + "]");
}

double norm(const std::vector<double>& values) {
    double result = 0.0;
    for (const auto value : values) result = std::hypot(result, value);
    return result;
}

struct ResidualCheck {
    double residual;
    double threshold;
};

// Independent from the solver's iteration/control flow: evaluate the returned
// state using the public elastic assembly, mass matrix, and external force.
ResidualCheck final_residual(cloth::Mesh mesh, const cloth::Material& material,
                              const cloth::ImplicitEulerSolver& solver,
                              const cloth::SimulationState& previous,
                              const cloth::SimulationState& current,
                              const cloth::DofMap& map, double dt) {
    mesh.positions = current.positions;
    const auto elastic = cloth::assemble_elastic(mesh, material);
    std::vector<double> delta(map.full_size());
    for (std::size_t i = 0; i < current.positions.size(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            delta[3 * i + k] = current.positions[i][k] -
                (previous.positions[i][k] + dt * previous.velocities[i][k]);
    const auto inertia = solver.mass_matrix().multiply(delta);
    std::vector<double> gravity(delta.size()), elastic_force(delta.size()), residual(delta.size());
    for (std::size_t i = 0; i < delta.size(); ++i) {
        gravity[i] = dt * dt * solver.external_force()[i];
        elastic_force[i] = dt * dt * elastic.gradient[i];
        residual[i] = inertia[i] - gravity[i] + elastic_force[i];
    }
    const double scale = std::max({1.0,
        norm(cloth::restrict_to_free(inertia, map)),
        norm(cloth::restrict_to_free(gravity, map)),
        norm(cloth::restrict_to_free(elastic_force, map))});
    return {norm(cloth::restrict_to_free(residual, map)),
            solver.options().residual_tolerance * scale};
}

} // namespace

int main(int argc, char** argv) {
    try {
        check(argc == 2, "Expected task-03 fixture path");
        Reader reader(argv[1]);
        reader.expect("CLOTH_TASK03_REFERENCE");
        check(reader.get<int>() == 1, "Unsupported task-03 format version");
        reader.expect("cases");
        const auto cases = reader.get<std::size_t>();
        check(cases == 6, "Expected all six task-03 trajectories");
        std::set<std::string> remaining{"rest_lumped", "rest_consistent", "fall_lumped",
            "fall_consistent", "pinned_lumped", "pinned_consistent"};
        std::size_t total_steps = 0;
        double maximum_residual = 0.0;
        for (std::size_t c = 0; c < cases; ++c) {
            reader.expect("case");
            const auto name = reader.get<std::string>();
            check(remaining.erase(name) == 1, "Duplicate/unknown trajectory " + name);
            reader.expect("vertices");
            const auto vertices = reader.get<std::size_t>();
            reader.expect("triangles");
            const auto triangles = reader.get<std::size_t>();
            check(vertices > 0 && vertices <= 100 && triangles > 0 && triangles <= 200,
                  "Invalid reference mesh size");
            cloth::Mesh mesh;
            mesh.reference_positions.resize(vertices);
            mesh.triangles.resize(triangles);
            reader.expect("reference");
            for (auto& point : mesh.reference_positions)
                for (auto& value : point) value = reader.get<double>();
            reader.expect("topology");
            for (auto& triangle : mesh.triangles)
                for (auto& vertex : triangle) vertex = reader.get<std::size_t>();
            reader.expect("material");
            const cloth::Material material{reader.get<double>(), reader.get<double>(),
                reader.get<double>(), reader.get<double>()};
            reader.expect("density");
            const auto density = reader.get<double>();
            reader.expect("mass_mode");
            const auto mode_name = reader.get<std::string>();
            check(mode_name == "lumped" || mode_name == "consistent", "Unknown mass mode");
            const auto mode = mode_name == "lumped" ? cloth::MassMode::lumped : cloth::MassMode::consistent;
            reader.expect("gravity");
            const cloth::Vec3 gravity{reader.get<double>(), reader.get<double>(), reader.get<double>()};
            reader.expect("dt");
            const auto dt = reader.get<double>();
            reader.expect("newton");
            cloth::NewtonOptions options;
            options.residual_tolerance = reader.get<double>();
            options.max_iterations = reader.get<std::size_t>();
            reader.expect("fixed");
            const auto fixed_count = reader.get<std::size_t>();
            check(fixed_count <= vertices, "Invalid fixed vertex count");
            std::vector<std::size_t> fixed(fixed_count);
            for (auto& vertex : fixed) vertex = reader.get<std::size_t>();
            cloth::SimulationState state{read_vectors(reader, "positions", vertices),
                                          read_vectors(reader, "velocities", vertices)};
            mesh.positions = state.positions;
            const cloth::ImplicitEulerSolver solver(mesh, material, density, mode, gravity, options);
            const cloth::DofMap map(vertices, fixed);
            reader.expect("steps");
            const auto steps = reader.get<std::size_t>();
            check(steps == 8, "Expected all eight steps in " + name);
            for (std::size_t index = 1; index <= steps; ++index) {
                reader.expect("step");
                check(reader.get<std::size_t>() == index, "Invalid reference step order");
                const auto expected_x = read_vectors(reader, "positions", vertices);
                const auto expected_v = read_vectors(reader, "velocities", vertices);
                reader.expect("residual_norm");
                const auto reference_residual = reader.get<double>();
                reader.expect("residual_threshold");
                const auto reference_threshold = reader.get<double>();
                reader.expect("iterations");
                const auto reference_iterations = reader.get<std::size_t>();
                reader.expect("end_step");
                check(std::isfinite(reference_residual) && reference_residual >= 0.0 &&
                      std::isfinite(reference_threshold) && reference_threshold > 0.0 &&
                      reference_residual <= reference_threshold && reference_iterations > 0 &&
                      reference_iterations <= options.max_iterations,
                      "Invalid unconverged Python reference in " + name);
                const auto previous = state;
                const auto result = solver.step(state, fixed, dt);
                const auto label = name + " step " + std::to_string(index);
                check(result.info.committed && result.info.solver_converged &&
                      result.info.finite_state && result.info.failure_reason == cloth::StepFailure::none,
                      label + " failed: " + cloth::failure_name(result.info.failure_reason) +
                      " " + result.info.message);
                check(state.positions == previous.positions && state.velocities == previous.velocities,
                      label + " mutated input state");
                compare_vectors(result.state.positions, expected_x, 2.e-9, label + " positions");
                compare_vectors(result.state.velocities, expected_v, 2.e-7, label + " velocities");
                for (const auto vertex : fixed) {
                    check(result.state.positions[vertex] == previous.positions[vertex],
                          label + " changed fixed position");
                    check(result.state.velocities[vertex] == cloth::Vec3{0.0, 0.0, 0.0},
                          label + " changed fixed velocity");
                }
                const auto evaluated = final_residual(mesh, material, solver, previous,
                                                     result.state, map, dt);
                check(std::isfinite(evaluated.residual) &&
                      evaluated.residual <= evaluated.threshold * (1.0 + 1.e-8),
                      label + " independent final residual fails convergence");
                compare(result.info.residual_norm, evaluated.residual, 1.e-14,
                        label + " reported final residual");
                compare(result.info.residual_threshold, evaluated.threshold, 1.e-15,
                        label + " reported threshold");
                compare(evaluated.threshold, reference_threshold, 1.e-15,
                        label + " Python threshold");
                maximum_residual = std::max(maximum_residual, evaluated.residual);
                state = result.state;
                ++total_steps;
            }
            reader.expect("end_case");
        }
        reader.expect("end");
        reader.finish();
        check(remaining.empty() && total_steps == 48, "Incomplete trajectory coverage");
        std::cout << "6 Python trajectories / 48 committed steps passed; max final residual="
                  << maximum_residual << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
