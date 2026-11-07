#include "cloth/constraints.hpp"
#include "cloth/xpbd.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

class Reader {
public:
    explicit Reader(const std::string& path) {
        std::ifstream file(path);
        check(static_cast<bool>(file), "Cannot open task-05A reference");
        std::string line, all;
        while (std::getline(file, line))
            if (!line.empty() && line.front() != '#') all += line + '\n';
        stream_.str(all);
    }
    template<class T> T get() {
        T value{};
        if (!(stream_ >> value)) throw std::runtime_error("Truncated or invalid task-05A reference");
        return value;
    }
    void expect(const std::string& expected) {
        const auto actual = get<std::string>();
        check(actual == expected, "Expected " + expected + ", got " + actual);
    }
    void finish() {
        std::string extra;
        check(!(stream_ >> extra), "Trailing task-05A reference data");
    }
private:
    std::istringstream stream_;
};

void near(double actual, double expected, double absolute, const std::string& label) {
    const double tolerance = absolute + 5.e-9*std::abs(expected);
    if (!std::isfinite(actual) || !std::isfinite(expected) || std::abs(actual-expected) > tolerance) {
        std::ostringstream message;
        message.precision(17);
        message << label << ": actual=" << actual << " expected=" << expected
                << " error=" << std::abs(actual-expected) << " tolerance=" << tolerance;
        throw std::runtime_error(message.str());
    }
}

std::vector<double> scalars(Reader& reader, const std::string& tag, std::size_t count) {
    reader.expect(tag);
    std::vector<double> result(count);
    for (auto& value : result) value = reader.get<double>();
    return result;
}

std::vector<cloth::Vec3> points(Reader& reader, const std::string& tag, std::size_t count) {
    reader.expect(tag);
    std::vector<cloth::Vec3> result(count);
    for (auto& point : result) for (auto& value : point) value = reader.get<double>();
    return result;
}

void compare_points(const std::vector<cloth::Vec3>& actual,
                    const std::vector<cloth::Vec3>& expected, double absolute,
                    const std::string& label) {
    check(actual.size() == expected.size(), label + " count");
    for (std::size_t node = 0; node < actual.size(); ++node)
        for (std::size_t axis = 0; axis < 3; ++axis)
            near(actual[node][axis], expected[node][axis], absolute,
                 label + "[" + std::to_string(node) + "][" + std::to_string(axis) + "]");
}

cloth::Material material(Reader& reader) {
    const auto values = scalars(reader, "material", 4);
    return {values[0], values[1], values[2], values[3]};
}

double scalar(Reader& reader, const std::string& tag) {
    reader.expect(tag);
    return reader.get<double>();
}

void evaluation(const cloth::Constraint& constraint, const std::vector<cloth::Vec3>& x,
                double value, const std::vector<cloth::Vec3>& gradient, const std::string& label) {
    const auto result = cloth::evaluate_constraint(constraint, x);
    const auto& energy = std::get<cloth::NeoHookeanEnergyConstraint>(constraint);
    check(result.valid_geometry && result.count == 3, label + " valid energy geometry");
    check(cloth::constraint_group(constraint) == cloth::ConstraintGroup::neohookean_energy,
          label + " energy group");
    near(result.value, value, 2.e-11, label + " C=E/scale");
    for (std::size_t local = 0; local < 3; ++local) {
        check(result.nodes[local] == energy.nodes[local], label + " local node ordering");
        for (std::size_t axis = 0; axis < 3; ++axis)
            near(result.gradients[local][axis], gradient[local][axis], 2.e-10,
                 label + " scaled gradient");
    }
}

void locals(Reader& reader) {
    reader.expect("locals");
    check(reader.get<std::size_t>() == 6, "Expected six energy local cases");
    std::set<std::string> names{"rest", "nonlinear", "scaled_constraint", "nonzero_multiplier",
                                "hard_energy", "log_continuation"};
    for (std::size_t index = 0; index < 6; ++index) {
        reader.expect("local");
        const auto name = reader.get<std::string>();
        check(names.erase(name) == 1, "Unknown or duplicate local case");
        cloth::NeoHookeanEnergyConstraint energy;
        reader.expect("nodes");
        for (auto& node : energy.nodes) node = reader.get<std::size_t>();
        reader.expect("reference");
        for (auto& point : energy.reference) for (auto& value : point) value = reader.get<double>();
        energy.material = material(reader);
        energy.compliance = scalar(reader, "compliance");
        energy.energy_scale = scalar(reader, "energy_scale");
        const auto dt = scalar(reader, "dt");
        energy.lambda_value = scalar(reader, "initial_lambda");
        cloth::Constraint constraint = energy;
        const auto inverse = scalars(reader, "inverse_masses", 3);
        auto x = points(reader, "positions", 3);
        const auto value = scalar(reader, "value");
        const auto gradient = points(reader, "gradient", 3);
        evaluation(constraint, x, value, gradient, name);
        reader.expect("projections");
        check(reader.get<std::size_t>() == 2, "Expected two accumulated projections");
        for (std::size_t projection = 1; projection <= 2; ++projection) {
            reader.expect("projection");
            check(reader.get<std::size_t>() == projection, "Invalid local projection order");
            const auto delta = scalar(reader, "delta_lambda");
            const auto before = scalar(reader, "value_before");
            const auto multiplier = scalar(reader, "lambda");
            const auto expected_x = points(reader, "positions", 3);
            const auto after = scalar(reader, "value_after");
            const auto after_gradient = points(reader, "gradient_after", 3);
            reader.expect("end_projection");
            const auto previous = x;
            const auto result = cloth::project_constraint(constraint, x, inverse, dt);
            check(result.valid_geometry, name + " projection geometry");
            near(result.delta_lambda, delta, 2.e-12, name + " delta lambda");
            near(result.value_before, before, 2.e-11, name + " pre-projection value");
            near(cloth::constraint_multiplier(constraint), multiplier, 2.e-12, name + " accumulated lambda");
            compare_points(x, expected_x, 2.e-10, name + " projection positions");
            for (std::size_t node = 0; node < 3; ++node)
                if (inverse[node] == 0.) check(x[node] == previous[node], name + " fixed node drift");
            evaluation(constraint, x, after, after_gradient, name + " projected");
        }
        const auto reset = scalar(reader, "reset_lambda");
        cloth::reset_constraint(constraint);
        near(cloth::constraint_multiplier(constraint), reset, 0., name + " reset");
        reader.expect("end_local");
    }
    check(names.empty(), "Missing local case");
}

cloth::Mesh mesh(Reader& reader) {
    reader.expect("vertices");
    const auto vertices = reader.get<std::size_t>();
    reader.expect("triangles");
    const auto triangles = reader.get<std::size_t>();
    check(vertices == 6 && triangles == 4, "Expected six-node/four-element reference mesh");
    cloth::Mesh result;
    result.reference_positions.resize(vertices);
    result.triangles.resize(triangles);
    reader.expect("reference");
    for (auto& point : result.reference_positions) for (auto& value : point) value = reader.get<double>();
    result.positions = points(reader, "rest_positions", vertices);
    reader.expect("topology");
    for (auto& triangle : result.triangles) for (auto& node : triangle) node = reader.get<std::size_t>();
    return result;
}

void trajectories(Reader& reader) {
    reader.expect("trajectories");
    check(reader.get<std::size_t>() == 4, "Expected four energy trajectories");
    std::set<std::string> names{"rest", "freefall", "pinned_deformed", "pinned_scaled_equivalent"};
    std::size_t count_steps = 0;
    for (std::size_t trajectory = 0; trajectory < 4; ++trajectory) {
        reader.expect("trajectory");
        const auto name = reader.get<std::string>();
        check(names.erase(name) == 1, "Unknown or duplicate energy trajectory");
        auto geometry = mesh(reader);
        const auto mat = material(reader);
        const auto density = scalar(reader, "density");
        const auto masses = scalars(reader, "physical_masses", geometry.positions.size());
        const auto assembled = cloth::lumped_masses(geometry, density, mat.thickness);
        for (std::size_t node = 0; node < masses.size(); ++node)
            near(assembled[node], masses[node], 1.e-13, name + " physical mass");
        const auto compliance = scalar(reader, "compliance");
        const auto scale = scalar(reader, "energy_scale");
        reader.expect("iterations");
        cloth::XpbdOptions options;
        options.iterations = reader.get<std::size_t>();
        check(options.iterations > 0 && options.iterations <= 100, "Invalid reference sweep count");
        const auto gravity_data = scalars(reader, "gravity", 3);
        const cloth::Vec3 gravity{gravity_data[0], gravity_data[1], gravity_data[2]};
        const auto dt = scalar(reader, "dt");
        reader.expect("fixed");
        const auto fixed_count = reader.get<std::size_t>();
        check(fixed_count <= masses.size(), "Invalid reference fixed count");
        std::vector<std::size_t> fixed(fixed_count);
        for (auto& node : fixed) node = reader.get<std::size_t>();
        cloth::SimulationState state{points(reader, "positions", masses.size()),
                                     points(reader, "velocities", masses.size())};
        const auto constraints = cloth::build_neohookean_constraints(geometry, mat, compliance, scale);
        check(constraints.size() == geometry.triangles.size(), name + " one constraint per element");
        for (std::size_t element = 0; element < constraints.size(); ++element) {
            const auto& c = std::get<cloth::NeoHookeanEnergyConstraint>(constraints[element]);
            check(c.nodes == geometry.triangles[element], name + " ordered element nodes");
            for (std::size_t local = 0; local < 3; ++local)
                check(c.reference[local] == geometry.reference_positions[c.nodes[local]], name + " reference copy");
            near(c.energy_scale, scale, 0., name + " scale");
            near(c.compliance, compliance, 0., name + " compliance");
        }
        const cloth::XpbdSolver solver(masses, constraints, gravity, options);
        reader.expect("steps");
        check(reader.get<std::size_t>() == 8, "Expected eight physical steps");
        for (std::size_t step = 1; step <= 8; ++step) {
            reader.expect("step");
            check(reader.get<std::size_t>() == step, "Invalid trajectory step order");
            const auto expected_x = points(reader, "positions", masses.size());
            const auto expected_v = points(reader, "velocities", masses.size());
            const auto multipliers = scalars(reader, "multipliers", constraints.size());
            const auto raw = scalar(reader, "raw_error");
            const auto compliant = scalar(reader, "compliant_error");
            reader.expect("end_step");
            const auto previous = state;
            const auto result = solver.step(state, fixed, dt);
            const auto label = name + " step " + std::to_string(step);
            check(result.info.committed && result.info.finite_state &&
                  result.info.failure_reason == cloth::XpbdFailure::none, label + " failed");
            check(state.positions == previous.positions && state.velocities == previous.velocities,
                  label + " input state changed");
            check(result.info.iterations == options.iterations, label + " sweep count");
            compare_points(result.state.positions, expected_x, 2.e-10, label + " positions");
            compare_points(result.state.velocities, expected_v, 4.e-8, label + " velocities");
            check(result.info.multipliers.size() == multipliers.size(), label + " multiplier count");
            double independent_raw = 0., independent_compliant = 0.;
            for (std::size_t element = 0; element < constraints.size(); ++element) {
                near(result.info.multipliers[element], multipliers[element], 2.e-11, label + " multiplier");
                const auto evaluated = cloth::evaluate_constraint(constraints[element], result.state.positions);
                check(evaluated.valid_geometry, label + " final valid geometry");
                const auto residual = evaluated.value + compliance/(dt*dt)*result.info.multipliers[element];
                independent_raw = std::max(independent_raw, std::abs(evaluated.value));
                independent_compliant = std::max(independent_compliant, std::abs(residual));
            }
            const auto group = static_cast<std::size_t>(cloth::ConstraintGroup::neohookean_energy);
            check(group == 3 && result.info.groups.size() == 4, label + " group layout");
            check(result.info.groups[group].count == constraints.size(), label + " energy count");
            near(result.info.groups[group].max_violation, raw, 2.e-10, label + " source raw C");
            near(result.info.groups[group].max_compliant_residual, compliant, 2.e-10, label + " source compliant R");
            near(result.info.groups[group].max_violation, independent_raw, 1.e-14, label + " independent raw C");
            near(result.info.groups[group].max_compliant_residual, independent_compliant, 1.e-14,
                 label + " independent compliant R");
            for (std::size_t empty_group = 0; empty_group < group; ++empty_group)
                check(result.info.groups[empty_group].count == 0, label + " unrelated constraint group");
            check(result.info.projection_converged == (independent_compliant <= options.energy_tolerance),
                  label + " convergence means compliant residual threshold");
            for (const auto node : fixed) {
                check(result.state.positions[node] == previous.positions[node], label + " fixed drift");
                check(result.state.velocities[node] == cloth::Vec3{0., 0., 0.}, label + " fixed velocity");
            }
            state = result.state;
            ++count_steps;
        }
        reader.expect("end_trajectory");
    }
    check(names.empty() && count_steps == 32, "Incomplete trajectory coverage");
}
} // namespace

int main(int argc, char** argv) {
    try {
        check(argc == 2, "Expected task-05A reference path");
        Reader reader(argv[1]);
        reader.expect("CLOTH_TASK05A_REFERENCE");
        check(reader.get<int>() == 1, "Unsupported task-05A fixture version");
        locals(reader);
        trajectories(reader);
        reader.expect("end");
        reader.finish();
        std::cout << "6 local energy cases / 12 projections / 4 trajectories / 32 Python energy-XPBD steps passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
