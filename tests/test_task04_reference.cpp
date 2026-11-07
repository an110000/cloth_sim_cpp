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
        check(static_cast<bool>(file), "Cannot open task-04 reference");
        std::string line, all;
        while (std::getline(file, line)) if (!line.empty() && line.front() != '#') all += line + '\n';
        stream_.str(all);
    }
    template<class T> T get() {
        T value{};
        if (!(stream_ >> value)) throw std::runtime_error("Truncated or invalid task-04 reference");
        return value;
    }
    void expect(const std::string& tag) {
        const auto actual = get<std::string>();
        check(actual == tag, "Expected " + tag + ", got " + actual);
    }
    void finish() {
        std::string extra;
        check(!(stream_ >> extra), "Trailing task-04 reference data");
    }
private:
    std::istringstream stream_;
};

void near(double actual, double expected, double absolute, const std::string& label) {
    const double tolerance = absolute + 1.e-8 * std::abs(expected);
    if (!std::isfinite(actual) || !std::isfinite(expected) || std::abs(actual - expected) > tolerance) {
        std::ostringstream out;
        out.precision(17);
        out << label << ": actual=" << actual << " expected=" << expected
            << " error=" << std::abs(actual - expected) << " tolerance=" << tolerance;
        throw std::runtime_error(out.str());
    }
}

std::vector<double> read_scalars(Reader& reader, const std::string& tag, std::size_t count) {
    reader.expect(tag);
    std::vector<double> result(count);
    for (auto& value : result) value = reader.get<double>();
    return result;
}

std::vector<cloth::Vec3> read_points(Reader& reader, const std::string& tag, std::size_t count) {
    reader.expect(tag);
    std::vector<cloth::Vec3> result(count);
    for (auto& point : result) for (auto& value : point) value = reader.get<double>();
    return result;
}

void compare_points(const std::vector<cloth::Vec3>& actual,
                    const std::vector<cloth::Vec3>& expected, double tolerance,
                    const std::string& label) {
    check(actual.size() == expected.size(), label + " vertex count");
    for (std::size_t i = 0; i < actual.size(); ++i)
        for (std::size_t k = 0; k < 3; ++k)
            near(actual[i][k], expected[i][k], tolerance,
                 label + "[" + std::to_string(i) + "][" + std::to_string(k) + "]");
}

cloth::ConstraintGroup parse_group(const std::string& name) {
    if (name == "structural") return cloth::ConstraintGroup::structural;
    if (name == "shear") return cloth::ConstraintGroup::shear;
    if (name == "bending") return cloth::ConstraintGroup::bending;
    throw std::runtime_error("Unknown reference group " + name);
}

cloth::Constraint read_constraint(Reader& reader) {
    reader.expect("constraint");
    const auto type = reader.get<std::string>();
    const auto group = parse_group(reader.get<std::string>());
    if (type == "distance") {
        check(group != cloth::ConstraintGroup::bending, "Distance cannot use bending group");
        const auto i = reader.get<std::size_t>();
        const auto j = reader.get<std::size_t>();
        const auto rest = reader.get<double>();
        const auto compliance = reader.get<double>();
        return cloth::DistanceConstraint{i, j, rest, compliance, 0.0, group};
    }
    check(type == "bending" && group == cloth::ConstraintGroup::bending, "Invalid reference constraint type");
    std::array<std::size_t, 4> nodes{};
    for (auto& node : nodes) node = reader.get<std::size_t>();
    const auto rest = reader.get<double>();
    const auto compliance = reader.get<double>();
    return cloth::DihedralConstraint{nodes, rest, compliance, 0.0};
}

void compare_constraints(const std::vector<cloth::Constraint>& actual,
                         const std::vector<cloth::Constraint>& expected) {
    check(actual.size() == expected.size(), "Ordered builder constraint count");
    for (std::size_t index = 0; index < actual.size(); ++index) {
        const auto& a = actual[index];
        const auto& b = expected[index];
        const auto label = "Ordered constraint " + std::to_string(index);
        check(a.index() == b.index() && cloth::constraint_group(a) == cloth::constraint_group(b), label + " type/group");
        if (const auto* distance = std::get_if<cloth::DistanceConstraint>(&a)) {
            const auto& expected_distance = std::get<cloth::DistanceConstraint>(b);
            check(distance->i == expected_distance.i && distance->j == expected_distance.j, label + " nodes");
            near(distance->rest_length, expected_distance.rest_length, 1.e-13, label + " length");
        } else {
            const auto& angle = std::get<cloth::DihedralConstraint>(a);
            const auto& expected_angle = std::get<cloth::DihedralConstraint>(b);
            check(angle.nodes == expected_angle.nodes, label + " oriented nodes");
            near(angle.rest_angle, expected_angle.rest_angle, 1.e-13, label + " rest angle");
        }
        near(cloth::constraint_compliance(a), cloth::constraint_compliance(b), 1.e-16, label + " compliance");
    }
}

std::vector<cloth::Constraint> read_constraints(Reader& reader) {
    reader.expect("constraints");
    const auto count = reader.get<std::size_t>();
    check(count > 0 && count <= 1000, "Invalid constraint count");
    std::vector<cloth::Constraint> result;
    for (std::size_t i = 0; i < count; ++i) result.push_back(read_constraint(reader));
    return result;
}

struct Grid {
    cloth::Mesh mesh;
    std::size_t nx;
    std::size_t ny;
};

Grid read_mesh(Reader& reader) {
    reader.expect("grid");
    const auto nx = reader.get<std::size_t>(), ny = reader.get<std::size_t>();
    check(nx > 0 && nx < 10 && ny > 0 && ny < 10, "Invalid reference grid");
    reader.expect("vertices");
    const auto vertices = reader.get<std::size_t>();
    reader.expect("triangles");
    const auto triangles = reader.get<std::size_t>();
    check(vertices == (nx+1)*(ny+1) && triangles == 2*nx*ny, "Reference grid count mismatch");
    cloth::Mesh mesh;
    mesh.reference_positions.resize(vertices);
    mesh.triangles.resize(triangles);
    reader.expect("reference");
    for (auto& point : mesh.reference_positions) for (auto& value : point) value = reader.get<double>();
    mesh.positions = read_points(reader, "rest_positions", vertices);
    reader.expect("topology");
    for (auto& triangle : mesh.triangles) for (auto& vertex : triangle) vertex = reader.get<std::size_t>();
    return {std::move(mesh), nx, ny};
}

void check_evaluation(const cloth::Constraint& constraint, const std::vector<cloth::Vec3>& points,
                      double expected_value, const std::vector<cloth::Vec3>& gradients,
                      const std::string& label) {
    const auto evaluation = cloth::evaluate_constraint(constraint, points);
    check(evaluation.valid_geometry && evaluation.count == gradients.size(), label + " geometry/count");
    near(evaluation.value, expected_value, 2.e-11, label + " value");
    for (std::size_t local = 0; local < gradients.size(); ++local)
        for (std::size_t axis = 0; axis < 3; ++axis)
            near(evaluation.gradients[local][axis], gradients[local][axis], 2.e-10, label + " gradient");
}

} // namespace

int main(int argc, char** argv) {
    try {
        check(argc == 2, "Expected task-04 reference path");
        Reader reader(argv[1]);
        reader.expect("CLOTH_TASK04_REFERENCE");
        check(reader.get<int>() == 1, "Unsupported task-04 reference version");
        reader.expect("locals");
        check(reader.get<std::size_t>() == 7, "Expected seven local reference cases");
        std::set<std::string> local_names;
        for (std::size_t local_index = 0; local_index < 7; ++local_index) {
            reader.expect("local");
            const auto name = reader.get<std::string>();
            check(local_names.insert(name).second, "Duplicate local case");
            reader.expect("type");
            const auto kind = reader.get<std::string>();
            reader.expect("vertices");
            const auto vertices = reader.get<std::size_t>();
            check((kind == "distance" && vertices == 2) || (kind == "bending" && vertices == 4), "Invalid local case kind");
            auto constraint = read_constraint(reader);
            reader.expect("dt");
            const auto dt = reader.get<double>();
            reader.expect("initial_lambda");
            const auto initial_lambda = reader.get<double>();
            std::visit([&](auto& c) { c.lambda_value = initial_lambda; }, constraint);
            const auto inverse = read_scalars(reader, "inverse_masses", vertices);
            auto x = read_points(reader, "positions", vertices);
            reader.expect("value");
            const auto value = reader.get<double>();
            const auto gradient = read_points(reader, "gradient", vertices);
            check_evaluation(constraint, x, value, gradient, name);
            reader.expect("projections");
            check(reader.get<std::size_t>() == 2, "Expected two accumulated local projections");
            for (std::size_t index = 1; index <= 2; ++index) {
                reader.expect("projection");
                check(reader.get<std::size_t>() == index, "Invalid projection order");
                reader.expect("delta_lambda");
                const auto delta = reader.get<double>();
                reader.expect("value_before");
                const auto before = reader.get<double>();
                reader.expect("lambda");
                const auto multiplier = reader.get<double>();
                const auto expected_x = read_points(reader, "positions", vertices);
                reader.expect("value_after");
                const auto after = reader.get<double>();
                const auto after_gradient = read_points(reader, "gradient_after", vertices);
                reader.expect("end_projection");
                const auto projected = cloth::project_constraint(constraint, x, inverse, dt);
                check(projected.valid_geometry, name + " projection invalid geometry");
                near(projected.delta_lambda, delta, 2.e-11, name + " delta lambda");
                near(projected.value_before, before, 2.e-11, name + " pre-projection value");
                near(cloth::constraint_multiplier(constraint), multiplier, 2.e-11, name + " accumulated lambda");
                compare_points(x, expected_x, 2.e-10, name + " projected positions");
                check_evaluation(constraint, x, after, after_gradient, name + " after projection");
            }
            reader.expect("reset_lambda");
            const auto reset_value = reader.get<double>();
            cloth::reset_constraint(constraint);
            near(cloth::constraint_multiplier(constraint), reset_value, 0.0, name + " reset");
            reader.expect("end_local");
        }
        reader.expect("builders");
        check(reader.get<std::size_t>() == 2, "Expected two ordered builders");
        for (std::size_t b = 0; b < 2; ++b) {
            reader.expect("builder");
            const auto name = reader.get<std::string>();
            auto grid = read_mesh(reader);
            const auto compliance = read_scalars(reader, "compliances", 3);
            const auto expected = read_constraints(reader);
            const auto actual = cloth::build_rectangular_constraints(grid.mesh, grid.nx, grid.ny,
                compliance[0], compliance[1], compliance[2]);
            compare_constraints(actual, expected);
            reader.expect("end_builder");
            check(!name.empty(), "Empty builder label");
        }
        reader.expect("trajectories");
        check(reader.get<std::size_t>() == 4, "Expected four trajectories");
        std::set<std::string> expected_names{"rest", "freefall", "pinned_deformed", "bending_folded"};
        std::size_t total_steps = 0;
        for (std::size_t t = 0; t < 4; ++t) {
            reader.expect("trajectory");
            const auto name = reader.get<std::string>();
            check(expected_names.erase(name) == 1, "Duplicate or unknown trajectory");
            auto grid = read_mesh(reader);
            const auto density_thickness = read_scalars(reader, "density_thickness", 2);
            const auto masses = read_scalars(reader, "physical_masses", grid.mesh.positions.size());
            const auto assembled_masses = cloth::lumped_masses(grid.mesh, density_thickness[0], density_thickness[1]);
            for (std::size_t i = 0; i < masses.size(); ++i) near(assembled_masses[i], masses[i], 1.e-13, "Nodal mass");
            const auto compliance = read_scalars(reader, "compliances", 3);
            reader.expect("iterations");
            cloth::XpbdOptions options;
            options.iterations = reader.get<std::size_t>();
            reader.expect("gravity");
            const cloth::Vec3 gravity{reader.get<double>(), reader.get<double>(), reader.get<double>()};
            reader.expect("dt");
            const auto dt = reader.get<double>();
            reader.expect("fixed");
            const auto fixed_count = reader.get<std::size_t>();
            check(fixed_count <= masses.size(), "Invalid fixed count");
            std::vector<std::size_t> fixed(fixed_count);
            for (auto& vertex : fixed) vertex = reader.get<std::size_t>();
            cloth::SimulationState state{read_points(reader, "positions", masses.size()),
                                          read_points(reader, "velocities", masses.size())};
            const auto expected_constraints = read_constraints(reader);
            const auto constraints = cloth::build_rectangular_constraints(grid.mesh, grid.nx, grid.ny,
                compliance[0], compliance[1], compliance[2]);
            compare_constraints(constraints, expected_constraints);
            const cloth::XpbdSolver solver(masses, constraints, gravity, options);
            reader.expect("steps");
            check(reader.get<std::size_t>() == 8, "Expected eight time steps");
            for (std::size_t step = 1; step <= 8; ++step) {
                reader.expect("step");
                check(reader.get<std::size_t>() == step, "Invalid trajectory step order");
                const auto expected_x = read_points(reader, "positions", masses.size());
                const auto expected_v = read_points(reader, "velocities", masses.size());
                const auto multipliers = read_scalars(reader, "multipliers", constraints.size());
                const auto raw = read_scalars(reader, "raw_errors", 3);
                const auto compliant = read_scalars(reader, "compliant_errors", 3);
                reader.expect("end_step");
                const auto previous = state;
                const auto result = solver.step(state, fixed, dt);
                const auto label = name + " step " + std::to_string(step);
                check(result.info.finite_state && result.info.committed &&
                      result.info.failure_reason == cloth::XpbdFailure::none, label + " failed");
                check(state.positions == previous.positions && state.velocities == previous.velocities, label + " input mutated");
                check(result.info.iterations == options.iterations, label + " sweep count");
                compare_points(result.state.positions, expected_x, 2.e-9, label + " positions");
                compare_points(result.state.velocities, expected_v, 2.e-7, label + " velocities");
                check(result.info.multipliers.size() == multipliers.size(), label + " multiplier count");
                std::array<double, 3> independent_raw{}, independent_compliant{};
                std::array<std::size_t, 3> counts{};
                for (std::size_t index = 0; index < constraints.size(); ++index) {
                    near(result.info.multipliers[index], multipliers[index], 2.e-10, label + " multiplier");
                    const auto evaluation = cloth::evaluate_constraint(constraints[index], result.state.positions);
                    check(evaluation.valid_geometry, label + " invalid final geometry");
                    const auto group = static_cast<std::size_t>(cloth::constraint_group(constraints[index]));
                    const auto residual = evaluation.value + cloth::constraint_compliance(constraints[index]) /
                        (dt * dt) * result.info.multipliers[index];
                    independent_raw[group] = std::max(independent_raw[group], std::abs(evaluation.value));
                    independent_compliant[group] = std::max(independent_compliant[group], std::abs(residual));
                    ++counts[group];
                }
                bool expected_converged = true;
                for (std::size_t group = 0; group < 3; ++group) {
                    check(result.info.groups[group].count == counts[group], label + " group count");
                    near(result.info.groups[group].max_violation, raw[group], 2.e-9, label + " raw residual");
                    near(result.info.groups[group].max_compliant_residual, compliant[group], 2.e-9, label + " compliant residual");
                    near(result.info.groups[group].max_violation, independent_raw[group], 1.e-14, label + " independent raw");
                    near(result.info.groups[group].max_compliant_residual, independent_compliant[group], 1.e-14, label + " independent compliant");
                    expected_converged = expected_converged && independent_compliant[group] <=
                        (group == 2 ? options.angle_tolerance : options.distance_tolerance);
                }
                check(result.info.projection_converged == expected_converged, label + " projection convergence label");
                for (const auto fixed_vertex : fixed) {
                    check(result.state.positions[fixed_vertex] == previous.positions[fixed_vertex], label + " fixed drift");
                    check(result.state.velocities[fixed_vertex] == cloth::Vec3{0.,0.,0.}, label + " fixed velocity");
                }
                state = result.state;
                ++total_steps;
            }
            reader.expect("end_trajectory");
        }
        reader.expect("end");
        reader.finish();
        check(expected_names.empty() && total_steps == 32, "Incomplete trajectory coverage");
        std::cout << "7 local cases / 14 projections / 2 ordered builders / 32 Python XPBD steps passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
