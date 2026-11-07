#include "cloth/implicit.hpp"

#include "cloth/assembly.hpp"
#include "cloth/dofs.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const std::string& message) {
    check(std::isfinite(actual) && std::isfinite(expected) &&
              std::abs(actual - expected) <= tolerance, message);
}

template<class F> void invalid(F action, const std::string& message) {
    bool caught = false;
    try { action(); } catch (const std::invalid_argument&) { caught = true; }
    check(caught, message);
}

cloth::SimulationState at_rest(const cloth::Mesh& mesh) {
    return {mesh.positions, std::vector<cloth::Vec3>(mesh.positions.size(), {0.0, 0.0, 0.0})};
}

void check_unchanged(const cloth::SimulationState& actual,
                     const cloth::SimulationState& expected, const std::string& message) {
    check(actual.positions == expected.positions && actual.velocities == expected.velocities, message);
}

std::vector<double> flatten(const std::vector<cloth::Vec3>& points) {
    std::vector<double> values;
    for (const auto& point : points) values.insert(values.end(), point.begin(), point.end());
    return values;
}

double final_residual(const cloth::ImplicitEulerSolver& solver,
                      const cloth::Mesh& geometry, const cloth::Material& material,
                      const cloth::SimulationState& initial,
                      const cloth::SimulationState& final,
                      const std::vector<std::size_t>& fixed, double dt) {
    auto mesh = geometry;
    mesh.positions = final.positions;
    const auto elastic = cloth::assemble_elastic(mesh, material);
    const auto x = flatten(final.positions);
    const auto x0 = flatten(initial.positions);
    const auto v0 = flatten(initial.velocities);
    std::vector<double> displacement(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) displacement[i] = x[i] - x0[i] - dt * v0[i];
    auto residual = solver.mass_matrix().multiply(displacement);
    for (std::size_t i = 0; i < residual.size(); ++i) {
        residual[i] += dt * dt * (elastic.gradient[i] - solver.external_force()[i]);
    }
    double length = 0.0;
    for (const auto value : cloth::restrict_to_free(residual, cloth::DofMap(initial.positions.size(), fixed))) {
        length = std::hypot(length, value);
    }
    return length;
}

void static_and_all_fixed() {
    const auto mesh = cloth::make_rectangular_mesh(1, 1, 1.0, 1.0);
    const cloth::ImplicitEulerSolver solver(mesh, {}, 1000.0, cloth::MassMode::lumped, {0.0, 0.0, 0.0});
    auto state = at_rest(mesh);
    const auto original = state;
    const auto result = solver.step(state, {2, 3}, 0.01);
    check(result.info.solver_converged && result.info.committed && result.info.finite_state, "Static convergence");
    check(result.info.failure_reason == cloth::StepFailure::none && result.info.iterations == 0, "Static zero updates");
    check_unchanged(result.state, original, "Static state preserved");
    check_unchanged(state, original, "Input immutable on success");

    for (auto& velocity : state.velocities) velocity = {0.5, -0.3, 0.7};
    const auto frozen = solver.step(state, {3, 1, 0, 2}, 0.01);
    check(frozen.info.committed && frozen.info.residual_norm == 0.0, "All-fixed empty solve");
    check(frozen.info.iterations == 0 && frozen.state.positions == state.positions, "All-fixed no movement");
    for (const auto& velocity : frozen.state.velocities) check(velocity == cloth::Vec3{0.0, 0.0, 0.0}, "All-fixed zero velocities");
}

void ballistic(cloth::MassMode mode) {
    const auto mesh = cloth::make_rectangular_mesh(2, 1, 1.3, 0.7);
    const cloth::Vec3 gravity{0.2, -9.81, 0.1};
    const cloth::Vec3 velocity0{0.3, 0.8, -0.2};
    const cloth::ImplicitEulerSolver solver(mesh, {}, 830.0, mode, gravity);
    auto state = at_rest(mesh);
    for (auto& velocity : state.velocities) velocity = velocity0;
    const double dt = 0.02;
    for (int step = 1; step <= 10; ++step) {
        const auto result = solver.step(state, {}, dt);
        check(result.info.committed && result.info.iterations == 0, "Uniform free fall solves at predictor");
        state = result.state;
        for (std::size_t i = 0; i < state.positions.size(); ++i) {
            for (std::size_t axis = 0; axis < 3; ++axis) {
                const double expected_position = mesh.positions[i][axis] + step * dt * velocity0[axis] +
                    0.5 * step * (step + 1) * dt * dt * gravity[axis];
                near(state.positions[i][axis], expected_position, 2.0e-10, "Implicit Euler ballistic position");
                near(state.velocities[i][axis], velocity0[axis] + step * dt * gravity[axis], 2.0e-9,
                     "Implicit Euler ballistic velocity; no double gravity");
            }
        }
    }
}

cloth::Mesh deformed() {
    auto mesh = cloth::make_rectangular_mesh(1, 1, 1.0, 0.8);
    for (auto& point : mesh.positions) point = {1.12 * point[0] + 0.06 * point[2],
                                               0.03 * point[0], 0.91 * point[2]};
    return mesh;
}

void elastic_and_fixed(cloth::MassMode mode) {
    const auto mesh = deformed();
    cloth::Material material;
    material.young = 2.0e5;
    const cloth::ImplicitEulerSolver solver(mesh, material, 900.0, mode);
    auto initial = at_rest(mesh);
    initial.velocities[2] = {0.7, 0.4, -0.2};
    initial.velocities[3] = {0.2, -0.1, 0.5};
    const auto original = initial;
    const auto result = solver.step(initial, {2, 3}, 0.01);
    check(result.info.committed && result.info.accepted_steps > 0, "Deformed step requires Newton");
    check(result.info.residual_norm <= result.info.residual_threshold, "Final residual convergence");
    near(final_residual(solver, mesh, material, initial, result.state, {2, 3}, 0.01),
         result.info.residual_norm, 2.0e-14, "Reported residual belongs to final accepted state");
    check_unchanged(initial, original, "Deformed input immutable");
    for (const auto fixed : {2U, 3U}) {
        check(result.state.positions[fixed] == initial.positions[fixed], "Fixed absolute positions preserved");
        check(result.state.velocities[fixed] == cloth::Vec3{0.0, 0.0, 0.0}, "Fixed velocities zero");
    }
    const auto& objectives = result.info.accepted_objectives;
    check(objectives.size() == result.info.accepted_steps + 1, "Objective trace count");
    for (std::size_t i = 1; i < objectives.size(); ++i) check(objectives[i] < objectives[i - 1], "Accepted potential strictly decreases");
    near(solver.incremental_potential(initial, result.state.positions, 0.01), result.info.final_objective,
         2.0e-14, "Final objective diagnostic");
}

void objective_derivative() {
    const auto mesh = deformed();
    const cloth::Material material;
    const cloth::ImplicitEulerSolver solver(mesh, material, 1000.0, cloth::MassMode::consistent);
    auto state = at_rest(mesh);
    state.velocities[0] = {0.1, -0.2, 0.3};
    const double dt = 0.013;
    auto candidate = state.positions;
    candidate[0][1] += 0.004;
    auto working = mesh;
    working.positions = candidate;
    const auto elastic = cloth::assemble_elastic(working, material);
    auto delta = flatten(candidate);
    const auto x0 = flatten(state.positions);
    const auto v0 = flatten(state.velocities);
    for (std::size_t i = 0; i < delta.size(); ++i) delta[i] -= x0[i] + dt * v0[i];
    auto expected_gradient = solver.mass_matrix().multiply(delta);
    for (std::size_t i = 0; i < delta.size(); ++i) {
        expected_gradient[i] += dt * dt * (elastic.gradient[i] - solver.external_force()[i]);
    }
    const double epsilon = 1.0e-6;
    for (std::size_t dof = 0; dof < delta.size(); ++dof) {
        auto plus = candidate;
        auto minus = candidate;
        plus[dof / 3][dof % 3] += epsilon;
        minus[dof / 3][dof % 3] -= epsilon;
        const double finite_difference = (solver.incremental_potential(state, plus, dt) -
                                           solver.incremental_potential(state, minus, dt)) / (2.0 * epsilon);
        near(finite_difference, expected_gradient[dof], 2.0e-9, "Incremental potential derivative");
    }
}

void final_update_and_indefinite_shift() {
    auto mesh = cloth::make_rectangular_mesh(1, 1, 1.0, 1.0);
    mesh.positions[0][0] -= 0.0001;
    mesh.positions[1][0] *= 1.0001;
    const auto initial = at_rest(mesh);
    cloth::NewtonOptions options;
    options.max_iterations = 0;
    const cloth::ImplicitEulerSolver check_initial(mesh, {}, 1000.0, cloth::MassMode::lumped,
                                                   {0.0, 0.0, 0.0}, options);
    const auto before = check_initial.step(initial, {2, 3}, 0.01);
    check(before.info.residual_norm > before.info.residual_threshold, "One-update fixture is initially unconverged");
    options.max_iterations = 1;
    const cloth::ImplicitEulerSolver one_update(mesh, {}, 1000.0, cloth::MassMode::lumped,
                                               {0.0, 0.0, 0.0}, options);
    const auto after = one_update.step(initial, {2, 3}, 0.01);
    check(after.info.committed && after.info.iterations == 1 && after.info.accepted_steps == 1,
          "The final allowed Newton update can converge");
    check(after.info.residual_norm <= after.info.residual_threshold,
          "Final-budget convergence uses the updated residual");

    // A strongly compressed membrane has negative geometric stiffness.
    // The asymmetric out-of-plane perturbation excites that part of the
    // system, requiring an actually accepted shifted descent direction.
    auto compressed = cloth::make_rectangular_mesh(1, 1, 1.0, 1.0);
    for (auto& point : compressed.positions) { point[0] *= 0.1; point[2] *= 0.1; }
    compressed.positions[0][1] = 0.02;
    compressed.positions[1][1] = -0.03;
    const auto compressed_initial = at_rest(compressed);
    options.max_iterations = 50;
    const cloth::ImplicitEulerSolver nonlinear(compressed, {}, 1000.0, cloth::MassMode::lumped,
                                               {0.0, 0.0, 0.0}, options);
    const auto recovered = nonlinear.step(compressed_initial, {2, 3}, 0.05);
    check(recovered.info.committed && recovered.info.shifted_steps > 0,
          "Nonconvex compressed membrane accepts a shifted descent step");
    near(final_residual(nonlinear, compressed, {}, compressed_initial, recovered.state, {2, 3}, 0.05),
         recovered.info.residual_norm, 2.0e-13, "Shifted solve final residual independently checked");
    check(recovered.info.residual_norm <= recovered.info.residual_threshold, "Shifted solve convergence");
}

void rollback_and_fallback() {
    const auto mesh = deformed();
    auto initial = at_rest(mesh);
    initial.velocities[0] = {0.1, 0.2, -0.3};
    cloth::NewtonOptions options;
    options.max_iterations = 0;
    const cloth::ImplicitEulerSolver no_updates(mesh, {}, 1000.0, cloth::MassMode::lumped,
                                                {0.0, 0.0, 0.0}, options);
    const auto failure = no_updates.step(initial, {2, 3}, 0.01);
    check(!failure.info.committed && !failure.info.solver_converged && failure.info.finite_state, "Failed finite state not converged");
    check(failure.info.failure_reason == cloth::StepFailure::iteration_limit, "Budget failure explicit");
    check_unchanged(failure.state, initial, "Failure rolls back positions and velocities");

    auto collapsed = initial;
    collapsed.positions[1] = collapsed.positions[0];
    const auto invalid_initial = no_updates.step(collapsed, {}, 0.01);
    check(invalid_initial.info.failure_reason == cloth::StepFailure::invalid_initial_configuration,
          "Collapsed initial state rejected structurally");
    check_unchanged(invalid_initial.state, collapsed, "Invalid geometry rollback");

    const auto regular = cloth::make_rectangular_mesh(1, 1, 1.0, 1.0);
    auto predictor_collapse = at_rest(regular);
    for (std::size_t i = 0; i < regular.positions.size(); ++i) {
        predictor_collapse.velocities[i][0] = -regular.positions[i][0] / 0.01;
    }
    const cloth::ImplicitEulerSolver fallback_solver(regular, {}, 1000.0, cloth::MassMode::lumped,
                                                     {0.0, 0.0, 0.0}, options);
    const auto fallback = fallback_solver.step(predictor_collapse, {}, 0.01);
    check(fallback.info.start_fallback && fallback.info.failure_reason == cloth::StepFailure::iteration_limit,
          "Collapsed predictor falls back before checking update budget");
    check_unchanged(fallback.state, predictor_collapse, "Fallback failure still rolls back input");

    options.max_iterations = 30;
    options.line_search_max_iterations = 1;
    options.armijo = 0.999999;
    options.diagonal_shifts = {0.0};
    const cloth::ImplicitEulerSolver no_backtracking(mesh, {}, 1000.0, cloth::MassMode::lumped,
                                                    {0.0, 0.0, 0.0}, options);
    const auto line_failure = no_backtracking.step(initial, {2, 3}, 0.01);
    check(line_failure.info.failure_reason == cloth::StepFailure::line_search_failed && !line_failure.info.committed,
          "Line-search failure explicit");
    check(line_failure.info.line_search_rejections > 0, "Rejected trial counted");
    check_unchanged(line_failure.state, initial, "Line-search failure rolls back");

    auto overflow = at_rest(regular);
    for (auto& velocity : overflow.velocities) velocity = {std::numeric_limits<double>::max(), 0.0, 0.0};
    const auto arithmetic_failure = fallback_solver.step(overflow, {}, 2.0);
    check(arithmetic_failure.info.failure_reason == cloth::StepFailure::non_finite_arithmetic,
          "Finite inputs with overflowing prediction fail structurally");
    check_unchanged(arithmetic_failure.state, overflow, "Arithmetic failure rollback");
}

void invalid_inputs() {
    const auto mesh = cloth::make_rectangular_mesh(1, 1, 1.0, 1.0);
    const cloth::ImplicitEulerSolver solver(mesh, {});
    const auto state = at_rest(mesh);
    invalid([&] { (void)solver.step(state, {}, 0.0); }, "Zero dt rejected");
    invalid([&] { (void)solver.step(state, {}, std::numeric_limits<double>::infinity()); }, "Infinite dt rejected");
    invalid([&] { (void)solver.step(state, {}, 1.0e-300); }, "Squared dt underflow rejected");
    invalid([&] { (void)solver.step(state, {0, 0}, 0.01); }, "Duplicate fixed rejected");
    invalid([&] { (void)solver.step(state, {4}, 0.01); }, "Out of range fixed rejected");
    auto bad = state;
    bad.velocities.pop_back();
    invalid([&] { (void)solver.step(bad, {}, 0.01); }, "Velocity size rejected");
    bad = state;
    bad.positions[0][1] = std::numeric_limits<double>::quiet_NaN();
    invalid([&] { (void)solver.step(bad, {}, 0.01); }, "NaN position rejected");
    cloth::NewtonOptions options;
    options.armijo = 1.0;
    invalid([&] { cloth::ImplicitEulerSolver bad_solver(mesh, {}, 1000.0, cloth::MassMode::lumped,
                                                      {0.0, 0.0, 0.0}, options); }, "Invalid Armijo rejected");
    options = {};
    options.diagonal_shifts = {0.0, 0.0};
    invalid([&] { cloth::ImplicitEulerSolver bad_solver(mesh, {}, 1000.0, cloth::MassMode::lumped,
                                                      {0.0, 0.0, 0.0}, options); }, "Duplicate shifts rejected");
}

} // namespace

int main() {
    try {
        static_and_all_fixed();
        ballistic(cloth::MassMode::lumped);
        ballistic(cloth::MassMode::consistent);
        elastic_and_fixed(cloth::MassMode::lumped);
        elastic_and_fixed(cloth::MassMode::consistent);
        objective_derivative();
        final_update_and_indefinite_shift();
        rollback_and_fallback();
        invalid_inputs();
        std::cout << "Implicit Euler residuals, ballistic motion, constraints and rollback passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
