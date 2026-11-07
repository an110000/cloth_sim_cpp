#include "cloth/implicit.hpp"

#include "cloth/assembly.hpp"
#include "cloth/dofs.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
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

std::vector<double> flatten(const std::vector<cloth::Vec3>& points) {
    std::vector<double> values;
    for (const auto& point : points) values.insert(values.end(), point.begin(), point.end());
    return values;
}

double kinetic_energy(const cloth::SparseMatrix& mass, const std::vector<cloth::Vec3>& velocities) {
    const auto v = flatten(velocities);
    const auto momentum = mass.multiply(v);
    long double energy = 0.0L;
    for (std::size_t i = 0; i < v.size(); ++i) energy += 0.5L * v[i] * momentum[i];
    return static_cast<double>(energy);
}

cloth::SimulationState uniform_state(const cloth::Mesh& mesh, double y, cloth::Vec3 velocity) {
    cloth::SimulationState state{mesh.positions, std::vector<cloth::Vec3>(mesh.positions.size(), velocity)};
    for (auto& point : state.positions) point[1] = y;
    return state;
}

void same_state(const cloth::SimulationState& a, const cloth::SimulationState& b, const std::string& message) {
    check(a.positions == b.positions && a.velocities == b.velocities, message);
}

double position_residual(const cloth::ImplicitEulerSolver& solver, const cloth::Mesh& geometry,
                         const cloth::Material& material, const cloth::SimulationState& initial,
                         const cloth::SimulationState& final, const std::vector<std::size_t>& fixed,
                         const std::vector<double>& masses, double dt) {
    auto mesh = geometry;
    mesh.positions = final.positions;
    const auto elastic = cloth::assemble_elastic(mesh, material);
    const auto penalty = cloth::evaluate_ground_penalty(final.positions, masses, fixed, dt,
                                                        solver.ground_options());
    const auto x = flatten(final.positions);
    const auto x0 = flatten(initial.positions);
    const auto v0 = flatten(initial.velocities);
    std::vector<double> displacement(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) displacement[i] = x[i] - x0[i] - dt * v0[i];
    auto residual = solver.mass_matrix().multiply(displacement);
    for (std::size_t i = 0; i < residual.size(); ++i) {
        const double ground_gradient = i % 3 == 1 ? penalty.gradient_y[i / 3] : 0.0;
        residual[i] += dt * dt * (elastic.gradient[i] + ground_gradient - solver.external_force()[i]);
    }
    double length = 0.0;
    for (const auto value : cloth::restrict_to_free(residual, cloth::DofMap(initial.positions.size(), fixed))) {
        length = std::hypot(length, value);
    }
    return length;
}

// Every node undergoes the same translation. Elasticity is exactly unchanged,
// and M*1=m_lumped for both mass modes, reducing the active contact solve to
// y=(y0+dt*vn+dt^2*g+eta*ground)/(1+eta).
void uniform_analytic_solution(cloth::MassMode mode) {
    const auto mesh = cloth::make_rectangular_mesh(1, 1, 1.0, 1.0);
    cloth::Material material;
    material.young = 1000.0;
    material.thickness = 0.01;
    cloth::NewtonOptions newton;
    newton.residual_tolerance = 1.0e-10;
    cloth::GroundContactOptions ground;
    ground.enabled = true;
    ground.penalty_eta = 20.0;
    ground.friction = {0.2, 0.3};
    const cloth::ImplicitEulerSolver solver(mesh, material, 4.0, mode, {0.0, -10.0, 0.0}, newton, ground);
    const auto initial = uniform_state(mesh, 0.02, {2.0, -1.0, 1.0});
    const auto original = initial;
    const double dt = 0.05;
    const double free_y = 0.02 - dt - 10.0 * dt * dt;
    const double expected_y = free_y / (1.0 + ground.penalty_eta);
    const double penetration = -expected_y;
    const double impulse_per_mass = ground.penalty_eta * penetration / dt;
    const double speed = std::hypot(2.0, 1.0);
    const double factor = 1.0 - ground.friction.kinetic * impulse_per_mass / speed;
    const auto masses = cloth::lumped_masses(mesh, 4.0, material.thickness);
    const auto result = solver.step(initial, {}, dt);
    check(result.info.committed && result.info.solver_converged, "Uniform penalty solve converges");
    same_state(initial, original, "Successful penalty step leaves input immutable");
    check(result.info.contact.active_count == initial.positions.size(), "Every penetrating node is active");
    check(result.info.contact.sliding_count == initial.positions.size(), "Uniform sliding classified");
    check(result.info.contact.sticking_count == 0, "Uniform sliding does not stick");
    double total_mass = 0.0;
    for (std::size_t i = 0; i < initial.positions.size(); ++i) {
        total_mass += masses[i];
        near(result.state.positions[i][0], initial.positions[i][0] + 2.0 * dt, 1.0e-10,
             "Tangential positions follow pre-friction implicit solve");
        near(result.state.positions[i][1], expected_y, 1.0e-10, "Penalty scalar position");
        near(result.state.positions[i][2], initial.positions[i][2] + dt, 1.0e-10, "Uniform Z translation");
        near(result.state.velocities[i][0], 2.0 * factor, 1.0e-9, "Sliding X velocity");
        near(result.state.velocities[i][2], factor, 1.0e-9, "Sliding Z velocity");
        near(result.state.velocities[i][1], (expected_y - 0.02) / dt, 1.0e-9,
             "Friction preserves the solved normal velocity");
        near(result.info.contact.normal_impulses[i], masses[i] * impulse_per_mass, 1.0e-10,
             "Penalty impulse dt*kappa*penetration uses physical nodal mass");
        near(result.info.contact.penetration[i], penetration, 1.0e-10, "Final penetration reported");
        check(result.info.contact.position_corrections[i] == cloth::Vec3{0.0, 0.0, 0.0},
              "Implicit contact does not project positions");
    }
    near(result.info.contact.penalty_energy,
         0.5 * ground.penalty_eta * total_mass * penetration * penetration / (dt * dt),
         1.0e-10, "Penalty energy in joules");
    near(result.info.contact.normal_impulse_total, total_mass * impulse_per_mass, 1.0e-10,
         "Total penalty impulse");
    near(result.info.contact.max_penetration_before, penetration, 1.0e-10, "Friction input penetration");
    near(result.info.contact.max_penetration_after, penetration, 1.0e-10, "Friction leaves penetration unchanged");
    near(result.info.contact.normal_kinetic_change, 0.0, 0.0, "No separate normal postprocessing stage");
    near(result.info.final_objective, solver.incremental_potential(initial, result.state.positions, dt),
         1.0e-12, "Reported final objective includes contact energy");
    near(position_residual(solver, mesh, material, initial, result.state, {}, masses, dt),
         result.info.residual_norm, 1.0e-12, "Residual contains contact gradient");
    check(result.info.residual_norm <= result.info.residual_threshold, "Pre-friction position equation converged");
    for (std::size_t i = 1; i < result.info.accepted_objectives.size(); ++i) {
        check(result.info.accepted_objectives[i] < result.info.accepted_objectives[i - 1],
              "Line search decreases the full contact objective");
    }

    // A uniform direction checks that the diagnostic potential really includes
    // the same penalty gradient as Newton (away from its branch boundary).
    auto candidate = initial.positions;
    for (auto& point : candidate) point[1] = -0.01;
    auto plus = candidate, minus = candidate;
    const double epsilon = 1.0e-6;
    for (auto& point : plus) point[1] += epsilon;
    for (auto& point : minus) point[1] -= epsilon;
    const double derivative = (solver.incremental_potential(initial, plus, dt) -
                               solver.incremental_potential(initial, minus, dt)) / (2.0 * epsilon);
    const double expected_derivative = total_mass * (-0.01 - (0.02 - dt) +
                                                    10.0 * dt * dt - 0.01 * ground.penalty_eta);
    near(derivative, expected_derivative, 1.0e-10, "Contact objective directional derivative");
}

void fixed_and_disabled() {
    const auto mesh = cloth::make_rectangular_mesh(1, 1, 1.0, 1.0);
    const auto state = uniform_state(mesh, -0.02, {1.0, -0.4, 0.3});
    cloth::GroundContactOptions ground;
    ground.enabled = true;
    const cloth::ImplicitEulerSolver solver(mesh, {}, 1000.0, cloth::MassMode::consistent,
                                             {0.0, 0.0, 0.0}, {}, ground);
    const std::vector<std::size_t> fixed{0, 1, 2, 3};
    const auto result = solver.step(state, fixed, 0.01);
    check(result.info.committed && result.info.iterations == 0, "All-fixed contact is an empty solve");
    check(result.state.positions == state.positions, "Fixed penetrations are not projected");
    for (const auto& velocity : result.state.velocities) check(velocity == cloth::Vec3{0.0, 0.0, 0.0}, "Fixed velocities zero");
    check(result.info.contact.fixed_penetration_count == 4 && result.info.contact.active_count == 0,
          "Fixed penetration has a separate count");
    near(result.info.contact.max_fixed_penetration, 0.02, 1.0e-15, "Fixed penetration depth");
    near(result.info.contact.penalty_energy, 0.0, 0.0, "Fixed nodes do not contribute penalty energy");
    near(result.info.contact.normal_impulse_total, 0.0, 0.0, "Fixed nodes receive no penalty impulse");
    for (const auto penetration : result.info.contact.penetration) near(penetration, 0.0, 0.0, "Free-only penetration array");
    near(result.info.final_objective, solver.incremental_potential(state, result.state.positions, 0.01, fixed),
         1.0e-12, "Diagnostic potential excludes specified fixed penalty nodes");
    check(solver.incremental_potential(state, result.state.positions, 0.01) > result.info.final_objective,
          "Diagnostic caller must supply fixed nodes to exclude their penalty");

    const cloth::ImplicitEulerSolver old_api(mesh, {}, 1000.0, cloth::MassMode::consistent);
    ground.enabled = false;
    ground.ground_y = 20.0;
    ground.restitution = 0.7;
    ground.stabilization_beta = 0.8;
    const cloth::ImplicitEulerSolver disabled(mesh, {}, 1000.0, cloth::MassMode::consistent,
                                               {0.0, -9.81, 0.0}, {}, ground);
    const auto a = old_api.step(state, {}, 0.01);
    const auto b = disabled.step(state, {}, 0.01);
    check(a.info.committed && b.info.committed, "Disabled regression steps commit");
    same_state(a.state, b.state, "Disabled contact preserves the original implicit trajectory exactly");
    check(a.info.accepted_objectives == b.info.accepted_objectives &&
          a.info.residual_norm == b.info.residual_norm, "Disabled Newton diagnostics unchanged");
    near(b.info.contact.normal_impulse_total, 0.0, 0.0, "Disabled contact impulse is zero");
}

void consistent_friction_energy_is_separate() {
    // Opposite neighboring velocities make consistent kinetic energy include
    // negative cross terms. Reducing just node 0 can increase this energy even
    // while the physical nodal (lumped) friction dissipation is positive.
    cloth::Mesh mesh;
    mesh.reference_positions = {{0.0, 0.0}, {1.0, 0.0}, {0.0, 1.0}};
    mesh.positions = {{0.0, -0.1, 0.0}, {1.0, 0.1, 0.0}, {0.0, 0.1, 1.0}};
    mesh.triangles = {{0, 1, 2}};
    cloth::Material material;
    material.young = 1.0e-5;
    material.thickness = 1.0;
    cloth::GroundContactOptions ground;
    ground.enabled = true;
    ground.penalty_eta = 1.0;
    ground.friction = {0.05, std::nullopt};
    const cloth::ImplicitEulerSolver solver(mesh, material, 1.0, cloth::MassMode::consistent,
                                             {0.0, 0.0, 0.0}, {}, ground);
    cloth::SimulationState initial{mesh.positions, {{1.0, 0.0, 0.0}, {-4.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
    const double dt = 0.01;
    const auto result = solver.step(initial, {}, dt);
    check(result.info.committed, "Nonuniform consistent-mass contact converges");
    auto before = result.state.velocities;
    for (std::size_t i = 0; i < before.size(); ++i) {
        for (std::size_t axis = 0; axis < 3; ++axis) {
            before[i][axis] = (result.state.positions[i][axis] - initial.positions[i][axis]) / dt;
        }
    }
    const double actual_change = kinetic_energy(solver.mass_matrix(), result.state.velocities) -
                                 kinetic_energy(solver.mass_matrix(), before);
    near(result.info.contact.friction_kinetic_change, actual_change, 1.0e-12,
         "Friction kinetic change uses the actual consistent mass matrix");
    check(result.info.contact.friction_dissipation > 0.0, "Nodal Coulomb friction dissipates lumped energy");
    check(result.info.contact.friction_kinetic_change > 0.0,
          "Nodal friction may increase consistent kinetic energy through cross terms");
    const auto masses = cloth::lumped_masses(mesh, 1.0, material.thickness);
    double nodal_dissipation = 0.0;
    for (std::size_t i = 0; i < before.size(); ++i) {
        nodal_dissipation += 0.5 * masses[i] *
            (before[i][0] * before[i][0] + before[i][2] * before[i][2] -
             result.state.velocities[i][0] * result.state.velocities[i][0] -
             result.state.velocities[i][2] * result.state.velocities[i][2]);
        near(result.state.velocities[i][1], before[i][1], 0.0, "Consistent friction leaves solved normal velocity");
    }
    near(result.info.contact.friction_dissipation, nodal_dissipation, 1.0e-12,
         "Reported friction dissipation has the declared nodal mass metric");
}

void proximity_and_rollback() {
    const auto mesh = cloth::make_rectangular_mesh(1, 1, 1.0, 1.0);
    cloth::GroundContactOptions ground;
    ground.enabled = true;
    ground.proximity_tolerance = 1.0;
    const cloth::ImplicitEulerSolver solver(mesh, {}, 1000.0, cloth::MassMode::lumped,
                                             {0.0, 0.0, 0.0}, {}, ground);
    const auto above = uniform_state(mesh, 1.0e-6, {0.0, 0.0, 0.0});
    const auto inactive = solver.step(above, {}, 0.01);
    check(inactive.info.committed && inactive.info.contact.active_count == 0,
          "Proximity does not create an implicit penalty above the plane");
    same_state(inactive.state, above, "Inactive plane leaves static state unchanged");

    cloth::NewtonOptions zero_budget;
    zero_budget.max_iterations = 0;
    const cloth::ImplicitEulerSolver limited(mesh, {}, 1000.0, cloth::MassMode::lumped,
                                              {0.0, 0.0, 0.0}, zero_budget, ground);
    const auto initial = uniform_state(mesh, -0.01, {1.0, -0.2, 0.0});
    const auto failed = limited.step(initial, {}, 0.01);
    check(!failed.info.committed && !failed.info.solver_converged && failed.info.finite_state,
          "Contact update budget failure rolls back without claiming convergence");
    check(failed.info.failure_reason == cloth::StepFailure::iteration_limit,
          "Contact budget failure has explicit reason");
    same_state(failed.state, initial, "Contact failure returns full original state including velocity");

    ground.restitution = 0.1;
    invalid([&] { cloth::ImplicitEulerSolver invalid_solver(mesh, {}, 1000.0, cloth::MassMode::lumped,
                                                            {0.0, 0.0, 0.0}, {}, ground); },
            "Implicit penalty rejects unsupported restitution instead of ignoring it");
    ground.restitution = 0.0;
    ground.stabilization_beta = 0.1;
    invalid([&] { cloth::ImplicitEulerSolver invalid_solver(mesh, {}, 1000.0, cloth::MassMode::lumped,
                                                            {0.0, 0.0, 0.0}, {}, ground); },
            "Implicit penalty rejects unsupported stabilization beta");
}

} // namespace

int main() {
    try {
        uniform_analytic_solution(cloth::MassMode::lumped);
        uniform_analytic_solution(cloth::MassMode::consistent);
        fixed_and_disabled();
        consistent_friction_energy_is_separate();
        proximity_and_rollback();
        std::cout << "Implicit ground penalty, split friction, mass metrics and rollback passed\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
