#pragma once

#include "cloth/contact.hpp"
#include "cloth/material.hpp"
#include "cloth/mass.hpp"
#include "cloth/state.hpp"

#include <cstddef>
#include <limits>
#include <string>
#include <vector>

namespace cloth {

struct NewtonOptions {
    // Maximum attempted Newton updates; convergence is also checked after
    // the last accepted update. Zero permits only the initial residual check.
    std::size_t max_iterations = 30;
    double residual_tolerance = 1.0e-8;
    std::size_t line_search_max_iterations = 24;
    double backtrack_factor = 0.5;
    double armijo = 1.0e-4;
    double linear_relative_tolerance = 1.0e-10;
    // Factors multiply max(1, max(abs(diag(A_FF)))). A finite shift is not a
    // proof of positive definiteness; each direction must pass r.dot(d) < 0.
    std::vector<double> diagonal_shifts{0.0, 1.0e-10, 1.0e-8, 1.0e-6,
                                         1.0e-4, 1.0e-2, 1.0, 100.0};
};

enum class StepFailure {
    none,
    invalid_initial_configuration,
    non_finite_arithmetic,
    linear_solve_failed,
    line_search_failed,
    iteration_limit
};

const char* failure_name(StepFailure reason) noexcept;

struct StepInfo {
    // Describes the returned state (the original state after rollback), not
    // whether an intermediate Newton trial was geometrically admissible.
    bool finite_state = true;
    // Convergence of the discrete position equation before the split friction
    // update. Friction subsequently changes tangential velocity only.
    bool solver_converged = false;
    bool committed = false;
    StepFailure failure_reason = StepFailure::none;
    std::string message;
    std::size_t iterations = 0;
    std::size_t accepted_steps = 0;
    std::size_t line_search_rejections = 0;
    std::size_t degenerate_rejections = 0;
    std::size_t shifted_steps = 0;
    bool start_fallback = false;
    double residual_norm = std::numeric_limits<double>::infinity();
    double residual_threshold = 0.0;
    double final_objective = std::numeric_limits<double>::infinity();
    // Includes the feasible starting objective, then every accepted iterate.
    // With no free DOFs it contains the initial objective only.
    std::vector<double> accepted_objectives;
    ContactInfo contact;
};

struct StepResult {
    SimulationState state;
    StepInfo info;
};

// Serial implicit Euler with optional in-Newton ground penalty. Reference topology, material, mass and
// constant gravity are owned by the solver. The input state is never modified.
// Fixed positions are their values at the start of this step; fixed output
// velocities are zero. A failed solve returns the complete original state.
class ImplicitEulerSolver {
public:
    ImplicitEulerSolver(const Mesh& mesh, const Material& material,
                        double density = 1000.0,
                        MassMode mass_mode = MassMode::lumped,
                        Vec3 gravity = {0.0, -9.81, 0.0},
                        NewtonOptions options = {},
                        GroundContactOptions ground = {});

    StepResult step(const SimulationState& state,
                    const std::vector<std::size_t>& fixed_vertices,
                    double dt) const;

    const SparseMatrix& mass_matrix() const noexcept { return mass_; }
    const std::vector<double>& external_force() const noexcept { return force_; }
    const NewtonOptions& options() const noexcept { return options_; }
    const GroundContactOptions& ground_options() const noexcept { return ground_; }

    // Phi(z) = .5 (z-y)^T M (z-y) + dt^2 (E(z) + Uc(z) - f_g^T z), y=x+dt*v.
    // Pass the same fixed vertices as step() to exclude their ground penalty.
    // The caller remains responsible for preserving their candidate positions.
    // Units are kg*m^2, i.e. dt^2 times energy, rather than joules.
    double incremental_potential(const SimulationState& state,
                                 const std::vector<Vec3>& candidate,
                                 double dt,
                                 const std::vector<std::size_t>& fixed_vertices = {}) const;

private:
    Mesh geometry_;
    Material material_;
    Vec3 gravity_;
    NewtonOptions options_;
    GroundContactOptions ground_;
    // Physical row-sum nodal masses, also when mass_ is consistent.
    std::vector<double> lumped_masses_;
    SparseMatrix mass_;
    std::vector<double> force_;
};

} // namespace cloth
