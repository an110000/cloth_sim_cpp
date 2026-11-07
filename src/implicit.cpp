#include "cloth/implicit.hpp"

#include "cloth/assembly.hpp"
#include "cloth/dofs.hpp"
#include "cloth/linear_solver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace cloth {
namespace {

bool finite_vector(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(),
                       [](double value) { return std::isfinite(value); });
}

std::vector<double> flatten(const std::vector<Vec3>& values) {
    std::vector<double> result;
    result.reserve(3 * values.size());
    for (const auto& point : values) result.insert(result.end(), point.begin(), point.end());
    return result;
}

std::vector<Vec3> unflatten(const std::vector<double>& values) {
    std::vector<Vec3> result(values.size() / 3);
    for (std::size_t i = 0; i < values.size(); ++i) result[i / 3][i % 3] = values[i];
    return result;
}

void validate_points(const std::vector<Vec3>& points, std::size_t count,
                     const char* description) {
    if (points.size() != count) throw std::invalid_argument(description);
    for (const auto& point : points) {
        for (const double value : point) {
            if (!std::isfinite(value)) throw std::invalid_argument(description);
        }
    }
}

double validate_step(const SimulationState& state, std::size_t count, double dt) {
    if (!std::isfinite(dt) || dt <= 0.0) {
        throw std::invalid_argument("Time step must be finite and positive");
    }
    const double dt2 = dt * dt;
    if (!std::isfinite(dt2) || dt2 == 0.0) {
        throw std::invalid_argument("Squared time step is outside double precision range");
    }
    validate_points(state.positions, count, "State positions have wrong size or non-finite values");
    validate_points(state.velocities, count, "State velocities have wrong size or non-finite values");
    return dt2;
}

double checked_dot(const std::vector<double>& a, const std::vector<double>& b) {
    long double sum = 0.0L;
    for (std::size_t i = 0; i < a.size(); ++i) {
        sum += static_cast<long double>(a[i]) * b[i];
    }
    const double value = static_cast<double>(sum);
    if (!std::isfinite(value)) throw std::overflow_error("Non-finite vector dot product");
    return value;
}

double norm(const std::vector<double>& values) {
    double result = 0.0;
    for (const double value : values) result = std::hypot(result, value);
    if (!std::isfinite(result)) throw std::overflow_error("Non-finite vector norm");
    return result;
}

struct FirstOrder {
    double energy = 0.0;
    std::vector<double> gradient;
};

// evaluate_triangle deliberately avoids the 9x9 Hessian. Line-search trials
// reuse this lightweight path instead of assembling a sparse matrix.
FirstOrder elastic_first_order(const Mesh& geometry, const Material& material,
                              const std::vector<double>& positions) {
    FirstOrder result{0.0, std::vector<double>(positions.size(), 0.0)};
    for (const auto& triangle : geometry.triangles) {
        std::array<Vec2, 3> reference{};
        std::array<Vec3, 3> current{};
        for (std::size_t vertex = 0; vertex < 3; ++vertex) {
            reference[vertex] = geometry.reference_positions[triangle[vertex]];
            for (std::size_t axis = 0; axis < 3; ++axis) {
                current[vertex][axis] = positions[3 * triangle[vertex] + axis];
            }
        }
        const auto local = evaluate_triangle(reference, current, material);
        result.energy += local.energy;
        for (std::size_t vertex = 0; vertex < 3; ++vertex) {
            for (std::size_t axis = 0; axis < 3; ++axis) {
                result.gradient[3 * triangle[vertex] + axis] += local.gradient[vertex][axis];
            }
        }
    }
    if (!std::isfinite(result.energy) || !finite_vector(result.gradient)) {
        throw std::overflow_error("Non-finite elastic energy or gradient");
    }
    return result;
}

std::vector<double> prediction(const std::vector<double>& x0,
                               const std::vector<double>& v0, double dt) {
    auto y = x0;
    for (std::size_t i = 0; i < y.size(); ++i) y[i] += dt * v0[i];
    if (!finite_vector(y)) throw std::overflow_error("Non-finite inertial prediction");
    return y;
}

double objective(const SparseMatrix& mass, const std::vector<double>& force,
                 const std::vector<double>& y, const std::vector<double>& x,
                 double elastic_energy, double contact_energy, double dt2) {
    auto displacement = x;
    for (std::size_t i = 0; i < x.size(); ++i) displacement[i] -= y[i];
    if (!finite_vector(displacement)) throw std::overflow_error("Non-finite inertial displacement");
    const auto momentum = mass.multiply(displacement);
    const long double value = 0.5L * checked_dot(displacement, momentum) +
        static_cast<long double>(dt2) * (static_cast<long double>(elastic_energy) +
                                        contact_energy - checked_dot(force, x));
    const double result = static_cast<double>(value);
    if (!std::isfinite(result)) throw std::overflow_error("Non-finite incremental potential");
    return result;
}

// Symmetrize only the Newton system to remove roundoff asymmetry, leaving the
// physical Hessian API and its raw derivative tests unchanged.
SparseMatrix newton_system(const SparseMatrix& mass, const SparseMatrix& hessian,
                           const std::vector<double>& contact_diagonal,
                           const DofMap& map, double dt2) {
    std::vector<Triplet> entries;
    auto append = [&](const SparseMatrix& matrix, double scale) {
        for (std::size_t row = 0; row < matrix.rows(); ++row) {
            for (std::size_t k = matrix.row_offsets()[row]; k < matrix.row_offsets()[row + 1]; ++k) {
                const auto column = matrix.column_indices()[k];
                const double value = 0.5 * scale * matrix.values()[k];
                if (!std::isfinite(value)) throw std::overflow_error("Non-finite Newton matrix entry");
                entries.push_back({row, column, value});
                entries.push_back({column, row, value});
            }
        }
    };
    append(mass, 1.0);
    append(hessian, dt2);
    for (std::size_t i = 0; i < contact_diagonal.size(); ++i) {
        const double value = dt2 * contact_diagonal[i];
        if (!std::isfinite(value)) throw std::overflow_error("Non-finite ground Newton diagonal");
        if (value != 0.0) entries.push_back({3 * i + 1, 3 * i + 1, value});
    }
    return restrict_to_free(SparseMatrix::from_triplets(mass.rows(), mass.cols(),
                                                       std::move(entries)), map);
}

SparseMatrix shifted(const SparseMatrix& matrix, double diagonal_shift) {
    std::vector<Triplet> entries;
    entries.reserve(matrix.non_zeros() + matrix.rows());
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        for (std::size_t k = matrix.row_offsets()[row]; k < matrix.row_offsets()[row + 1]; ++k) {
            entries.push_back({row, matrix.column_indices()[k], matrix.values()[k]});
        }
        entries.push_back({row, row, diagonal_shift});
    }
    return SparseMatrix::from_triplets(matrix.rows(), matrix.cols(), std::move(entries));
}

} // namespace

const char* failure_name(StepFailure reason) noexcept {
    switch (reason) {
    case StepFailure::none: return "none";
    case StepFailure::invalid_initial_configuration: return "invalid_initial_configuration";
    case StepFailure::non_finite_arithmetic: return "non_finite_arithmetic";
    case StepFailure::linear_solve_failed: return "linear_solve_failed";
    case StepFailure::line_search_failed: return "line_search_failed";
    case StepFailure::iteration_limit: return "iteration_limit";
    }
    return "unknown";
}

ImplicitEulerSolver::ImplicitEulerSolver(const Mesh& mesh, const Material& material,
                                       double density, MassMode mass_mode,
                                       Vec3 gravity, NewtonOptions options,
                                       GroundContactOptions ground)
    : geometry_(mesh), material_(material), gravity_(gravity), options_(std::move(options)),
      ground_(std::move(ground)), lumped_masses_(lumped_masses(mesh, density, material.thickness)),
      mass_(assemble_mass_matrix(mesh, density, material.thickness, mass_mode)) {
    validate_ground_contact_options(ground_);
    if (ground_.enabled && (ground_.restitution != 0.0 || ground_.stabilization_beta != 0.0)) {
        throw std::invalid_argument("Implicit penalty contact does not use restitution or stabilization beta");
    }
    if (!std::isfinite(material.young) || material.young <= 0.0 ||
        !std::isfinite(material.poisson) || material.poisson <= -1.0 || material.poisson >= 0.5 ||
        !std::isfinite(material.j_min) || material.j_min < 0.0) {
        throw std::invalid_argument("Invalid implicit-solver material");
    }
    for (const double value : gravity_) {
        if (!std::isfinite(value)) throw std::invalid_argument("Gravity must be finite");
    }
    if (!std::isfinite(options_.residual_tolerance) || options_.residual_tolerance <= 0.0 ||
        options_.line_search_max_iterations == 0 ||
        !std::isfinite(options_.backtrack_factor) || options_.backtrack_factor <= 0.0 ||
        options_.backtrack_factor >= 1.0 || !std::isfinite(options_.armijo) ||
        options_.armijo <= 0.0 || options_.armijo >= 1.0 ||
        !std::isfinite(options_.linear_relative_tolerance) || options_.linear_relative_tolerance <= 0.0 ||
        options_.diagonal_shifts.empty()) {
        throw std::invalid_argument("Invalid Newton or line-search options");
    }
    double previous = -1.0;
    for (const double value : options_.diagonal_shifts) {
        if (!std::isfinite(value) || value < 0.0 || value <= previous) {
            throw std::invalid_argument("Diagonal shifts must be finite, nonnegative and increasing");
        }
        previous = value;
    }
    std::vector<double> acceleration(mass_.rows());
    for (std::size_t i = 0; i < acceleration.size(); ++i) acceleration[i] = gravity_[i % 3];
    force_ = mass_.multiply(acceleration);
}

double ImplicitEulerSolver::incremental_potential(const SimulationState& state,
                                                 const std::vector<Vec3>& candidate,
                                                 double dt,
                                                 const std::vector<std::size_t>& fixed_vertices) const {
    const double dt2 = validate_step(state, geometry_.reference_positions.size(), dt);
    validate_points(candidate, state.positions.size(), "Candidate positions have wrong size or non-finite values");
    const auto y = prediction(flatten(state.positions), flatten(state.velocities), dt);
    const auto x = flatten(candidate);
    const auto contact = evaluate_ground_penalty(candidate, lumped_masses_, fixed_vertices, dt, ground_);
    return objective(mass_, force_, y, x, elastic_first_order(geometry_, material_, x).energy,
                     contact.energy, dt2);
}

StepResult ImplicitEulerSolver::step(const SimulationState& state,
                                    const std::vector<std::size_t>& fixed_vertices,
                                    double dt) const {
    const double dt2 = validate_step(state, geometry_.reference_positions.size(), dt);
    const DofMap map(state.positions.size(), fixed_vertices);
    const auto x0 = flatten(state.positions);
    const auto v0 = flatten(state.velocities);
    StepResult result{state, {}};
    auto& info = result.info;
    info.contact.penetration.assign(state.positions.size(), 0.0);
    info.contact.normal_impulses.assign(state.positions.size(), 0.0);
    info.contact.tangential_impulses.assign(state.positions.size(), 0.0);
    info.contact.friction_dissipation_per_node.assign(state.positions.size(), 0.0);
    info.contact.position_corrections.assign(state.positions.size(), Vec3{0.0, 0.0, 0.0});
    auto fail = [&](StepFailure reason, const std::string& message) {
        info.failure_reason = reason;
        info.message = message;
        // result.state remains the original input, including its velocities.
        return result;
    };

    try {
        // Reject an already collapsed input instead of hiding it with a predictor.
        try { (void)elastic_first_order(geometry_, material_, x0); }
        catch (const std::domain_error& error) {
            return fail(StepFailure::invalid_initial_configuration, error.what());
        }
        const auto y = prediction(x0, v0, dt);
        auto x = y;
        for (std::size_t i = 0; i < x.size(); ++i) x[i] += dt2 * gravity_[i % 3];
        for (const auto vertex : fixed_vertices) {
            for (std::size_t axis = 0; axis < 3; ++axis) x[3 * vertex + axis] = x0[3 * vertex + axis];
        }
        FirstOrder elastic;
        if (!finite_vector(x)) {
            info.start_fallback = true;
            x = x0;
        }
        try { elastic = elastic_first_order(geometry_, material_, x); }
        catch (const std::domain_error&) {
            info.start_fallback = true;
            x = x0;
            elastic = elastic_first_order(geometry_, material_, x);
        }
        auto contact = evaluate_ground_penalty(unflatten(x), lumped_masses_, fixed_vertices, dt, ground_);
        double phi = objective(mass_, force_, y, x, elastic.energy, contact.energy, dt2);
        info.accepted_objectives.push_back(phi);

        while (true) {
            auto displacement = x;
            for (std::size_t i = 0; i < x.size(); ++i) displacement[i] -= y[i];
            if (!finite_vector(displacement)) {
                return fail(StepFailure::non_finite_arithmetic, "Non-finite inertial displacement");
            }
            const auto inertia = mass_.multiply(displacement);
            std::vector<double> residual(x.size()), scaled_internal(x.size()), scaled_force(x.size());
            for (std::size_t i = 0; i < x.size(); ++i) {
                const double ground_gradient = i % 3 == 1 ? contact.gradient_y[i / 3] : 0.0;
                scaled_internal[i] = dt2 * (elastic.gradient[i] + ground_gradient);
                scaled_force[i] = dt2 * force_[i];
                residual[i] = inertia[i] + scaled_internal[i] - scaled_force[i];
            }
            if (!finite_vector(residual) || !finite_vector(scaled_internal) || !finite_vector(scaled_force)) {
                return fail(StepFailure::non_finite_arithmetic, "Non-finite implicit residual");
            }
            const auto residual_free = restrict_to_free(residual, map);
            info.residual_norm = norm(residual_free);
            const double scale = std::max({1.0, norm(restrict_to_free(inertia, map)),
                norm(restrict_to_free(scaled_force, map)), norm(restrict_to_free(scaled_internal, map))});
            info.residual_threshold = options_.residual_tolerance * scale;
            info.final_objective = phi;
            if (!std::isfinite(info.residual_threshold)) {
                return fail(StepFailure::non_finite_arithmetic, "Non-finite residual threshold");
            }
            if (info.residual_norm <= info.residual_threshold) {
                auto velocity = x;
                for (std::size_t i = 0; i < x.size(); ++i) velocity[i] = (x[i] - x0[i]) / dt;
                for (const auto vertex : fixed_vertices) {
                    for (std::size_t axis = 0; axis < 3; ++axis) velocity[3 * vertex + axis] = 0.0;
                }
                if (!finite_vector(velocity)) {
                    return fail(StepFailure::non_finite_arithmetic, "Non-finite updated velocity");
                }
                auto final_velocity = unflatten(velocity);
                if (ground_.enabled) {
                    auto friction = apply_coulomb_friction(final_velocity, contact.normal_impulses,
                                                           lumped_masses_, fixed_vertices, ground_.friction);
                    auto& contact_info = friction.info;
                    contact_info.penetration = contact.penetration;
                    contact_info.penalty_energy = contact.energy;
                    // No second normal impulse or position projection occurs.
                    // "Before/after" bound this friction-only stage.
                    contact_info.active_count = 0;
                    for (std::size_t i = 0; i < state.positions.size(); ++i) {
                        if (map.full_to_free()[3 * i] == DofMap::npos) {
                            const long double depth = std::max(0.0L,
                                static_cast<long double>(ground_.ground_y) - x[3 * i + 1]);
                            const double penetration = static_cast<double>(depth);
                            if (!std::isfinite(penetration)) throw std::overflow_error("Fixed ground penetration overflows");
                            if (penetration > 0.0) ++contact_info.fixed_penetration_count;
                            contact_info.max_fixed_penetration = std::max(contact_info.max_fixed_penetration, penetration);
                        } else {
                            if (contact.penetration[i] > 0.0) ++contact_info.active_count;
                            contact_info.max_penetration_before = std::max(contact_info.max_penetration_before,
                                                                          contact.penetration[i]);
                        }
                    }
                    contact_info.max_penetration_after = contact_info.max_penetration_before;
                    const auto after = flatten(friction.velocities);
                    std::vector<double> delta(velocity.size()), average(velocity.size());
                    for (std::size_t i = 0; i < velocity.size(); ++i) {
                        delta[i] = after[i] - velocity[i];
                        average[i] = 0.5 * after[i] + 0.5 * velocity[i];
                    }
                    if (!finite_vector(delta) || !finite_vector(average)) {
                        throw std::overflow_error("Non-finite friction kinetic energy arguments");
                    }
                    // delta^T M average equals K(after)-K(before). With a
                    // consistent M it need not be -the nodal dissipation,
                    // and the split nodal friction need not dissipate that K.
                    contact_info.friction_kinetic_change = checked_dot(delta, mass_.multiply(average));
                    info.contact = std::move(contact_info);
                    final_velocity = std::move(friction.velocities);
                }
                result.state = {unflatten(x), std::move(final_velocity)};
                info.solver_converged = true;
                info.committed = true;
                return result;
            }
            if (info.iterations == options_.max_iterations) {
                return fail(StepFailure::iteration_limit, "Newton update budget exhausted");
            }
            ++info.iterations;
            auto working_mesh = geometry_;
            working_mesh.positions = unflatten(x);
            const auto assembled = assemble_elastic(working_mesh, material_);
            const auto system = newton_system(mass_, assembled.hessian, contact.diagonal_y, map, dt2);
            double diagonal_scale = 1.0;
            for (std::size_t i = 0; i < system.rows(); ++i) {
                diagonal_scale = std::max(diagonal_scale, std::abs(system.coeff(i, i)));
            }
            auto rhs = residual_free;
            for (auto& value : rhs) value = -value;
            bool solved_any = false;
            bool accepted = false;
            for (const double factor : options_.diagonal_shifts) {
                const double diagonal_shift = factor * diagonal_scale;
                if (!std::isfinite(diagonal_shift)) continue;
                const auto solve = solve_sparse_lu(shifted(system, diagonal_shift), rhs,
                                                    options_.linear_relative_tolerance);
                if (!solve.success || solve.solution.size() != rhs.size() || !finite_vector(solve.solution)) continue;
                solved_any = true;
                const double slope = checked_dot(residual_free, solve.solution);
                if (!(slope < 0.0)) continue;
                const auto direction = lift_from_free(solve.solution, map);
                double alpha = 1.0;
                for (std::size_t trial_index = 0; trial_index < options_.line_search_max_iterations; ++trial_index) {
                    auto trial = x;
                    for (std::size_t i = 0; i < x.size(); ++i) trial[i] += alpha * direction[i];
                    bool valid = finite_vector(trial);
                    double trial_phi = std::numeric_limits<double>::infinity();
                    FirstOrder trial_elastic;
                    GroundPenaltyEvaluation trial_contact;
                    if (valid) {
                        try {
                            trial_elastic = elastic_first_order(geometry_, material_, trial);
                            trial_contact = evaluate_ground_penalty(unflatten(trial), lumped_masses_,
                                                                    fixed_vertices, dt, ground_);
                            trial_phi = objective(mass_, force_, y, trial, trial_elastic.energy,
                                                  trial_contact.energy, dt2);
                        } catch (const std::domain_error&) {
                            ++info.degenerate_rejections;
                            valid = false;
                        } catch (const std::overflow_error&) {
                            valid = false;
                        }
                    }
                    // Strict decrease prevents a rounded-to-identical trial
                    // from being mistaken for progress at a large objective.
                    if (valid && trial_phi < phi &&
                        trial_phi <= phi + options_.armijo * alpha * slope) {
                        x = std::move(trial);
                        elastic = std::move(trial_elastic);
                        contact = std::move(trial_contact);
                        phi = trial_phi;
                        info.accepted_objectives.push_back(phi);
                        ++info.accepted_steps;
                        if (factor > 0.0) ++info.shifted_steps;
                        accepted = true;
                        break;
                    }
                    ++info.line_search_rejections;
                    alpha *= options_.backtrack_factor;
                }
                if (accepted) break;
            }
            if (!accepted) {
                return fail(solved_any ? StepFailure::line_search_failed : StepFailure::linear_solve_failed,
                    solved_any ? "No finite descent direction passed line search" : "All shifted sparse solves failed");
            }
        }
    } catch (const std::overflow_error& error) {
        return fail(StepFailure::non_finite_arithmetic, error.what());
    } catch (const std::domain_error& error) {
        return fail(StepFailure::invalid_initial_configuration, error.what());
    }
}

} // namespace cloth
