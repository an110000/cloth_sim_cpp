#include "cloth/xpbd.hpp"

#include "cloth/dofs.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace cloth {
namespace {

void validate_points(const std::vector<Vec3>& points, std::size_t count) {
    if (points.size() != count) throw std::invalid_argument("XPBD state size does not match masses");
    for (const auto& point : points) {
        for (const auto value : point) {
            if (!std::isfinite(value)) throw std::invalid_argument("XPBD state must be finite");
        }
    }
}

bool finite_points(const std::vector<Vec3>& points) {
    for (const auto& point : points)
        for (const auto value : point)
            if (!std::isfinite(value)) return false;
    return true;
}

void measure(const std::vector<Constraint>& constraints, const std::vector<Vec3>& positions,
             double dt2, const XpbdOptions& options, XpbdStepInfo& info) {
    info.groups = {};
    info.invalid_constraints = 0;
    info.multipliers.clear();
    info.multipliers.reserve(constraints.size());
    for (const auto& constraint : constraints) {
        const auto group = constraint_group(constraint);
        auto& metrics = info.groups[static_cast<std::size_t>(group)];
        ++metrics.count;
        const auto evaluation = evaluate_constraint(constraint, positions);
        const auto multiplier = constraint_multiplier(constraint);
        info.multipliers.push_back(multiplier);
        if (!evaluation.valid_geometry) {
            ++info.invalid_constraints;
            metrics.max_violation = metrics.max_compliant_residual = std::numeric_limits<double>::infinity();
            continue;
        }
        const long double value = static_cast<long double>(evaluation.value) +
            static_cast<long double>(constraint_compliance(constraint)) / dt2 * multiplier;
        const double residual = static_cast<double>(std::abs(value));
        if (!std::isfinite(residual)) throw std::overflow_error("XPBD compliant residual is non-finite");
        metrics.max_violation = std::max(metrics.max_violation, std::abs(evaluation.value));
        metrics.max_compliant_residual = std::max(metrics.max_compliant_residual, residual);
    }
    info.projection_converged = info.invalid_constraints == 0 &&
        info.groups[0].max_compliant_residual <= options.distance_tolerance &&
        info.groups[1].max_compliant_residual <= options.distance_tolerance &&
        info.groups[2].max_compliant_residual <= options.angle_tolerance &&
        info.groups[3].max_compliant_residual <= options.energy_tolerance;
}

} // namespace

const char* failure_name(XpbdFailure reason) noexcept {
    switch (reason) {
    case XpbdFailure::none: return "none";
    case XpbdFailure::invalid_geometry: return "invalid_geometry";
    case XpbdFailure::non_finite_arithmetic: return "non_finite_arithmetic";
    case XpbdFailure::projection_not_converged: return "projection_not_converged";
    }
    return "unknown";
}

XpbdSolver::XpbdSolver(std::vector<double> masses, std::vector<Constraint> constraints,
                       Vec3 gravity, XpbdOptions options, GroundContactOptions ground)
    : masses_(std::move(masses)), constraints_(std::move(constraints)), gravity_(gravity),
      options_(options), ground_(ground) {
    validate_ground_contact_options(ground_);
    if (masses_.empty() || options_.iterations == 0 ||
        !std::isfinite(options_.distance_tolerance) || options_.distance_tolerance <= 0.0 ||
        !std::isfinite(options_.angle_tolerance) || options_.angle_tolerance <= 0.0 ||
        !std::isfinite(options_.energy_tolerance) || options_.energy_tolerance <= 0.0) {
        throw std::invalid_argument("Invalid XPBD masses, iteration budget or tolerances");
    }
    for (const auto mass : masses_) {
        if (!std::isfinite(mass) || mass <= 0.0 || !std::isfinite(1.0 / mass) || 1.0 / mass == 0.0) {
            throw std::invalid_argument("XPBD masses require finite positive inverse masses");
        }
    }
    for (const auto value : gravity_) {
        if (!std::isfinite(value)) throw std::invalid_argument("XPBD gravity must be finite");
    }
    for (auto& constraint : constraints_) {
        validate_constraint(constraint, masses_.size());
        reset_constraint(constraint);
    }
}

XpbdStepResult XpbdSolver::step(const SimulationState& state,
                               const std::vector<std::size_t>& fixed_vertices, double dt) const {
    validate_points(state.positions, masses_.size());
    validate_points(state.velocities, masses_.size());
    if (!std::isfinite(dt) || dt <= 0.0 || !std::isfinite(dt * dt) || dt * dt == 0.0) {
        throw std::invalid_argument("XPBD time step and its square must be finite and positive");
    }
    const DofMap map(masses_.size(), fixed_vertices); // Validates duplicates and bounds.
    (void)map;
    const double dt2 = dt * dt;
    XpbdStepResult result{state, {}};
    auto& info = result.info;
    auto constraints = constraints_;
    for (auto& constraint : constraints) reset_constraint(constraint);
    auto fail = [&](XpbdFailure reason, const char* message) {
        info.failure_reason = reason;
        info.projection_converged = false;
        info.message = message;
        return result;
    };
    try {
        measure(constraints, state.positions, dt2, options_, info);
        if (info.invalid_constraints != 0) return fail(XpbdFailure::invalid_geometry, "Invalid initial constraint geometry");
        auto candidate = state;
        std::vector<double> inverse_masses(masses_.size());
        std::vector<bool> fixed(masses_.size(), false);
        for (std::size_t i = 0; i < masses_.size(); ++i) inverse_masses[i] = 1.0 / masses_[i];
        for (const auto vertex : fixed_vertices) { inverse_masses[vertex] = 0.0; fixed[vertex] = true; }
        for (std::size_t vertex = 0; vertex < masses_.size(); ++vertex) {
            for (std::size_t axis = 0; axis < 3; ++axis) {
                candidate.velocities[vertex][axis] = fixed[vertex] ? 0.0 :
                    state.velocities[vertex][axis] + dt * gravity_[axis];
                candidate.positions[vertex][axis] = state.positions[vertex][axis] +
                    dt * candidate.velocities[vertex][axis];
            }
        }
        if (!finite_points(candidate.positions) || !finite_points(candidate.velocities)) {
            return fail(XpbdFailure::non_finite_arithmetic, "XPBD prediction overflow");
        }
        for (std::size_t sweep = 0; sweep < options_.iterations; ++sweep) {
            for (auto& constraint : constraints) {
                const auto projection = project_constraint(constraint, candidate.positions, inverse_masses, dt);
                ++info.projections;
                if (!projection.valid_geometry) {
                    measure(constraints, candidate.positions, dt2, options_, info);
                    return fail(XpbdFailure::invalid_geometry, "Invalid XPBD projection geometry");
                }
            }
            ++info.iterations;
        }
        measure(constraints, candidate.positions, dt2, options_, info);
        if (info.invalid_constraints != 0) return fail(XpbdFailure::invalid_geometry, "Invalid final constraint geometry");
        info.pre_contact_groups = info.groups;
        info.pre_contact_projection_converged = info.projection_converged;
        for (std::size_t vertex = 0; vertex < masses_.size(); ++vertex) {
            for (std::size_t axis = 0; axis < 3; ++axis) {
                candidate.velocities[vertex][axis] = fixed[vertex] ? 0.0 :
                    (candidate.positions[vertex][axis] - state.positions[vertex][axis]) / dt;
            }
        }
        if (!finite_points(candidate.positions) || !finite_points(candidate.velocities)) {
            return fail(XpbdFailure::non_finite_arithmetic, "XPBD corrected state is not finite");
        }
        if (ground_.enabled) {
            auto contact = apply_ground_contact(candidate, masses_, fixed_vertices, dt, ground_);
            candidate = std::move(contact.state);
            info.contact = std::move(contact.info);
            measure(constraints, candidate.positions, dt2, options_, info);
            if (info.invalid_constraints != 0)
                return fail(XpbdFailure::invalid_geometry, "Ground projection invalidated material geometry");
        }
        if (options_.require_convergence && !info.projection_converged) {
            return fail(XpbdFailure::projection_not_converged, "Final post-contact compliant residual exceeds the stated tolerance");
        }
        result.state = std::move(candidate);
        info.committed = true;
        info.message = info.projection_converged ? "compliant residual test passed" : "finite fixed-budget approximation";
        return result;
    } catch (const std::overflow_error& error) {
        info.failure_reason = XpbdFailure::non_finite_arithmetic;
        info.projection_converged = false;
        info.message = error.what();
        return result;
    }
}

} // namespace cloth
