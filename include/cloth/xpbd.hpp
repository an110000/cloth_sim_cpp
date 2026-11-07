#pragma once

#include "cloth/constraints.hpp"
#include "cloth/contact.hpp"
#include "cloth/state.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace cloth {

struct XpbdOptions {
    std::size_t iterations = 30;
    double distance_tolerance = 1.0e-6; // Compliant residual, metres.
    double angle_tolerance = 1.0e-6;    // Compliant residual, radians.
    // Fixed-budget XPBD normally commits a valid finite approximation. Opt-in
    // strict mode rejects a step whose final compliant residual is too large.
    bool require_convergence = false;
    // Added after the old fields to retain four-field aggregate initializers.
    double energy_tolerance = 1.0e-6; // Compliant residual for dimensionless C=E/scale.
};

struct XpbdGroupMetrics {
    std::size_t count = 0;
    double max_violation = 0.0;          // max |C|; m, rad, or dimensionless by group.
    double max_compliant_residual = 0.0; // max |C + compliance/dt^2 * lambda|.
};

enum class XpbdFailure { none, invalid_geometry, non_finite_arithmetic, projection_not_converged };
const char* failure_name(XpbdFailure reason) noexcept;

struct XpbdStepInfo {
    bool finite_state = true; // Describes the returned state, including rollback.
    bool committed = false;
    // This is the stated per-constraint compliant residual test. It is not a
    // full nonlinear stationarity test or a velocity-level constraint test.
    bool projection_converged = false;
    XpbdFailure failure_reason = XpbdFailure::none;
    std::string message;
    std::size_t iterations = 0; // Completed Gauss-Seidel sweeps.
    std::size_t projections = 0;
    std::size_t invalid_constraints = 0;
    std::array<XpbdGroupMetrics, 4> groups{}; // structural, shear, bending, energy.
    // Contact splitting can disturb an already projected material constraint.
    std::array<XpbdGroupMetrics, 4> pre_contact_groups{};
    bool pre_contact_projection_converged = false;
    ContactInfo contact;
    std::vector<double> multipliers; // Final candidate, even on strict rejection.
};

struct XpbdStepResult { SimulationState state; XpbdStepInfo info; };

class XpbdSolver {
public:
    // Physical nodal masses must be finite and positive. step() converts fixed
    // nodes to zero inverse mass; callers need not encode infinity in mass.
    XpbdSolver(std::vector<double> masses, std::vector<Constraint> constraints,
               Vec3 gravity = {0.0, -9.81, 0.0}, XpbdOptions options = {},
               GroundContactOptions ground = {});

    XpbdStepResult step(const SimulationState& state,
                        const std::vector<std::size_t>& fixed_vertices, double dt) const;

    const std::vector<Constraint>& constraints() const noexcept { return constraints_; }
    const std::vector<double>& masses() const noexcept { return masses_; }
    const XpbdOptions& options() const noexcept { return options_; }

private:
    std::vector<double> masses_;
    std::vector<Constraint> constraints_;
    Vec3 gravity_;
    XpbdOptions options_;
    GroundContactOptions ground_;
};

} // namespace cloth
