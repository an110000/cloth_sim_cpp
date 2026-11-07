#pragma once

#include "cloth/state.hpp"

#include <cstddef>
#include <optional>
#include <vector>

namespace cloth {

struct FrictionOptions {
    double kinetic = 0.3;
    std::optional<double> static_coefficient = std::nullopt;
};

struct GroundContactOptions {
    bool enabled = false;
    double ground_y = 0.0;
    FrictionOptions friction{};
    double restitution = 0.0;
    double stabilization_beta = 0.0;
    double proximity_tolerance = 1.0e-12;
    double penalty_eta = 1.0e3;
};

void validate_ground_contact_options(const GroundContactOptions& options);

struct ContactInfo {
    std::size_t active_count = 0;
    std::size_t sticking_count = 0;
    std::size_t sliding_count = 0;
    std::size_t fixed_penetration_count = 0;
    double max_penetration_before = 0.0;
    double max_penetration_after = 0.0;
    double max_fixed_penetration = 0.0;
    double normal_impulse_total = 0.0;
    double tangential_impulse_total = 0.0;
    double friction_dissipation = 0.0;
    double friction_kinetic_change = 0.0;
    double normal_kinetic_change = 0.0;
    double penalty_energy = 0.0;
    // Per-node values exclude fixed vertices. Tangential impulses are magnitudes.
    std::vector<double> penetration;
    std::vector<double> normal_impulses;
    std::vector<double> tangential_impulses;
    std::vector<double> friction_dissipation_per_node;
    std::vector<Vec3> position_corrections;
};

struct FrictionResult {
    std::vector<Vec3> velocities;
    ContactInfo info;
};

// The plane normal is +Y. Only the X/Z velocities receive friction impulses.
// Positive finite physical masses are required, including fixed vertices.
FrictionResult apply_coulomb_friction(
    const std::vector<Vec3>& velocities, const std::vector<double>& normal_impulses,
    const std::vector<double>& masses, const std::vector<std::size_t>& fixed_vertices,
    const FrictionOptions& options = {});

struct GroundContactResult {
    SimulationState state;
    ContactInfo info;
};

// Enabled contact projects only free positions and zeros fixed velocities.
// Disabled contact validates the input and returns its state unchanged.
// This value-returning operation never partially modifies its inputs on failure.
GroundContactResult apply_ground_contact(
    const SimulationState& state, const std::vector<double>& masses,
    const std::vector<std::size_t>& fixed_vertices, double dt,
    const GroundContactOptions& options = {});

struct GroundPenaltyEvaluation {
    double energy = 0.0;
    std::vector<double> gradient_y;
    std::vector<double> diagonal_y;
    std::vector<double> normal_impulses;
    std::vector<double> penetration;
};

// E=.5*sum(k_i*pen_i^2), k_i=eta*m_i/dt^2, excluding fixed nodes.
// gradient_y/diagonal_y are derivatives in world Y; Jn=dt*k_i*pen_i.
// At the C1 plane boundary the inactive (zero Hessian) branch is selected.
GroundPenaltyEvaluation evaluate_ground_penalty(
    const std::vector<Vec3>& positions, const std::vector<double>& masses,
    const std::vector<std::size_t>& fixed_vertices, double dt,
    const GroundContactOptions& options = {});

} // namespace cloth
