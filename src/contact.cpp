#include "cloth/contact.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cloth {
namespace {

constexpr double impulse_tolerance = 1.0e-14;
constexpr double speed_tolerance = 1.0e-14;

double finite_result(long double value, const char* message) {
    if (!std::isfinite(value) ||
        std::abs(value) > static_cast<long double>(std::numeric_limits<double>::max())) {
        throw std::overflow_error(message);
    }
    return static_cast<double>(value);
}

void add_finite(double& total, long double value, const char* message) {
    total = finite_result(static_cast<long double>(total) + value, message);
}

void validate_friction_options(const FrictionOptions& options) {
    const double static_mu = options.static_coefficient.value_or(options.kinetic);
    if (!std::isfinite(options.kinetic) || options.kinetic < 0.0 ||
        !std::isfinite(static_mu) || static_mu < options.kinetic) {
        throw std::invalid_argument("Friction requires finite coefficients with static >= kinetic >= 0");
    }
}

void validate_points(const std::vector<Vec3>& points) {
    for (const auto& point : points) {
        for (const auto value : point) {
            if (!std::isfinite(value)) throw std::invalid_argument("Contact coordinates must be finite");
        }
    }
}

std::vector<bool> validate_mass_and_fixed(std::size_t count, const std::vector<double>& masses,
                                        const std::vector<std::size_t>& fixed_vertices) {
    if (masses.size() != count) throw std::invalid_argument("Contact masses must match the node count");
    for (const auto mass : masses) {
        if (!std::isfinite(mass) || mass <= 0.0) {
            throw std::invalid_argument("Contact requires finite positive physical masses");
        }
    }
    std::vector<bool> fixed(count, false);
    for (const auto node : fixed_vertices) {
        if (node >= count || fixed[node]) throw std::invalid_argument("Invalid or duplicate fixed contact node");
        fixed[node] = true;
    }
    return fixed;
}

void validate_dt(double dt) {
    if (!std::isfinite(dt) || dt <= 0.0 || !std::isfinite(dt * dt) || dt * dt == 0.0) {
        throw std::invalid_argument("Contact time step and its square must be finite and positive");
    }
}

ContactInfo zero_info(std::size_t count) {
    ContactInfo info;
    info.penetration.assign(count, 0.0);
    info.normal_impulses.assign(count, 0.0);
    info.tangential_impulses.assign(count, 0.0);
    info.friction_dissipation_per_node.assign(count, 0.0);
    info.position_corrections.assign(count, Vec3{0.0, 0.0, 0.0});
    return info;
}

double penetration_at(double y, double ground_y) {
    return finite_result(std::max(0.0L, static_cast<long double>(ground_y) - y),
                         "Ground penetration overflows");
}

} // namespace

void validate_ground_contact_options(const GroundContactOptions& options) {
    validate_friction_options(options.friction);
    if (!std::isfinite(options.ground_y) ||
        !std::isfinite(options.restitution) || options.restitution < 0.0 || options.restitution > 1.0 ||
        !std::isfinite(options.stabilization_beta) || options.stabilization_beta < 0.0 ||
        !std::isfinite(options.proximity_tolerance) || options.proximity_tolerance < 0.0 ||
        !std::isfinite(options.penalty_eta) || options.penalty_eta <= 0.0) {
        throw std::invalid_argument("Invalid ground contact options");
    }
}

FrictionResult apply_coulomb_friction(
    const std::vector<Vec3>& velocities, const std::vector<double>& normal_impulses,
    const std::vector<double>& masses, const std::vector<std::size_t>& fixed_vertices,
    const FrictionOptions& options) {
    validate_friction_options(options);
    validate_points(velocities);
    const auto fixed = validate_mass_and_fixed(velocities.size(), masses, fixed_vertices);
    if (normal_impulses.size() != velocities.size()) {
        throw std::invalid_argument("Normal impulses must match the contact node count");
    }
    for (const auto impulse : normal_impulses) {
        if (!std::isfinite(impulse) || impulse < 0.0) {
            throw std::invalid_argument("Normal impulses must be finite and nonnegative");
        }
    }
    FrictionResult result{velocities, zero_info(velocities.size())};
    const double static_mu = options.static_coefficient.value_or(options.kinetic);
    for (std::size_t i = 0; i < velocities.size(); ++i) {
        if (fixed[i]) {
            result.velocities[i] = {0.0, 0.0, 0.0};
            continue;
        }
        auto& info = result.info;
        const double normal = normal_impulses[i];
        info.normal_impulses[i] = normal;
        add_finite(info.normal_impulse_total, normal, "Normal impulse sum overflows");
        if (normal <= impulse_tolerance) continue;
        ++info.active_count;
        const long double vx = velocities[i][0], vz = velocities[i][2];
        const long double speed = std::hypot(vx, vz);
        if (!std::isfinite(speed)) throw std::overflow_error("Tangential speed overflows");
        if (speed <= speed_tolerance) {
            ++info.sticking_count;
            continue;
        }
        const long double mass = masses[i];
        const long double required = mass * speed;
        const long double static_limit = static_cast<long double>(static_mu) * normal;
        if (!std::isfinite(required) || !std::isfinite(static_limit)) {
            throw std::overflow_error("Friction impulse comparison overflows");
        }
        long double tangent = 0.0;
        if (required <= static_limit) {
            tangent = required;
            result.velocities[i][0] = result.velocities[i][2] = 0.0;
            ++info.sticking_count;
        } else {
            tangent = static_cast<long double>(options.kinetic) * normal;
            // static >= kinetic makes the sliding factor nonnegative.
            const long double factor = std::max(0.0L, 1.0L - tangent / required);
            result.velocities[i][0] = finite_result(vx * factor, "Friction velocity overflows");
            result.velocities[i][2] = finite_result(vz * factor, "Friction velocity overflows");
            ++info.sliding_count;
        }
        info.tangential_impulses[i] = finite_result(tangent, "Tangential impulse overflows");
        add_finite(info.tangential_impulse_total, tangent, "Tangential impulse sum overflows");
        const long double after_x = result.velocities[i][0], after_z = result.velocities[i][2];
        // This is the energy difference of the actual returned velocities;
        // the factored form avoids subtracting two large nearly equal squares.
        const long double difference = mass *
            ((vx - after_x) * (0.5L * vx + 0.5L * after_x) +
             (vz - after_z) * (0.5L * vz + 0.5L * after_z));
        if (!std::isfinite(difference)) throw std::overflow_error("Friction dissipation overflows");
        const long double dissipation = std::max(0.0L, difference);
        info.friction_dissipation_per_node[i] = finite_result(dissipation, "Friction dissipation overflows");
        add_finite(info.friction_dissipation, dissipation, "Friction dissipation sum overflows");
    }
    result.info.friction_kinetic_change = -result.info.friction_dissipation;
    return result;
}

GroundContactResult apply_ground_contact(
    const SimulationState& state, const std::vector<double>& masses,
    const std::vector<std::size_t>& fixed_vertices, double dt,
    const GroundContactOptions& options) {
    validate_ground_contact_options(options);
    validate_dt(dt);
    validate_points(state.positions);
    validate_points(state.velocities);
    if (state.velocities.size() != state.positions.size()) {
        throw std::invalid_argument("Contact position and velocity counts must match");
    }
    const auto fixed = validate_mass_and_fixed(state.positions.size(), masses, fixed_vertices);
    GroundContactResult result{state, zero_info(state.positions.size())};
    if (!options.enabled) return result;
    std::size_t active = 0;
    double normal_change = 0.0;
    for (std::size_t i = 0; i < state.positions.size(); ++i) {
        auto& info = result.info;
        const double penetration = penetration_at(state.positions[i][1], options.ground_y);
        if (fixed[i]) {
            result.state.velocities[i] = {0.0, 0.0, 0.0};
            if (penetration > 0.0) ++info.fixed_penetration_count;
            info.max_fixed_penetration = std::max(info.max_fixed_penetration, penetration);
            continue;
        }
        info.penetration[i] = penetration;
        info.max_penetration_before = std::max(info.max_penetration_before, penetration);
        if (penetration > 0.0) {
            result.state.positions[i][1] = options.ground_y;
            info.position_corrections[i][1] = penetration;
        }
        info.max_penetration_after = std::max(info.max_penetration_after,
            penetration_at(result.state.positions[i][1], options.ground_y));
        const long double distance = std::abs(static_cast<long double>(state.positions[i][1]) - options.ground_y);
        if (penetration == 0.0 && distance > options.proximity_tolerance) continue;
        ++active;
        const long double before = state.velocities[i][1];
        const long double target = -static_cast<long double>(options.restitution) * std::min(before, 0.0L) +
            static_cast<long double>(options.stabilization_beta) * penetration / dt;
        if (!std::isfinite(target)) throw std::overflow_error("Normal target velocity overflows");
        if (target <= before) continue; // Preserve a body leaving the plane.
        const double after = finite_result(target, "Normal velocity overflows");
        result.state.velocities[i][1] = after;
        // Report the impulse of the actual rounded velocity change.
        const long double impulse = static_cast<long double>(masses[i]) * (static_cast<long double>(after) - before);
        info.normal_impulses[i] = finite_result(impulse, "Normal impulse overflows");
        const long double change = impulse * (0.5L * before + 0.5L * after);
        add_finite(normal_change, change, "Normal kinetic energy change overflows");
    }
    const auto friction = apply_coulomb_friction(result.state.velocities, result.info.normal_impulses,
                                                masses, fixed_vertices, options.friction);
    result.state.velocities = friction.velocities;
    result.info.active_count = active;
    result.info.sticking_count = friction.info.sticking_count;
    result.info.sliding_count = friction.info.sliding_count;
    result.info.normal_impulse_total = friction.info.normal_impulse_total;
    result.info.tangential_impulse_total = friction.info.tangential_impulse_total;
    result.info.friction_dissipation = friction.info.friction_dissipation;
    result.info.friction_kinetic_change = friction.info.friction_kinetic_change;
    result.info.normal_kinetic_change = normal_change;
    result.info.tangential_impulses = friction.info.tangential_impulses;
    result.info.friction_dissipation_per_node = friction.info.friction_dissipation_per_node;
    return result;
}

GroundPenaltyEvaluation evaluate_ground_penalty(
    const std::vector<Vec3>& positions, const std::vector<double>& masses,
    const std::vector<std::size_t>& fixed_vertices, double dt,
    const GroundContactOptions& options) {
    validate_ground_contact_options(options);
    validate_dt(dt);
    validate_points(positions);
    const auto fixed = validate_mass_and_fixed(positions.size(), masses, fixed_vertices);
    GroundPenaltyEvaluation result;
    result.gradient_y.assign(positions.size(), 0.0);
    result.diagonal_y.assign(positions.size(), 0.0);
    result.normal_impulses.assign(positions.size(), 0.0);
    result.penetration.assign(positions.size(), 0.0);
    if (!options.enabled) return result;
    const long double dt2 = static_cast<long double>(dt) * dt;
    for (std::size_t i = 0; i < positions.size(); ++i) {
        if (fixed[i]) continue;
        const double penetration = penetration_at(positions[i][1], options.ground_y);
        result.penetration[i] = penetration;
        if (penetration == 0.0) continue;
        const long double stiffness = static_cast<long double>(options.penalty_eta) * masses[i] / dt2;
        result.diagonal_y[i] = finite_result(stiffness, "Ground penalty stiffness overflows");
        const long double force = stiffness * penetration;
        result.gradient_y[i] = finite_result(-force, "Ground penalty gradient overflows");
        result.normal_impulses[i] = finite_result(dt * force, "Ground penalty impulse overflows");
        add_finite(result.energy, 0.5L * force * penetration, "Ground penalty energy overflows");
    }
    return result;
}

} // namespace cloth
