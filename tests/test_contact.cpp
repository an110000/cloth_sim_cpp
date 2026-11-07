#include "cloth/contact.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

void check(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}

void near(double actual, double expected, double tolerance, const char* message) {
    check(std::isfinite(actual) && std::isfinite(expected) &&
          std::abs(actual - expected) <= tolerance * std::max(1.0, std::abs(expected)), message);
}

template<class Exception, class Function>
void rejects(Function function, const char* message) {
    try { function(); } catch (const Exception&) { return; }
    throw std::runtime_error(message);
}

cloth::GroundContactOptions ground() {
    cloth::GroundContactOptions options;
    options.enabled = true;
    return options;
}

double tangent_energy(const cloth::Vec3& velocity, double mass) {
    return 0.5 * mass * (velocity[0] * velocity[0] + velocity[2] * velocity[2]);
}

void friction_hand_calculation() {
    const std::vector<cloth::Vec3> velocities{{3, 7, 4}, {.1, -2, 0}, {0, 9, 0}, {1, 2, 1}};
    const std::vector<double> impulses{2, 1, 2, 100}, masses{2, 2, 1, 1};
    const auto original = velocities;
    const auto result = cloth::apply_coulomb_friction(velocities, impulses, masses, {3}, {.3, .5});
    check(velocities == original, "Friction has const input state");
    near(result.velocities[0][0], 2.82, 1e-14, "Sliding X velocity");
    near(result.velocities[0][2], 3.76, 1e-14, "Sliding Z velocity");
    near(result.velocities[0][1], 7, 0, "Friction preserves normal velocity");
    near(result.info.tangential_impulses[0], .6, 1e-14, "Sliding impulse is mu_k Jn");
    near(result.info.friction_dissipation_per_node[0], 2.91, 1e-14, "Sliding dissipation is actual kinetic loss");
    check(result.velocities[1] == cloth::Vec3{0, -2, 0}, "Static friction stops tangent velocity");
    near(result.info.tangential_impulses[1], .2, 1e-14, "Static impulse may exceed kinetic limit");
    near(result.info.friction_dissipation_per_node[1], .01, 1e-14, "Sticking kinetic energy loss");
    check(result.info.active_count == 3 && result.info.sticking_count == 2 &&
          result.info.sliding_count == 1, "Positive-normal zero-tangent node counts as sticking");
    check(result.velocities[3] == cloth::Vec3{0, 0, 0}, "Fixed friction velocity is zero");
    near(result.info.normal_impulses[3], 0, 0, "Fixed normal impulse excluded");
    near(result.info.normal_impulse_total, 5, 1e-14, "Free normal impulse sum");
    near(result.info.tangential_impulse_total, .8, 1e-14, "Tangential impulse sum");
    near(result.info.friction_dissipation, 2.92, 1e-14, "Total friction dissipation");
    near(result.info.friction_kinetic_change, -2.92, 1e-14, "Friction kinetic change has correct sign");
    const auto tiny = cloth::apply_coulomb_friction({{1, 2, 3}}, {1e-15}, {1}, {}, {.3, .5});
    check(tiny.velocities[0] == cloth::Vec3{1, 2, 3} && tiny.info.active_count == 0,
          "Original small-impulse protection is explicit");
}

void friction_physical_properties() {
    for (const double mu : {0.0, .1, .7}) {
        for (const double normal : {0.0, .2, 4.0}) {
            for (const double mass : {.2, 1.0, 3.0}) {
                for (const double speed : {0.0, .01, .1, 1.0, 10.0}) {
                    const cloth::Vec3 before{-.6 * speed, -7, .8 * speed};
                    const auto result = cloth::apply_coulomb_friction({before}, {normal}, {mass}, {}, {mu, mu + .2});
                    const auto after = result.velocities[0];
                    near(after[1], before[1], 0, "Coulomb step never changes Y");
                    check(before[0] * after[0] >= 0 && before[2] * after[2] >= 0,
                          "Coulomb step never reverses tangential motion");
                    const bool sticking = result.info.sticking_count != 0;
                    const double limit = (sticking ? mu + .2 : mu) * normal;
                    check(result.info.tangential_impulses[0] <= limit + 1e-13,
                          "Stick/slip impulse satisfies its friction cone");
                    const double loss = tangent_energy(before, mass) - tangent_energy(after, mass);
                    check(loss >= -1e-13, "Coulomb friction does not increase tangential energy");
                    near(result.info.friction_dissipation, loss, 2e-13, "Dissipation equals returned kinetic energy difference");
                    near(result.info.friction_kinetic_change, -loss, 2e-13, "Signed friction energy ledger");
                    if (normal > 0 && speed > 0 && !sticking) {
                        near(result.info.tangential_impulses[0], mu * normal, 1e-14, "Sliding lies on kinetic cone");
                    }
                }
            }
        }
    }
}

void ground_impulses_and_energy() {
    auto options = ground();
    options.restitution = .5;
    options.friction.kinetic = .25;
    const cloth::SimulationState initial{{{0, -.1, 0}}, {{1, -2, 0}}};
    const auto result = cloth::apply_ground_contact(initial, {2}, {}, .1, options);
    check(initial.positions[0][1] == -.1 && initial.velocities[0][1] == -2, "Ground inputs unchanged");
    near(result.state.positions[0][1], 0, 0, "Free penetration projected to plane");
    near(result.info.position_corrections[0][1], .1, 1e-14, "Position correction recorded separately");
    near(result.info.max_penetration_before, .1, 1e-14, "Free penetration before");
    near(result.info.max_penetration_after, 0, 0, "Free penetration after");
    near(result.state.velocities[0][1], 1, 1e-14, "Restitution target");
    near(result.info.normal_impulses[0], 6, 1e-14, "Normal impulse is m times actual normal velocity change");
    near(result.state.velocities[0][0], .25, 1e-14, "Contact uses normal impulse for friction");
    near(result.info.tangential_impulse_total, 1.5, 1e-14, "Contact kinetic friction impulse");
    near(result.info.normal_kinetic_change, -3, 1e-14, "Normal restitution energy loss");
    near(result.info.friction_dissipation, .9375, 1e-14, "Contact tangent energy loss");
    check(result.info.active_count == 1 && result.info.sliding_count == 1, "Ground contact counts");

    const cloth::SimulationState separating{{{0, -.1, 0}}, {{2, 3, 0}}};
    const auto leaving = cloth::apply_ground_contact(separating, {2}, {}, .1, options);
    check(leaving.state.velocities == separating.velocities, "Separating normal velocity is preserved");
    near(leaving.info.normal_impulse_total, 0, 0, "No normal impulse for separating contact");
    near(leaving.info.friction_dissipation, 0, 0, "No friction without a normal impulse");

    options = ground();
    options.stabilization_beta = .2;
    options.friction.kinetic = .25;
    const cloth::SimulationState resting_penetration{{{0, -.1, 0}}, {{1, 0, 0}}};
    const auto stabilized = cloth::apply_ground_contact(resting_penetration, {2}, {}, .1, options);
    near(stabilized.state.velocities[0][1], .2, 1e-14, "Positive beta changes the actual normal velocity");
    near(stabilized.info.normal_impulse_total, .4, 1e-14, "Stabilization impulse agrees with actual m delta-v");
    near(stabilized.info.normal_kinetic_change, .04, 1e-14, "Stabilization energy injection explicitly recorded");
    near(stabilized.info.friction_dissipation, .0975, 1e-14, "Stabilization-driven friction uses actual impulse");

    options = ground();
    options.ground_y = 1000;
    const cloth::SimulationState distant{{{0, 1000 + 1e-5, 0}}, {{0, -2, 0}}};
    const auto no_relative_tolerance = cloth::apply_ground_contact(distant, {1}, {}, .1, options);
    check(no_relative_tolerance.info.active_count == 0 && no_relative_tolerance.state.velocities == distant.velocities,
          "Plane proximity is absolute and does not inherit numpy relative tolerance");
}

void fixed_disabled_and_empty() {
    const cloth::SimulationState input{{{0, -.1, 0}, {1, -.4, 0}}, {{0, -2, 1}, {4, -3, 2}}};
    const auto result = cloth::apply_ground_contact(input, {2, 3}, {1}, .1, ground());
    check(result.state.positions[1] == input.positions[1], "Contact never projects fixed positions");
    check(result.state.velocities[1] == cloth::Vec3{0, 0, 0}, "Enabled contact zeros fixed velocity");
    check(result.info.fixed_penetration_count == 1 && result.info.active_count == 1, "Fixed penetration reported separately");
    near(result.info.max_fixed_penetration, .4, 1e-14, "Fixed penetration maximum");
    near(result.info.max_penetration_before, .1, 1e-14, "Fixed penetration excluded from free maximum");
    near(result.info.penetration[1], 0, 0, "Fixed per-node free penetration excluded");
    near(result.info.normal_impulses[1], 0, 0, "Fixed per-node normal impulse excluded");
    check(result.info.position_corrections[1] == cloth::Vec3{0, 0, 0}, "Fixed position correction is zero");
    const cloth::GroundContactOptions disabled;
    check(disabled.stabilization_beta == 0, "Default stabilization does not inject normal energy");
    const auto unchanged = cloth::apply_ground_contact(input, {2, 3}, {1}, .1, disabled);
    check(unchanged.state.positions == input.positions && unchanged.state.velocities == input.velocities,
          "Disabled ground returns original state including fixed velocities");
    check(unchanged.info.active_count == 0 && unchanged.info.fixed_penetration_count == 0,
          "Disabled ground reports no contact");
    const auto empty = cloth::apply_ground_contact({}, {}, {}, .1, ground());
    check(empty.state.positions.empty() && empty.info.penetration.empty(), "Empty valid contact batch");
}

void penalty_derivatives_and_scaling() {
    auto options = ground();
    options.penalty_eta = 10;
    const std::vector<cloth::Vec3> positions{{0, -.2, 0}, {0, 0, 0}, {0, .3, 0}, {0, -.4, 0}};
    const std::vector<double> masses{2, 3, 4, 5};
    const auto result = cloth::evaluate_ground_penalty(positions, masses, {3}, .1, options);
    near(result.energy, 40, 1e-14, "Ground penalty energy");
    near(result.gradient_y[0], -400, 1e-14, "Penalty gradient points down, force up");
    near(result.diagonal_y[0], 2000, 1e-14, "Penalty stiffness eta m / dt squared");
    near(result.normal_impulses[0], 40, 1e-14, "Penalty impulse dt times upward force");
    near(result.penetration[0], .2, 1e-14, "Penalty free penetration");
    for (const std::size_t i : {std::size_t{1}, std::size_t{2}, std::size_t{3}}) {
        near(result.gradient_y[i], 0, 0, "Inactive or fixed penalty gradient zero");
        near(result.diagonal_y[i], 0, 0, "Inactive or fixed penalty Hessian zero");
        near(result.normal_impulses[i], 0, 0, "Inactive or fixed penalty impulse zero");
    }
    auto plus = positions, minus = positions;
    const double epsilon = 1e-6;
    plus[0][1] += epsilon;
    minus[0][1] -= epsilon;
    const auto positive = cloth::evaluate_ground_penalty(plus, masses, {3}, .1, options);
    const auto negative = cloth::evaluate_ground_penalty(minus, masses, {3}, .1, options);
    near((positive.energy - negative.energy) / (2 * epsilon), result.gradient_y[0], 2e-10,
         "Penalty energy finite difference");
    near((positive.gradient_y[0] - negative.gradient_y[0]) / (2 * epsilon), result.diagonal_y[0], 2e-10,
         "Penalty gradient finite difference Hessian");
    const auto half_dt = cloth::evaluate_ground_penalty(positions, masses, {3}, .05, options);
    near(half_dt.energy, 4 * result.energy, 1e-14, "Penalty energy changes with timestep by design");
    near(half_dt.normal_impulses[0], 2 * result.normal_impulses[0], 1e-14, "Penalty impulse inverse timestep scaling");
    auto doubled_mass = masses;
    for (auto& mass : doubled_mass) mass *= 2;
    const auto doubled = cloth::evaluate_ground_penalty(positions, doubled_mass, {3}, .1, options);
    near(doubled.energy, 2 * result.energy, 1e-14, "Penalty mass scaling");
    options.enabled = false;
    const auto disabled = cloth::evaluate_ground_penalty(positions, masses, {3}, .1, options);
    near(disabled.energy, 0, 0, "Disabled penalty is zero");
}

void invalid_inputs_and_overflow() {
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double largest = std::numeric_limits<double>::max();
    const auto zero_friction = cloth::apply_coulomb_friction({{largest,0,0}}, {1}, {1}, {}, {0,0});
    near(zero_friction.velocities[0][0],largest,0,"Zero friction preserves huge finite speed");
    near(zero_friction.info.friction_dissipation,0,0,"Zero impulse must not produce zero-times-infinity NaN");
    const cloth::SimulationState state{{{0, -.1, 0}}, {{1, -2, 0}}};
    for (const double invalid : {0.0, -1.0, infinity, nan}) {
        rejects<std::invalid_argument>([&] { (void)cloth::apply_ground_contact(state, {invalid}, {0}, .1, ground()); },
                                       "All masses, including fixed masses, must be valid");
    }
    for (const double invalid : {0.0, -1.0, infinity, nan, 1e-300, largest}) {
        rejects<std::invalid_argument>([&] { (void)cloth::apply_ground_contact(state, {1}, {}, invalid, ground()); },
                                       "Invalid timestep rejected");
    }
    rejects<std::invalid_argument>([&] { (void)cloth::apply_ground_contact(state, {}, {}, .1, ground()); }, "Wrong mass size");
    rejects<std::invalid_argument>([&] { (void)cloth::apply_ground_contact(state, {1}, {1}, .1, ground()); }, "Invalid fixed index");
    rejects<std::invalid_argument>([&] { (void)cloth::apply_ground_contact(state, {1}, {0, 0}, .1, ground()); }, "Duplicate fixed index");
    rejects<std::invalid_argument>([&] { (void)cloth::apply_ground_contact({state.positions, {}}, {1}, {}, .1, ground()); }, "Wrong velocity size");
    for (const double invalid : {-1.0, infinity, nan}) {
        rejects<std::invalid_argument>([&] { (void)cloth::apply_coulomb_friction(state.velocities, {invalid}, {1}, {}, {}); }, "Invalid normal impulse");
    }
    rejects<std::invalid_argument>([&] { (void)cloth::apply_coulomb_friction(state.velocities, {}, {1}, {}, {}); }, "Wrong impulse size");
    rejects<std::invalid_argument>([&] { (void)cloth::apply_coulomb_friction(state.velocities, {1}, {1}, {}, {.4, .3}); }, "Static friction smaller than kinetic");
    rejects<std::invalid_argument>([&] { (void)cloth::apply_coulomb_friction(state.velocities, {1}, {1}, {}, {-1, std::nullopt}); }, "Negative friction");
    auto bad_state = state;
    bad_state.positions[0][0] = nan;
    rejects<std::invalid_argument>([&] { (void)cloth::apply_ground_contact(bad_state, {1}, {}, .1, {}); }, "Disabled operation still validates state");
    for (int field = 0; field < 7; ++field) {
        auto options = ground();
        if (field == 0) options.ground_y = nan;
        if (field == 1) options.restitution = 1.1;
        if (field == 2) options.stabilization_beta = -.1;
        if (field == 3) options.proximity_tolerance = -.1;
        if (field == 4) options.penalty_eta = 0;
        if (field == 5) options.friction.kinetic = infinity;
        if (field == 6) options.friction.static_coefficient = nan;
        rejects<std::invalid_argument>([&] { cloth::validate_ground_contact_options(options); }, "Invalid ground option");
    }
    rejects<std::overflow_error>([&] { (void)cloth::apply_ground_contact(state, {largest}, {}, .1, ground()); }, "Normal impulse overflow");
    check(state.positions[0][1] == -.1 && state.velocities[0][1] == -2, "Overflow cannot modify contact input");
    auto extreme_ground = ground();
    extreme_ground.ground_y = largest;
    const cloth::SimulationState extreme{{{0, -largest, 0}}, {{0, 0, 0}}};
    rejects<std::overflow_error>([&] { (void)cloth::apply_ground_contact(extreme, {1}, {}, .1, extreme_ground); }, "Penetration overflow");
    rejects<std::overflow_error>([] { (void)cloth::apply_coulomb_friction({{1e200, 0, 0}}, {1e200}, {1}, {}, {1, 1}); }, "Friction energy overflow");
    auto penalty_options = ground();
    penalty_options.penalty_eta = largest;
    rejects<std::overflow_error>([&] { (void)cloth::evaluate_ground_penalty(state.positions, {2}, {}, .1, penalty_options); }, "Penalty coefficient overflow");
}

} // namespace

int main() {
    try {
        friction_hand_calculation();
        friction_physical_properties();
        ground_impulses_and_energy();
        fixed_disabled_and_empty();
        penalty_derivatives_and_scaling();
        invalid_inputs_and_overflow();
        std::cout << "ground contact, Coulomb friction and penalty tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "contact test failed: " << error.what() << '\n';
        return 1;
    }
}
