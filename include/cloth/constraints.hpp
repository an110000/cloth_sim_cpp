#pragma once

#include "cloth/material.hpp"

#include <array>
#include <cstddef>
#include <variant>
#include <vector>

namespace cloth {

// Keep the first three indices stable for the task-04 API and reports.
enum class ConstraintGroup { structural, shear, bending, neohookean_energy };

struct DistanceConstraint {
    std::size_t i = 0;
    std::size_t j = 0;
    double rest_length = 0.0;
    double compliance = 0.0;
    double lambda_value = 0.0;
    ConstraintGroup group = ConstraintGroup::structural;
};

struct DihedralConstraint {
    // Oriented adjacent faces are (i,j,k) and (j,i,l).
    std::array<std::size_t, 4> nodes{};
    double rest_angle = 0.0;
    double compliance = 1.0e-5;
    double lambda_value = 0.0;
};

struct NeoHookeanEnergyConstraint {
    std::array<std::size_t, 3> nodes{};
    std::array<Vec2, 3> reference{};
    Material material{};
    double compliance = 1.0e-8;
    double energy_scale = 1.0; // Joules; C=element_energy/energy_scale is dimensionless.
    double lambda_value = 0.0;
};

using Constraint = std::variant<DistanceConstraint, DihedralConstraint,
                                 NeoHookeanEnergyConstraint>;

struct ConstraintEvaluation {
    double value = 0.0;
    double angle = 0.0; // Signed angle for a dihedral; unused for distance.
    std::array<std::size_t, 4> nodes{};
    std::array<Vec3, 4> gradients{};
    std::size_t count = 0;
    bool valid_geometry = true;
};

struct ProjectionResult {
    double delta_lambda = 0.0;
    double value_before = 0.0;
    bool valid_geometry = true;
    bool applied = false;
};

double wrap_angle(double angle);
ConstraintGroup constraint_group(const Constraint& constraint) noexcept;
const char* group_name(ConstraintGroup group) noexcept;
double constraint_compliance(const Constraint& constraint) noexcept;
double constraint_multiplier(const Constraint& constraint) noexcept;
void reset_constraint(Constraint& constraint) noexcept;
void validate_constraint(const Constraint& constraint, std::size_t vertex_count);
ConstraintEvaluation evaluate_constraint(const Constraint& constraint,
                                          const std::vector<Vec3>& positions);

// In-place Gauss-Seidel projection. Fixed nodes have inverse_mass=0. Invalid
// geometry returns valid_geometry=false and value_before=+infinity without
// changing positions/multiplier. Numeric errors throw before any local write.
ProjectionResult project_constraint(Constraint& constraint,
                                     std::vector<Vec3>& positions,
                                     const std::vector<double>& inverse_masses,
                                     double dt);

// Interior hinges retain the Python dictionary's first-seen edge order.
// Boundary edges are omitted; inconsistent winding/nonmanifold edges rejected.
std::vector<Constraint> build_bending_constraints(const Mesh& mesh,
                                                 double compliance = 1.0e-5);
// One scalar energy constraint per triangle, in the original triangle order.
// Uses material.thickness as the single elastic thickness source.
std::vector<Constraint> build_neohookean_constraints(
    const Mesh& mesh, const Material& material,
    double compliance = 1.0e-8, double energy_scale = 1.0);
std::vector<Constraint> build_rectangular_constraints(
    const Mesh& mesh, std::size_t nx, std::size_t ny,
    double structural_compliance = 1.0e-8, double shear_compliance = 1.0e-8,
    double bending_compliance = 1.0e-5, bool include_bending = true);

} // namespace cloth
