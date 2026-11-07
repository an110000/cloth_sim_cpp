#include "cloth/material.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
using Reference = std::array<cloth::Vec2, 3>;
using Positions = std::array<cloth::Vec3, 3>;

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

double norm(const cloth::Mat9& matrix) {
    double total = 0.0;
    for (const auto& row : matrix)
        for (const auto value : row) total += value * value;
    return std::sqrt(total);
}

void relative_matrix(const cloth::Mat9& actual, const cloth::Mat9& expected,
                     double tolerance, const std::string& message) {
    cloth::Mat9 difference{};
    for (std::size_t i = 0; i < 9; ++i)
        for (std::size_t j = 0; j < 9; ++j) {
            check(std::isfinite(actual[i][j]), message + ": non-finite value");
            difference[i][j] = actual[i][j] - expected[i][j];
        }
    const double relative = norm(difference) / std::max(1.0, norm(expected));
    check(relative < tolerance, message + ": relative error " + std::to_string(relative));
}

template<class Exception, class Function>
void must_throw(Function function, const std::string& message) {
    try { function(); } catch (const Exception&) { return; }
    throw std::runtime_error(message);
}

void check_derivative_case(const Reference& reference, const Positions& positions,
                           const cloth::Material& material, double step,
                           double tolerance, const std::string& name) {
    const auto result = cloth::evaluate_triangle_with_hessian(reference, positions, material);
    const auto first = cloth::evaluate_triangle(reference, positions, material);
    check(result.value.energy == first.energy && result.value.area_ratio == first.area_ratio,
          name + ": first-order values changed");
    check(result.value.gradient == first.gradient, name + ": first-order gradient changed");
    cloth::Mat9 numeric{};
    cloth::Mat9 transpose{};
    for (std::size_t column = 0; column < 9; ++column) {
        auto plus = positions;
        auto minus = positions;
        plus[column / 3][column % 3] += step;
        minus[column / 3][column % 3] -= step;
        const auto gp = cloth::evaluate_triangle(reference, plus, material).gradient;
        const auto gm = cloth::evaluate_triangle(reference, minus, material).gradient;
        for (std::size_t row = 0; row < 9; ++row) {
            numeric[row][column] = (gp[row / 3][row % 3] - gm[row / 3][row % 3]) / (2.0 * step);
            transpose[row][column] = result.hessian[column][row];
        }
    }
    relative_matrix(result.hessian, numeric, tolerance, name + ": gradient difference");
    relative_matrix(result.hessian, transpose, 2.0e-12, name + ": raw symmetry");

    // Uniform translations leave each gradient component unchanged.
    for (std::size_t axis = 0; axis < 3; ++axis) {
        double null_norm2 = 0.0;
        for (std::size_t row = 0; row < 9; ++row) {
            const double entry = result.hessian[row][axis] +
                                 result.hessian[row][3 + axis] + result.hessian[row][6 + axis];
            null_norm2 += entry * entry;
        }
        check(std::sqrt(null_norm2) / std::max(1.0, norm(result.hessian)) < 1.0e-12,
              name + ": translation null direction");
    }
}

void test_finite_differences() {
    const Reference unit{{{0, 0}, {1, 0}, {0, 1}}};
    const Reference arbitrary{{{0.2, -0.1}, {1.3, 0.2}, {-0.2, 0.8}}};
    const Positions deformed{{{0, 0, 0}, {1.1, 0.1, 0.2}, {-0.1, 1.0, -0.05}}};
    check_derivative_case(unit, {{{0, 0, 0}, {1, 0, 0}, {0, 0, 1}}}, {}, 1.0e-6, 2.0e-8, "rest");
    check_derivative_case(unit, deformed, {}, 1.0e-6, 2.0e-8, "ordinary strain");
    check_derivative_case(arbitrary, deformed, {}, 1.0e-6, 2.0e-8, "arbitrary reference");
    auto reversed_reference = arbitrary;
    auto reversed_positions = deformed;
    std::swap(reversed_reference[1], reversed_reference[2]);
    std::swap(reversed_positions[1], reversed_positions[2]);
    check_derivative_case(reversed_reference, reversed_positions, {}, 1.0e-6, 2.0e-8, "reversed winding");
    check_derivative_case(unit, {{{0, 0, 0}, {0.015, 0, 0}, {0, 0, 1}}}, {},
                          1.0e-8, 2.0e-7, "regularized strain");
    check_derivative_case(unit, {{{0, 0, 0}, {0.049, 0, 0}, {0, 0, 1}}}, {},
                          1.0e-8, 2.0e-7, "below join");
    check_derivative_case(unit, {{{0, 0, 0}, {0.051, 0, 0}, {0, 0, 1}}}, {},
                          1.0e-8, 2.0e-7, "above join");
    // The extension is C2, but not C3. Across the join, central differentiation
    // of the gradient can have O(step) error rather than ordinary O(step^2).
    check_derivative_case(unit, {{{0, 0, 0}, {0.05, 0, 0}, {0, 0, 1}}}, {},
                          1.0e-8, 3.0e-6, "at C2 join");
    auto pure_log = cloth::Material{};
    pure_log.j_min = 0.0;
    check_derivative_case(unit, deformed, pure_log, 1.0e-6, 2.0e-8, "unregularized logarithm");
}

void test_covariance_and_scaling() {
    const Reference reference{{{0.2, -0.1}, {1.3, 0.2}, {-0.2, 0.8}}};
    const Positions positions{{{0, 0, 0}, {1.1, 0.1, 0.2}, {-0.1, 1.0, -0.05}}};
    const auto original = cloth::evaluate_triangle_with_hessian(reference, positions, {});
    auto transformed = positions;
    for (auto& vertex : transformed) {
        const auto old = vertex;
        vertex = {-old[1] + 2.0, old[0] - 3.0, old[2] + 4.0};
    }
    const auto rotated = cloth::evaluate_triangle_with_hessian(reference, transformed, {});
    // Q = [[0,-1,0],[1,0,0],[0,0,1]]. Compare H' = block(Q) H block(Q)^T.
    const std::array<std::size_t, 3> source_axis{{1, 0, 2}};
    const std::array<double, 3> sign{{-1, 1, 1}};
    cloth::Mat9 expected{};
    for (std::size_t i = 0; i < 9; ++i)
        for (std::size_t j = 0; j < 9; ++j)
            expected[i][j] = sign[i % 3] * sign[j % 3] *
                original.hessian[3 * (i / 3) + source_axis[i % 3]]
                                [3 * (j / 3) + source_axis[j % 3]];
    relative_matrix(rotated.hessian, expected, 3.0e-12, "rigid-motion covariance");

    auto thicker = cloth::Material{};
    thicker.thickness *= 2.0;
    const auto doubled = cloth::evaluate_triangle_with_hessian(reference, positions, thicker);
    for (std::size_t i = 0; i < 9; ++i)
        for (std::size_t j = 0; j < 9; ++j) expected[i][j] = 2.0 * original.hessian[i][j];
    relative_matrix(doubled.hessian, expected, 1.0e-14, "thickness scaling");

    // Swapping local vertices must permute both matrix indices, not change physics.
    auto reversed_reference = reference;
    auto reversed_positions = positions;
    std::swap(reversed_reference[1], reversed_reference[2]);
    std::swap(reversed_positions[1], reversed_positions[2]);
    const auto reversed = cloth::evaluate_triangle_with_hessian(reversed_reference, reversed_positions, {});
    const std::array<std::size_t, 3> permutation{{0, 2, 1}};
    for (std::size_t i = 0; i < 9; ++i)
        for (std::size_t j = 0; j < 9; ++j)
            expected[i][j] = original.hessian[3 * permutation[i / 3] + i % 3]
                                           [3 * permutation[j / 3] + j % 3];
    relative_matrix(reversed.hessian, expected, 3.0e-12, "winding permutation");
}

void test_invalid_inputs() {
    const Reference reference{{{0, 0}, {1, 0}, {0, 1}}};
    const Positions rest{{{0, 0, 0}, {1, 0, 0}, {0, 0, 1}}};
    auto collapsed = rest;
    collapsed[2] = collapsed[1];
    must_throw<std::domain_error>([&] {
        cloth::evaluate_triangle_with_hessian(reference, collapsed, {});
    }, "Collapsed current triangle accepted");
    auto bad_reference = reference;
    bad_reference[2] = bad_reference[1];
    must_throw<std::domain_error>([&] {
        cloth::evaluate_triangle_with_hessian(bad_reference, rest, {});
    }, "Degenerate reference accepted");
    auto bad = rest;
    bad[0][0] = std::numeric_limits<double>::quiet_NaN();
    must_throw<std::invalid_argument>([&] {
        cloth::evaluate_triangle_with_hessian(reference, bad, {});
    }, "NaN position accepted");
    auto material = cloth::Material{};
    material.j_min = -1.0;
    must_throw<std::invalid_argument>([&] {
        cloth::evaluate_triangle_with_hessian(reference, rest, material);
    }, "Invalid continuation accepted");
}
} // namespace

int main() {
    try {
        test_finite_differences();
        test_covariance_and_scaling();
        test_invalid_inputs();
        std::cout << "Hessian: analytic derivatives, raw symmetry, translations, rotations, scaling, invalid inputs passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
