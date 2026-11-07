#include "cloth/material.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
void close(double actual, double expected, double tolerance, const std::string& message) {
    check(std::isfinite(actual) && std::abs(actual - expected) <= tolerance, message);
}
template<class F> void must_throw(F function) {
    bool threw = false;
    try { function(); } catch (const std::exception&) { threw = true; }
    check(threw, "Expected invalid input to be rejected");
}
const std::array<cloth::Vec2, 3> reference{{{0, 0}, {1, 0}, {0, 1}}};

void finite_difference(std::array<cloth::Vec3, 3> x) {
    const cloth::Material material{};
    const auto analytic = cloth::evaluate_triangle(reference, x, material);
    double error2 = 0.0;
    double norm2 = 0.0;
    const double h = 1.0e-7;
    for (std::size_t vertex = 0; vertex < 3; ++vertex) {
        for (std::size_t axis = 0; axis < 3; ++axis) {
            auto plus = x;
            auto minus = x;
            plus[vertex][axis] += h;
            minus[vertex][axis] -= h;
            const double numeric =
                (cloth::evaluate_triangle(reference, plus, material).energy -
                 cloth::evaluate_triangle(reference, minus, material).energy) / (2.0 * h);
            const double delta = numeric - analytic.gradient[vertex][axis];
            error2 += delta * delta;
            norm2 += numeric * numeric;
        }
    }
    check(std::sqrt(error2) / std::max(1.0, std::sqrt(norm2)) < 2.0e-7,
          "Analytic gradient disagrees with energy finite difference");
}
} // namespace

int main() {
    try {
        const std::array<cloth::Vec3, 3> rest{{{0, 0, 0}, {1, 0, 0}, {0, 0, 1}}};
        const auto zero = cloth::evaluate_triangle(reference, rest, {});
        close(zero.energy, 0.0, 1.0e-12, "Rest energy");
        close(zero.area_ratio, 1.0, 1.0e-14, "Rest area ratio");
        for (const auto& vertex : zero.gradient)
            for (const auto value : vertex) close(value, 0.0, 1.0e-12, "Rest gradient");

        const std::array<cloth::Vec3, 3> deformed{{
            {0.0, 0.0, 0.0}, {1.1, 0.1, 0.2}, {-0.1, 1.0, -0.05}}};
        finite_difference(deformed);
        // Exercise the actual J < j_min extension, not only ordinary strains.
        finite_difference({{{0, 0, 0}, {0.015, 0, 0}, {0, 0, 1}}});
        const auto original = cloth::evaluate_triangle(reference, deformed, {});
        auto moved = deformed;
        for (auto& vertex : moved) {
            const auto old = vertex;
            vertex = {-old[1] + 2.0, old[0] - 3.0, old[2] + 4.0};
        }
        const auto rotated = cloth::evaluate_triangle(reference, moved, {});
        close(rotated.energy, original.energy, 1.0e-9, "Rigid-motion energy invariance");
        for (std::size_t i = 0; i < 3; ++i) {
            close(rotated.gradient[i][0], -original.gradient[i][1], 1.0e-9, "Rotated gradient x");
            close(rotated.gradient[i][1], original.gradient[i][0], 1.0e-9, "Rotated gradient y");
            close(rotated.gradient[i][2], original.gradient[i][2], 1.0e-9, "Rotated gradient z");
        }
        for (std::size_t axis = 0; axis < 3; ++axis)
            close(original.gradient[0][axis] + original.gradient[1][axis] +
                  original.gradient[2][axis], 0.0, 1.0e-10, "Internal forces must balance");
        auto thicker = cloth::Material{};
        thicker.thickness *= 2.0;
        const auto twice = cloth::evaluate_triangle(reference, deformed, thicker);
        close(twice.energy, 2.0 * original.energy, 1.0e-10, "Thickness energy scaling");
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t axis = 0; axis < 3; ++axis)
                close(twice.gradient[i][axis], 2.0 * original.gradient[i][axis],
                      1.0e-10, "Thickness gradient scaling");

        auto invalid = cloth::Material{};
        invalid.poisson = 0.5;
        must_throw([&] { cloth::evaluate_triangle(reference, rest, invalid); });
        invalid = {};
        invalid.thickness = std::numeric_limits<double>::quiet_NaN();
        must_throw([&] { cloth::evaluate_triangle(reference, rest, invalid); });
        auto collapsed = rest;
        collapsed[2] = collapsed[1];
        must_throw([&] { cloth::evaluate_triangle(reference, collapsed, {}); });
        auto bad_reference = reference;
        bad_reference[2] = bad_reference[1];
        must_throw([&] { cloth::evaluate_triangle(bad_reference, rest, {}); });
        auto bad_position = rest;
        bad_position[0][0] = std::numeric_limits<double>::infinity();
        must_throw([&] { cloth::evaluate_triangle(reference, bad_position, {}); });
        std::cout << "Material: rest, derivatives, rigid motion, thickness, invalid inputs passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
