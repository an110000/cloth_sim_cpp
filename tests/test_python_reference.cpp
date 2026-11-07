#include "cloth/material.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
void compare(double actual, double expected, const std::string& name) {
    const double tolerance = 1.0e-9 + 2.0e-10 * std::abs(expected);
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(name + ": C++ differs from Python reference");
    }
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw std::runtime_error("Expected fixture path");
        std::ifstream input(argv[1]);
        if (!input) throw std::runtime_error("Cannot open Python reference fixture");
        std::string line;
        std::size_t count = 0;
        std::set<std::string> remaining{
            "rest", "deformed", "thicker", "pure_log", "rigid_motion",
            "regularized", "below_join", "at_join", "above_join",
            "arbitrary_reference", "reversed_winding"};
        while (std::getline(input, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream row(line);
            std::string name;
            cloth::Material material;
            std::array<cloth::Vec2, 3> reference{};
            std::array<cloth::Vec3, 3> positions{};
            double energy = 0.0;
            std::array<cloth::Vec3, 3> gradient{};
            double area_ratio = 0.0;
            row >> name >> material.young >> material.poisson >> material.thickness >> material.j_min;
            if (remaining.erase(name) != 1)
                throw std::runtime_error("Duplicate or unknown fixture: " + name);
            for (auto& vertex : reference) for (auto& value : vertex) row >> value;
            for (auto& vertex : positions) for (auto& value : vertex) row >> value;
            row >> energy;
            for (auto& vertex : gradient) for (auto& value : vertex) row >> value;
            row >> area_ratio;
            if (!row) throw std::runtime_error("Malformed fixture: " + name);
            std::string extra;
            if (row >> extra) throw std::runtime_error("Unexpected fixture fields: " + name);
            const auto actual = cloth::evaluate_triangle(reference, positions, material);
            compare(actual.energy, energy, name + " energy");
            compare(actual.area_ratio, area_ratio, name + " J");
            for (std::size_t i = 0; i < 3; ++i)
                for (std::size_t axis = 0; axis < 3; ++axis)
                    compare(actual.gradient[i][axis], gradient[i][axis], name + " gradient");
            ++count;
        }
        if (!remaining.empty()) throw std::runtime_error("Missing reference cases");
        std::cout << count << " Python reference cases passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
