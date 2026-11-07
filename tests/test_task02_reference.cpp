#include "cloth/assembly.hpp"
#include "cloth/dofs.hpp"
#include "cloth/mass.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {
class Reader {
public:
    explicit Reader(const std::string& path) {
        std::ifstream file(path);
        if (!file) throw std::runtime_error("Cannot open task-02 reference data");
        std::string text, line;
        while (std::getline(file, line)) {
            if (!line.empty() && line.front() != '#') text += line + '\n';
        }
        input_.str(text);
    }
    template<class T> T get() {
        T value{};
        if (!(input_ >> value)) throw std::runtime_error("Truncated/malformed reference data");
        return value;
    }
    void expect(const std::string& expected) {
        const auto found = get<std::string>();
        if (found != expected) throw std::runtime_error("Expected " + expected + ", got " + found);
    }
    void finish() {
        std::string extra;
        if (input_ >> extra) throw std::runtime_error("Unexpected trailing reference data");
    }
private:
    std::istringstream input_;
};

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
void compare(double actual, double expected, double absolute, const std::string& name) {
    const double limit = absolute + 2.0e-10 * std::abs(expected);
    if (!std::isfinite(actual) || !std::isfinite(expected) || std::abs(actual - expected) > limit) {
        std::ostringstream detail;
        detail.precision(17);
        detail << name << ": actual=" << actual << " expected=" << expected
               << " error=" << std::abs(actual - expected) << " limit=" << limit;
        throw std::runtime_error(detail.str());
    }
}
cloth::Material read_material(Reader& reader) {
    reader.expect("material");
    return {reader.get<double>(), reader.get<double>(), reader.get<double>(), reader.get<double>()};
}
void compare_values(Reader& reader, const std::string& tag, const std::vector<double>& actual,
                    double absolute) {
    reader.expect(tag);
    for (std::size_t i = 0; i < actual.size(); ++i)
        compare(actual[i], reader.get<double>(), absolute, tag + "[" + std::to_string(i) + "]");
}
} // namespace

int main(int argc, char** argv) {
    try {
        check(argc == 2, "Expected task-02 fixture path");
        Reader reader(argv[1]);
        reader.expect("CLOTH_TASK02_REFERENCE");
        check(reader.get<int>() == 1, "Unsupported reference format version");
        reader.expect("elements");
        const auto element_count = reader.get<std::size_t>();
        check(element_count == 11, "Expected all 11 element reference cases");
        std::set<std::string> remaining{"rest", "deformed", "thicker", "pure_log", "rigid_motion",
            "regularized", "below_join", "at_join", "above_join", "arbitrary_reference", "reversed_winding"};
        for (std::size_t e = 0; e < element_count; ++e) {
            reader.expect("element");
            const auto name = reader.get<std::string>();
            check(remaining.erase(name) == 1, "Duplicate/unknown element case: " + name);
            const auto material = read_material(reader);
            std::array<cloth::Vec2, 3> X{};
            std::array<cloth::Vec3, 3> x{};
            reader.expect("reference");
            for (auto& point : X) for (auto& v : point) v = reader.get<double>();
            reader.expect("positions");
            for (auto& point : x) for (auto& v : point) v = reader.get<double>();
            const auto result = cloth::evaluate_triangle_with_hessian(X, x, material);
            reader.expect("energy");
            compare(result.value.energy, reader.get<double>(), 1.0e-9, name + " energy");
            reader.expect("gradient");
            for (const auto& point : result.value.gradient)
                for (const auto v : point) compare(v, reader.get<double>(), 1.0e-9, name + " gradient");
            reader.expect("area_ratio");
            compare(result.value.area_ratio, reader.get<double>(), 1.0e-12, name + " J");
            reader.expect("hessian");
            for (const auto& row : result.hessian)
                for (const auto v : row) compare(v, reader.get<double>(), 2.0e-8, name + " Hessian");
            reader.expect("end");
        }
        reader.expect("systems");
        const auto system_count = reader.get<std::size_t>();
        check(system_count == 3, "Expected all three global reference cases");
        remaining = {"two_triangles", "nonuniform", "regularized_strip"};
        for (std::size_t s = 0; s < system_count; ++s) {
            reader.expect("system");
            const auto name = reader.get<std::string>();
            check(remaining.erase(name) == 1, "Duplicate/unknown system case: " + name);
            reader.expect("vertices");
            const auto n = reader.get<std::size_t>();
            check(n > 0 && n <= 100, "Invalid reference mesh size");
            cloth::Mesh mesh;
            mesh.reference_positions.resize(n);
            mesh.positions.resize(n);
            reader.expect("reference");
            for (auto& p : mesh.reference_positions) for (auto& v : p) v = reader.get<double>();
            reader.expect("positions");
            for (auto& p : mesh.positions) for (auto& v : p) v = reader.get<double>();
            reader.expect("triangles");
            const auto triangles = reader.get<std::size_t>();
            check(triangles > 0 && triangles <= 200, "Invalid reference triangle count");
            mesh.triangles.resize(triangles);
            for (auto& tri : mesh.triangles) for (auto& i : tri) i = reader.get<std::size_t>();
            const auto material = read_material(reader);
            reader.expect("density");
            const auto density = reader.get<double>();
            reader.expect("fixed");
            const auto fixed_count = reader.get<std::size_t>();
            check(fixed_count <= n, "Invalid fixed vertex count");
            std::vector<std::size_t> fixed(fixed_count);
            for (auto& i : fixed) i = reader.get<std::size_t>();
            const auto result = cloth::assemble_elastic(mesh, material);
            reader.expect("energy");
            compare(result.energy, reader.get<double>(), 1.0e-9, name + " energy");
            compare_values(reader, "gradient", result.gradient, 1.0e-9);
            compare_values(reader, "hessian", result.hessian.dense_row_major(), 2.0e-8);
            const auto lumped = cloth::assemble_mass_matrix(mesh, density, material.thickness, cloth::MassMode::lumped);
            const auto consistent = cloth::assemble_mass_matrix(mesh, density, material.thickness, cloth::MassMode::consistent);
            compare_values(reader, "mass_lumped", lumped.dense_row_major(), 1.0e-13);
            compare_values(reader, "mass_consistent", consistent.dense_row_major(), 1.0e-13);
            const cloth::DofMap map(n, fixed);
            reader.expect("free_indices");
            check(reader.get<std::size_t>() == map.free_size(), "Free count differs from Python");
            for (const auto i : map.free_to_full())
                check(i == reader.get<std::size_t>(), "Free ordering differs from Python");
            compare_values(reader, "hessian_free", cloth::restrict_to_free(result.hessian, map).dense_row_major(), 2.0e-8);
            reader.expect("end");
        }
        reader.expect("done");
        reader.finish();
        std::cout << "11 analytic element Hessians and 3 assembled Python systems passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
