#include "cloth/constraints.hpp"
#include "cloth/experiment.hpp"
#include "cloth/scene.hpp"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct ParsedScene {
    cloth::SceneConfigInput input;
    std::map<std::size_t, cloth::Vec3> position_overrides;
    std::map<std::size_t, cloth::Vec3> velocity_overrides;
    std::map<std::size_t, cloth::Vec2> mesh_references;
    std::map<std::size_t, cloth::Vec3> mesh_positions;
    std::map<std::size_t, cloth::Triangle> mesh_triangles;
    bool geometry_changed = false;
};

std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

double number(const std::string& text, const char* name) {
    std::size_t used = 0;
    double value = 0.0;
    try { value = std::stod(text, &used); }
    catch (...) { throw std::invalid_argument(std::string(name) + " must be a finite number"); }
    if (used != text.size() || !std::isfinite(value))
        throw std::invalid_argument(std::string(name) + " must be a finite number");
    return value;
}

std::size_t integer(const std::string& text, const char* name) {
    std::size_t value = 0;
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    if (text.empty() || parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size())
        throw std::invalid_argument(std::string(name) + " must be a nonnegative integer");
    return value;
}

bool boolean(const std::string& text, const char* name) {
    if (text == "true" || text == "1") return true;
    if (text == "false" || text == "0") return false;
    throw std::invalid_argument(std::string(name) + " must be true or false");
}

std::vector<std::string> split(const std::string& text) {
    std::vector<std::string> values;
    std::string value;
    std::istringstream stream(text);
    while (std::getline(stream, value, ',')) values.push_back(trim(value));
    if (values.size() == 1 && values.front().empty()) values.clear();
    return values;
}

cloth::Vec3 vector3(const std::string& text, const char* name) {
    const auto values = split(text);
    if (values.size() != 3) throw std::invalid_argument(std::string(name) + " needs three comma-separated numbers");
    return {number(values[0], name), number(values[1], name), number(values[2], name)};
}

cloth::Vec2 vector2(const std::string& text, const char* name) {
    const auto values = split(text);
    if (values.size() != 2) throw std::invalid_argument(std::string(name) + " needs two comma-separated numbers");
    return {number(values[0], name), number(values[1], name)};
}

cloth::Triangle triangle(const std::string& text, const char* name) {
    const auto values = split(text);
    if (values.size() != 3) throw std::invalid_argument(std::string(name) + " needs three comma-separated indices");
    return {integer(values[0], name), integer(values[1], name), integer(values[2], name)};
}

std::vector<std::size_t> indices(const std::string& text, const char* name) {
    std::vector<std::size_t> values;
    for (const auto& token : split(text)) values.push_back(integer(token, name));
    return values;
}

void assign(ParsedScene& parsed, const std::string& key, const std::string& value) {
    auto& input = parsed.input;
    if (key == "base") {
        if (value != "teaching" && value != "generated") throw std::invalid_argument("base must be teaching or generated");
        if (value == "teaching") input = cloth::default_scene_config_input();
    } else if (key == "solver") {
        if (value == "implicit") input.solver = cloth::SceneSolver::implicit_euler;
        else if (value == "xpbd-distance") input.solver = cloth::SceneSolver::xpbd_distance;
        else if (value == "xpbd-energy") input.solver = cloth::SceneSolver::xpbd_energy;
        else throw std::invalid_argument("solver must be implicit, xpbd-distance or xpbd-energy");
    } else if (key == "mass_mode") {
        if (value == "lumped") input.mass_mode = cloth::MassMode::lumped;
        else if (value == "consistent") input.mass_mode = cloth::MassMode::consistent;
        else throw std::invalid_argument("mass_mode must be lumped or consistent");
    } else if (key == "nx") { input.rectangle.nx = integer(value, "nx"); parsed.geometry_changed = true; }
    else if (key == "ny") { input.rectangle.ny = integer(value, "ny"); parsed.geometry_changed = true; }
    else if (key == "width") { input.rectangle.width = number(value, "width"); parsed.geometry_changed = true; }
    else if (key == "height") { input.rectangle.height = number(value, "height"); parsed.geometry_changed = true; }
    else if (key == "density") input.density = number(value, "density");
    else if (key == "young") input.material.young = number(value, "young");
    else if (key == "poisson") input.material.poisson = number(value, "poisson");
    else if (key == "thickness") input.material.thickness = number(value, "thickness");
    else if (key == "j_min") input.material.j_min = number(value, "j_min");
    else if (key == "dt") input.dt = number(value, "dt");
    else if (key == "steps") input.steps = integer(value, "steps");
    else if (key == "run_seed") input.run_seed = static_cast<std::uint64_t>(integer(value, "run_seed"));
    else if (key == "gravity") input.gravity = vector3(value, "gravity");
    else if (key == "fixed_vertices") input.fixed_vertices = indices(value, "fixed_vertices");
    else if (key == "ground_enabled") input.contact.enabled = boolean(value, "ground_enabled");
    else if (key == "ground_y") input.contact.ground_y = number(value, "ground_y");
    else if (key == "friction") input.contact.friction.kinetic = number(value, "friction");
    else if (key == "static_friction") input.contact.friction.static_coefficient = number(value, "static_friction");
    else if (key == "restitution") input.contact.restitution = number(value, "restitution");
    else if (key == "contact_beta") input.contact.stabilization_beta = number(value, "contact_beta");
    else if (key == "contact_eta") input.contact.penalty_eta = number(value, "contact_eta");
    else if (key == "asset_sha256") input.asset_sha256 = value;
    else if (key == "scene_sha256") input.scene_sha256 = value;
    else if (key == "state_sha256") input.state_sha256 = value;
    else if (key == "newton_max_iterations") input.newton_options.max_iterations = integer(value, "newton_max_iterations");
    else if (key == "xpbd_iterations") input.xpbd_options.iterations = integer(value, "xpbd_iterations");
    else if (key.rfind("position.", 0) == 0) parsed.position_overrides[integer(key.substr(9), "position index")] = vector3(value, key.c_str());
    else if (key.rfind("velocity.", 0) == 0) parsed.velocity_overrides[integer(key.substr(9), "velocity index")] = vector3(value, key.c_str());
    else if (key.rfind("mesh_reference.", 0) == 0) parsed.mesh_references[integer(key.substr(15), "mesh reference index")] = vector2(value, key.c_str());
    else if (key.rfind("mesh_position.", 0) == 0) parsed.mesh_positions[integer(key.substr(14), "mesh position index")] = vector3(value, key.c_str());
    else if (key.rfind("mesh_triangle.", 0) == 0) parsed.mesh_triangles[integer(key.substr(14), "mesh triangle index")] = triangle(value, key.c_str());
    else throw std::invalid_argument("Unknown scene key: " + key);
}

template<class Map>
void require_contiguous(const Map& values, const char* name) {
    std::size_t expected = 0;
    for (const auto& [index, value] : values) {
        (void)value;
        if (index != expected++) throw std::invalid_argument(std::string(name) + " indices must be contiguous from zero");
    }
}

std::vector<cloth::Constraint> build_custom_distance_constraints(const cloth::Mesh& mesh) {
    std::set<std::pair<std::size_t, std::size_t>> edges;
    for (const auto& t : mesh.triangles) {
        for (const auto& edge : {std::pair{t[0], t[1]}, std::pair{t[1], t[2]}, std::pair{t[2], t[0]}}) {
            edges.insert(std::minmax(edge.first, edge.second));
        }
    }
    std::vector<cloth::Constraint> constraints;
    constraints.reserve(edges.size());
    for (const auto& [i, j] : edges) {
        const auto& a = mesh.positions[i]; const auto& b = mesh.positions[j];
        const double dx = a[0]-b[0], dy = a[1]-b[1], dz = a[2]-b[2];
        const double length = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (!std::isfinite(length) || length <= 0.0) throw std::invalid_argument("Custom mesh edge has invalid rest length");
        constraints.emplace_back(cloth::DistanceConstraint{i, j, length, 1.0e-8, 0.0,
                                                            cloth::ConstraintGroup::structural});
    }
    return constraints;
}

ParsedScene read_scene(const std::filesystem::path& path) {
    std::ifstream stream(path);
    if (!stream) throw std::invalid_argument("Cannot open scene file: " + path.string());
    ParsedScene parsed;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(stream, line)) {
        ++line_number;
        const auto comment = line.find('#');
        if (comment != std::string::npos) line.erase(comment);
        line = trim(line);
        if (line.empty()) continue;
        const auto separator = line.find('=');
        if (separator == std::string::npos) throw std::invalid_argument("Scene line lacks '=' at line " + std::to_string(line_number));
        try { assign(parsed, trim(line.substr(0, separator)), trim(line.substr(separator + 1))); }
        catch (const std::exception& error) { throw std::invalid_argument("Scene line " + std::to_string(line_number) + ": " + error.what()); }
    }
    if (!parsed.mesh_references.empty() || !parsed.mesh_positions.empty() || !parsed.mesh_triangles.empty()) {
        if (parsed.mesh_references.empty() || parsed.mesh_positions.empty() || parsed.mesh_triangles.empty() ||
            parsed.mesh_references.size() != parsed.mesh_positions.size())
            throw std::invalid_argument("Custom mesh needs matching reference, position and triangle sections");
        require_contiguous(parsed.mesh_references, "mesh_reference");
        require_contiguous(parsed.mesh_positions, "mesh_position");
        require_contiguous(parsed.mesh_triangles, "mesh_triangle");
        cloth::Mesh custom;
        custom.reference_positions.resize(parsed.mesh_references.size());
        custom.positions.resize(parsed.mesh_positions.size());
        custom.triangles.resize(parsed.mesh_triangles.size());
        for (const auto& [index, value] : parsed.mesh_references) custom.reference_positions[index] = value;
        for (const auto& [index, value] : parsed.mesh_positions) custom.positions[index] = value;
        for (const auto& [index, value] : parsed.mesh_triangles) custom.triangles[index] = value;
        parsed.input.custom_mesh = std::move(custom);
        parsed.input.initial_positions.clear();
        parsed.input.initial_velocities.clear();
    }
    const auto mesh = parsed.input.custom_mesh
        ? *parsed.input.custom_mesh
        : cloth::make_rectangular_mesh(parsed.input.rectangle.nx, parsed.input.rectangle.ny,
                                       parsed.input.rectangle.width, parsed.input.rectangle.height);
    const auto count = mesh.positions.size();
    if (parsed.geometry_changed) {
        parsed.input.initial_positions.clear();
        parsed.input.initial_velocities.clear();
    }
    if (!parsed.position_overrides.empty()) {
        parsed.input.initial_positions = mesh.positions;
        for (const auto& [index, value] : parsed.position_overrides) {
            if (index >= count) throw std::invalid_argument("position index is out of range");
            parsed.input.initial_positions[index] = value;
        }
    }
    if (!parsed.velocity_overrides.empty()) {
        parsed.input.initial_velocities.assign(count, cloth::Vec3{0.0, 0.0, 0.0});
        for (const auto& [index, value] : parsed.velocity_overrides) {
            if (index >= count) throw std::invalid_argument("velocity index is out of range");
            parsed.input.initial_velocities[index] = value;
        }
    }
    return parsed;
}

struct Arguments { std::filesystem::path scene; std::filesystem::path output; bool help = false; };

Arguments parse_arguments(int argc, char** argv) {
    Arguments arguments;
    for (int i = 1; i < argc; ++i) {
        const std::string option = argv[i];
        if (option == "--help") { arguments.help = true; continue; }
        if (option != "--scene" && option != "--output") throw std::invalid_argument("Unknown option: " + option);
        if (++i >= argc || std::string(argv[i]).empty()) throw std::invalid_argument("Missing value for " + option);
        if (option == "--scene") arguments.scene = argv[i]; else arguments.output = argv[i];
    }
    if (!arguments.help && (arguments.scene.empty() || arguments.output.empty()))
        throw std::invalid_argument("Both --scene and --output are required");
    return arguments;
}

std::string solver_name(cloth::SceneSolver solver) {
    if (solver == cloth::SceneSolver::implicit_euler) return "implicit_euler";
    if (solver == cloth::SceneSolver::xpbd_distance) return "xpbd_distance";
    return "xpbd_energy";
}

void fill_metadata(cloth::ExperimentRecord& record, const cloth::SceneConfig& scene) {
    auto& m = record.metadata;
    m.run_id = "scene-run-07b";
    m.source_sha256 = "<unset>";
    m.config_sha256 = "<unset>";
    m.scene_sha256 = scene.scene_sha256().empty() ? "<unset>" : scene.scene_sha256();
    m.model_identity = solver_name(scene.solver());
    m.mass_mode = scene.mass_mode() == cloth::MassMode::lumped ? "lumped" : "consistent";
    m.bending_model = scene.solver() == cloth::SceneSolver::xpbd_distance
                      ? (scene.has_custom_mesh() ? "none_custom_structural_edges" : "rectangular_distance_bending") :
                      (scene.solver() == cloth::SceneSolver::xpbd_energy ? "none" : "implicit_membrane");
    m.precision = "double"; m.backend = "cpu"; m.device = "host"; m.toolchain = "C++17";
    m.container_digest = "<unset>"; m.random_state = std::to_string(scene.run_seed());
    m.dt_seconds = scene.dt(); m.requested_steps = scene.steps();
    if (scene.solver() == cloth::SceneSolver::implicit_euler) {
        m.solver_budget = "newton_max_iterations=" + std::to_string(scene.newton_options().max_iterations);
        m.tolerances = "residual=" + std::to_string(scene.newton_options().residual_tolerance);
    } else {
        m.solver_budget = "xpbd_iterations=" + std::to_string(scene.xpbd_options().iterations);
        m.tolerances = "distance=" + std::to_string(scene.xpbd_options().distance_tolerance);
    }
}

int run(const Arguments& arguments) {
    const auto parsed = read_scene(arguments.scene);
    const cloth::SceneConfig scene(parsed.input);
    cloth::ExperimentRecord record; fill_metadata(record, scene);
    cloth::SimulationState state = scene.initial_state();
    const auto started = std::chrono::steady_clock::now();
    double solver_seconds = 0.0;
    auto fail = [&](const std::string& reason) {
        record.wall_clock_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        record.solver_seconds = solver_seconds;
        cloth::finish_experiment(record, cloth::RunStatus::failed_numerical, reason);
        cloth::write_experiment_text(arguments.output, record);
        return 2;
    };
    try {
        if (scene.solver() == cloth::SceneSolver::implicit_euler) {
            cloth::ImplicitEulerSolver solver(scene.mesh(), scene.material(), scene.density(), scene.mass_mode(),
                                              scene.gravity(), scene.newton_options(), scene.contact());
            for (std::size_t step = 1; step <= scene.steps(); ++step) {
                const auto before = std::chrono::steady_clock::now();
                auto result = solver.step(state, scene.fixed_vertices(), scene.dt());
                solver_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now()-before).count();
                if (!result.info.committed) return fail(result.info.message.empty() ? cloth::failure_name(result.info.failure_reason) : result.info.message);
                state = std::move(result.state);
                cloth::append_committed_frame(record, {step, static_cast<double>(step)*scene.dt(), state});
            }
        } else {
            const auto constraints = scene.solver() == cloth::SceneSolver::xpbd_distance
                ? (scene.has_custom_mesh() ? build_custom_distance_constraints(scene.mesh())
                                           : cloth::build_rectangular_constraints(scene.mesh(), scene.rectangle().nx, scene.rectangle().ny))
                : cloth::build_neohookean_constraints(scene.mesh(), scene.material());
            cloth::XpbdSolver solver(scene.masses(), constraints, scene.gravity(), scene.xpbd_options(), scene.contact());
            for (std::size_t step = 1; step <= scene.steps(); ++step) {
                const auto before = std::chrono::steady_clock::now();
                auto result = solver.step(state, scene.fixed_vertices(), scene.dt());
                solver_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now()-before).count();
                if (!result.info.committed) return fail(result.info.message.empty() ? cloth::failure_name(result.info.failure_reason) : result.info.message);
                state = std::move(result.state);
                cloth::append_committed_frame(record, {step, static_cast<double>(step)*scene.dt(), state});
            }
        }
        record.wall_clock_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        record.solver_seconds = solver_seconds;
        cloth::finish_experiment(record, cloth::RunStatus::committed);
        cloth::write_experiment_text(arguments.output, record);
        return 0;
    } catch (const std::exception& error) {
        if (record.status == cloth::RunStatus::not_run) {
            try { return fail(error.what()); } catch (...) {}
        }
        throw;
    }
}

} // namespace

int main(int argc, char** argv) {
    try {
        const auto arguments = parse_arguments(argc, argv);
        if (arguments.help) {
            std::cout << "cloth_scene_run --scene FILE --output REPORT\n"
                         "Scene files use key=value lines; comments begin with #.\n";
            return 0;
        }
        return run(arguments);
    } catch (const std::exception& error) {
        std::cerr << "cloth_scene_run: " << error.what() << '\n';
        return 1;
    }
}
