#include "cloth/implicit.hpp"
#include "contact_output.hpp"

#include <charconv>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>

namespace {

struct Arguments {
    std::size_t steps = 8;
    double dt = 0.01;
    cloth::MassMode mass_mode = cloth::MassMode::lumped;
    std::filesystem::path output = "output/task03";
    std::size_t max_iterations = 50;
    bool help = false;
    cloth::GroundContactOptions ground;
};

void usage() {
    std::cout << "Usage: cloth_simulate [--steps N] [--dt SECONDS]\n"
                 "                      [--mass lumped|consistent] [--output DIRECTORY]\n"
                 "                      [--max-iterations N] [--help]\n"
                 "                      [--ground-y M] [--friction MU] [--static-friction MU_S]\n"
                 "                      [--contact-eta ETA]\n"
                 "Runs the fixed pinned-cloth teaching scene; ground contact is optional.\n"
                 "Defaults: 8 steps, dt=0.01, lumped mass, 50 Newton updates,\n"
                 "          output=output/task03. A failed step is rolled back.\n"
                 "Exit codes: 0 completed, 2 numerical failure, 1 argument or I/O error.\n";
}

std::size_t unsigned_integer(const std::string& text, bool allow_zero, const char* name) {
    std::size_t value = 0;
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    if (text.empty() || parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size() ||
        (!allow_zero && value == 0)) {
        throw std::invalid_argument(std::string(name) + " must be " +
                                    (allow_zero ? "a nonnegative integer" : "a positive integer"));
    }
    return value;
}

double positive_double(const std::string& text, const char* name) {
    if (text.empty() || text.find_first_not_of("0123456789.eE+-") != std::string::npos) {
        throw std::invalid_argument(std::string(name) + " must be finite and positive");
    }
    std::size_t consumed = 0;
    double value;
    try { value = std::stod(text, &consumed); }
    catch (const std::exception&) {
        throw std::invalid_argument(std::string(name) + " must be finite and positive");
    }
    if (consumed != text.size() || !std::isfinite(value) || value <= 0.0 ||
        !std::isfinite(value * value) || value * value == 0.0) {
        throw std::invalid_argument(std::string(name) + " or its square is outside the positive finite range");
    }
    return value;
}

Arguments parse(int argc, char** argv) {
    Arguments arguments;
    std::set<std::string> seen;
    for (int i = 1; i < argc; ++i) {
        const std::string option = argv[i];
        if (option == "--help") { arguments.help = true; continue; }
        if (option != "--steps" && option != "--dt" && option != "--mass" &&
            option != "--output" && option != "--max-iterations" && !cloth_cli::contact_option(option,false)) {
            throw std::invalid_argument("Unknown option: " + option);
        }
        if (!seen.insert(option).second) throw std::invalid_argument("Duplicate option: " + option);
        if (++i >= argc) throw std::invalid_argument("Missing value for " + option);
        const std::string value = argv[i];
        if (cloth_cli::contact_option(option,false)) cloth_cli::set_contact(arguments.ground,option,value);
        else if (option == "--steps") arguments.steps = unsigned_integer(value, false, "--steps");
        else if (option == "--dt") arguments.dt = positive_double(value, "--dt");
        else if (option == "--max-iterations") {
            arguments.max_iterations = unsigned_integer(value, true, "--max-iterations");
        } else if (option == "--mass") {
            if (value == "lumped") arguments.mass_mode = cloth::MassMode::lumped;
            else if (value == "consistent") arguments.mass_mode = cloth::MassMode::consistent;
            else throw std::invalid_argument("--mass must be lumped or consistent");
        } else {
            if (value.empty() || value.rfind("--", 0) == 0) {
                throw std::invalid_argument("--output requires a directory path");
            }
            arguments.output = std::filesystem::path(value);
        }
    }
    const long double duration = static_cast<long double>(arguments.steps) * arguments.dt;
    if (duration > std::numeric_limits<double>::max()) {
        throw std::invalid_argument("Requested simulation duration is outside the finite range");
    }
    cloth::validate_ground_contact_options(arguments.ground);
    for(const auto& option:seen)
        if(cloth_cli::contact_option(option,false) && !arguments.ground.enabled)
            throw std::invalid_argument("Contact parameters require --ground-y");
    return arguments;
}

std::ofstream output_file(const std::filesystem::path& path) {
    std::ofstream file(path);
    file.exceptions(std::ios::failbit | std::ios::badbit);
    file << std::setprecision(17) << std::boolalpha;
    return file;
}

void json_number(std::ostream& file, double value) {
    if (std::isfinite(value)) file << value;
    else file << "null";
}

void json_string(std::ostream& file, const std::string& value) {
    file << '"';
    for (const unsigned char ch : value) {
        if (ch == '"' || ch == '\\') file << '\\' << static_cast<char>(ch);
        else if (ch == '\n') file << "\\n";
        else if (ch == '\r') file << "\\r";
        else if (ch == '\t') file << "\\t";
        else if (ch < 0x20) {
            static constexpr char digits[] = "0123456789abcdef";
            file << "\\u00" << digits[ch >> 4] << digits[ch & 15];
        } else file << static_cast<char>(ch);
    }
    file << '"';
}

template<class Points> void json_points(std::ostream& file, const Points& points) {
    file << '[';
    for (std::size_t i = 0; i < points.size(); ++i) {
        if (i != 0) file << ',';
        file << '[';
        for (std::size_t axis = 0; axis < points[i].size(); ++axis) {
            if (axis != 0) file << ',';
            file << points[i][axis];
        }
        file << ']';
    }
    file << ']';
}

void write_config(const Arguments& args, const cloth::Mesh& mesh,
                  const cloth::Material& material, const cloth::SimulationState& initial,
                  const cloth::NewtonOptions& options) {
    auto file = output_file(args.output / "config.json");
    file << "{\n  \"format_version\": 1,\n  \"scene\": \"pinned_cloth\",\n"
            "  \"units\": {\"length\":\"m\",\"time\":\"s\",\"mass\":\"kg\",\"young\":\"Pa\"},\n"
            "  \"mesh\": {\"nx_cells\":2,\"ny_cells\":1,\"width\":1.0,\"height\":0.7,\n"
            "    \"reference_positions\":";
    json_points(file, mesh.reference_positions);
    file << ",\n    \"triangles\":";
    json_points(file, mesh.triangles);
    file << "},\n  \"material\": {\"young\":" << material.young
         << ",\"poisson\":" << material.poisson << ",\"thickness\":" << material.thickness
         << ",\"j_min\":" << material.j_min << "},\n"
            "  \"density\":800.0,\n  \"gravity\":[0.0,-9.81,0.0],\n"
            "  \"fixed_vertices\":[3,4,5],\n  \"mass_mode\":\""
         << (args.mass_mode == cloth::MassMode::lumped ? "lumped" : "consistent")
         << "\",\n  \"dt\":" << args.dt << ",\n  \"steps\":" << args.steps
         << ",\n  ";
    cloth_cli::config(file,args.ground);
    file << ",\n  \"initial_positions\":";
    json_points(file, initial.positions);
    file << ",\n  \"initial_velocities\":";
    json_points(file, initial.velocities);
    file << ",\n  \"solver\":{\"max_iterations\":" << options.max_iterations
         << ",\"residual_tolerance\":" << options.residual_tolerance
         << ",\"linear_relative_tolerance\":" << options.linear_relative_tolerance
         << ",\"armijo\":" << options.armijo
         << ",\"backtrack_factor\":" << options.backtrack_factor
         << ",\"line_search_max_iterations\":" << options.line_search_max_iterations
         << ",\"diagonal_shifts\":[";
    for (std::size_t i = 0; i < options.diagonal_shifts.size(); ++i) {
        if (i != 0) file << ',';
        file << options.diagonal_shifts[i];
    }
    file << "]}\n}\n";
    file.close();
}

void write_frame(std::ostream& csv, std::size_t step, double time,
                 const cloth::SimulationState& state) {
    for (std::size_t vertex = 0; vertex < state.positions.size(); ++vertex) {
        csv << step << ',' << time << ',' << vertex;
        for (const auto coordinate : state.positions[vertex]) csv << ',' << coordinate;
        for (const auto velocity : state.velocities[vertex]) csv << ',' << velocity;
        csv << '\n';
    }
}

void write_summary(const Arguments& args, std::size_t completed,
                   const cloth::StepInfo& info) {
    auto file = output_file(args.output / "summary.json");
    file << "{\n  \"format_version\":1,\n  \"steps_requested\":" << args.steps
         << ",\n  \"steps_completed\":" << completed << ",\n  \"time_completed\":";
    json_number(file, static_cast<double>(completed) * args.dt);
    file << ",\n  \"finite_state\":" << info.finite_state
         << ",\n  \"solver_converged\":" << info.solver_converged
         << ",\n  \"committed\":" << info.committed
         << ",\n  \"failure_reason\":\"" << cloth::failure_name(info.failure_reason)
         << "\",\n  \"message\":";
    json_string(file, info.message);
    file << ",\n  \"failed_step\":";
    if (info.committed) file << "null";
    else file << completed + 1;
    file << ",\n  \"residual_norm\":";
    json_number(file, info.residual_norm);
    file << ",\n  \"residual_threshold\":";
    json_number(file, info.residual_threshold);
    file << ",\n  \"last_attempt\":{\"iterations\":" << info.iterations
         << ",\"accepted_steps\":" << info.accepted_steps
         << ",\"line_search_rejections\":" << info.line_search_rejections
         << ",\"degenerate_rejections\":" << info.degenerate_rejections
         << ",\"shifted_steps\":" << info.shifted_steps
         << ",\"start_fallback\":" << info.start_fallback << "},\n";
    cloth_cli::summary(file,info.contact);
    file << "\n}\n";
    file.close();
}

} // namespace

int main(int argc, char** argv) {
    try {
        const auto args = parse(argc, argv);
        if (args.help) { usage(); return 0; }
        auto mesh = cloth::make_rectangular_mesh(2, 1, 1.0, 0.7);
        cloth::SimulationState state{mesh.positions,
            std::vector<cloth::Vec3>(mesh.positions.size(), {0.0, 0.0, 0.0})};
        for (std::size_t vertex = 0; vertex < 3; ++vertex) {
            state.positions[vertex][0] = 1.06 * state.positions[vertex][0] + 0.025;
            state.positions[vertex][1] = -0.025 * (1.0 + mesh.reference_positions[vertex][0]);
            state.positions[vertex][2] -= 0.025;
            state.velocities[vertex] = {0.015, -0.02, 0.01};
        }
        const cloth::Material material{2.0e4, 0.3, 0.001, 0.05};
        cloth::NewtonOptions options;
        options.max_iterations = args.max_iterations;
        options.residual_tolerance = 1.0e-10;
        const cloth::ImplicitEulerSolver solver(mesh, material, 800.0, args.mass_mode,
                                                {0.0, -9.81, 0.0}, options, args.ground);
        std::filesystem::create_directories(args.output);
        write_config(args, mesh, material, state, options);
        auto csv = output_file(args.output / "trajectory.csv");
        csv << "step,time,vertex,x,y,z,vx,vy,vz\n";
        write_frame(csv, 0, 0.0, state);
        auto contact_csv=output_file(args.output/"contact.csv");cloth_cli::header(contact_csv);
        std::size_t completed = 0;
        cloth::StepInfo info;
        while (completed < args.steps) {
            auto result = solver.step(state, {3, 4, 5}, args.dt);
            info = std::move(result.info);
            if (!info.committed) break;
            state = std::move(result.state);
            ++completed;
            write_frame(csv, completed, static_cast<double>(completed) * args.dt, state);
            cloth_cli::row(contact_csv,completed,info.contact);
        }
        csv.close();contact_csv.close();
        write_summary(args, completed, info);
        std::cout << "steps_completed=" << completed << '/' << args.steps
                  << " converged=" << std::boolalpha << info.solver_converged
                  << " failure_reason=" << cloth::failure_name(info.failure_reason)
                  << " output=" << args.output.string() << '\n';
        return completed == args.steps ? 0 : 2;
    } catch (const std::exception& error) {
        std::cerr << "cloth_simulate: " << error.what() << '\n';
        return 1;
    }
}
