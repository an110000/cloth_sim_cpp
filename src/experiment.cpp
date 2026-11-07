#include "cloth/experiment.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <utility>

namespace cloth {
namespace {

void require_nonempty(const std::string& value, const char* name) {
    if (value.empty()) throw std::invalid_argument(std::string(name)+" must not be empty");
}

void require_finite(double value, const char* name) {
    if (!std::isfinite(value)) throw std::domain_error(std::string(name)+" must be finite");
}

void write_string(std::ostream& output, const char* key, const std::string& value) {
    // The format is line-oriented; reject control characters instead of
    // silently producing ambiguous records.
    for (const unsigned char character : value) {
        if (character == '\n' || character == '\r' || character == '\t')
            throw std::invalid_argument(std::string(key)+" contains a control character");
    }
    output << key << '=' << value << '\n';
}

void write_number(std::ostream& output, const char* key, double value) {
    require_finite(value, key);
    output << key << '=' << std::setprecision(std::numeric_limits<double>::max_digits10)
           << value << '\n';
}

void write_integer(std::ostream& output, const char* key, std::size_t value) {
    output << key << '=' << value << '\n';
}

} // namespace

const char* run_status_name(RunStatus status) noexcept {
    switch (status) {
    case RunStatus::committed: return "committed";
    case RunStatus::failed_numerical: return "failed_numerical";
    case RunStatus::failed_memory: return "failed_memory";
    case RunStatus::failed_environment: return "failed_environment";
    case RunStatus::timeout: return "timeout";
    case RunStatus::interrupted: return "interrupted";
    case RunStatus::unsupported: return "unsupported";
    case RunStatus::not_run: return "not_run";
    }
    return "unknown";
}

void append_committed_frame(ExperimentRecord& record, ExperimentFrame frame) {
    if (record.status != RunStatus::not_run)
        throw std::logic_error("Cannot append a frame after a failed experiment");
    const std::size_t expected = record.committed_frames.size()+1;
    if (frame.step != expected) throw std::invalid_argument("Committed frame steps must be contiguous from one");
    require_finite(frame.time_seconds, "frame.time_seconds");
    if (frame.time_seconds < 0.0) throw std::invalid_argument("frame.time_seconds must be nonnegative");
    if (frame.state.positions.empty() || frame.state.positions.size() != frame.state.velocities.size())
        throw std::invalid_argument("Committed frame state sizes must match and be nonzero");
    for (const auto& point : frame.state.positions) for (const auto value : point)
        require_finite(value, "frame.position");
    for (const auto& point : frame.state.velocities) for (const auto value : point)
        require_finite(value, "frame.velocity");
    record.committed_frames.push_back(std::move(frame));
    record.metadata.committed_steps = record.committed_frames.size();
    if (record.status == RunStatus::committed) record.status = RunStatus::not_run;
}

void finish_experiment(ExperimentRecord& record, RunStatus status, std::string reason) {
    if (status == RunStatus::committed) {
        if (record.metadata.committed_steps != record.metadata.requested_steps)
            throw std::invalid_argument("A committed run must contain every requested physical step");
        if (!reason.empty()) throw std::invalid_argument("A committed run cannot have a failure reason");
    } else if (reason.empty()) {
        throw std::invalid_argument("A non-committed run requires a failure or availability reason");
    }
    record.status = status;
    record.failure_reason = std::move(reason);
}

void validate_experiment_record(const ExperimentRecord& record) {
    require_nonempty(record.metadata.run_id, "metadata.run_id");
    require_nonempty(record.metadata.model_identity, "metadata.model_identity");
    if (!std::isfinite(record.metadata.dt_seconds) || record.metadata.dt_seconds <= 0.0)
        throw std::invalid_argument("metadata.dt_seconds must be finite and positive");
    if (record.metadata.committed_steps != record.committed_frames.size())
        throw std::invalid_argument("metadata.committed_steps does not match committed frames");
    if (record.metadata.committed_steps > record.metadata.requested_steps)
        throw std::invalid_argument("committed steps exceed requested steps");
    if (!std::isfinite(record.wall_clock_seconds) || record.wall_clock_seconds < 0.0)
        throw std::domain_error("wall_clock_seconds must be finite and nonnegative");
    if (!std::isfinite(record.solver_seconds) || record.solver_seconds < 0.0)
        throw std::domain_error("solver_seconds must be finite and nonnegative");
    if (record.status == RunStatus::committed) {
        if (record.metadata.committed_steps != record.metadata.requested_steps || !record.failure_reason.empty())
            throw std::invalid_argument("committed status does not describe a complete run");
    } else if (record.failure_reason.empty()) {
        throw std::invalid_argument("non-committed status requires failure_reason");
    }
    std::size_t expected_step = 1;
    double previous_time = -1.0;
    for (const auto& frame : record.committed_frames) {
        if (frame.step != expected_step++) throw std::invalid_argument("Frame steps are not contiguous");
        require_finite(frame.time_seconds, "frame.time_seconds");
        if (frame.time_seconds < previous_time) throw std::invalid_argument("Frame time is not monotone");
        previous_time = frame.time_seconds;
        if (frame.state.positions.empty() || frame.state.positions.size() != frame.state.velocities.size())
            throw std::invalid_argument("Frame state sizes are invalid");
        for (const auto& point : frame.state.positions) for (const auto value : point)
            require_finite(value, "frame.position");
        for (const auto& point : frame.state.velocities) for (const auto value : point)
            require_finite(value, "frame.velocity");
    }
}

void write_experiment_text(const std::filesystem::path& path, const ExperimentRecord& record) {
    validate_experiment_record(record);
    if (path.empty()) throw std::invalid_argument("Experiment output path must not be empty");
    if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());
    const auto temporary = std::filesystem::path(path.string()+".tmp");
    std::error_code error;
    if (std::filesystem::exists(path, error)) throw std::filesystem::filesystem_error(
        "Refusing to replace an existing experiment report", path, std::make_error_code(std::errc::file_exists));
    if (std::filesystem::exists(temporary, error)) throw std::filesystem::filesystem_error(
        "Temporary experiment report already exists", temporary, std::make_error_code(std::errc::file_exists));
    try {
        std::ofstream output(temporary, std::ios::binary|std::ios::trunc);
        output.exceptions(std::ios::badbit|std::ios::failbit);
        output << "CLOTH_EXPERIMENT_TEXT 1\n";
        write_string(output,"status",run_status_name(record.status));
        write_string(output,"failure_reason",record.failure_reason);
        write_string(output,"run_id",record.metadata.run_id);
        write_string(output,"source_sha256",record.metadata.source_sha256);
        write_string(output,"config_sha256",record.metadata.config_sha256);
        write_string(output,"scene_sha256",record.metadata.scene_sha256);
        write_string(output,"model_identity",record.metadata.model_identity);
        write_string(output,"mass_mode",record.metadata.mass_mode);
        write_string(output,"bending_model",record.metadata.bending_model);
        write_string(output,"precision",record.metadata.precision);
        write_string(output,"backend",record.metadata.backend);
        write_string(output,"device",record.metadata.device);
        write_string(output,"toolchain",record.metadata.toolchain);
        write_string(output,"container_digest",record.metadata.container_digest);
        write_string(output,"random_state",record.metadata.random_state);
        write_string(output,"solver_budget",record.metadata.solver_budget);
        write_string(output,"tolerances",record.metadata.tolerances);
        write_number(output,"dt_seconds",record.metadata.dt_seconds);
        write_integer(output,"requested_steps",record.metadata.requested_steps);
        write_integer(output,"committed_steps",record.metadata.committed_steps);
        write_number(output,"wall_clock_seconds",record.wall_clock_seconds);
        write_number(output,"solver_seconds",record.solver_seconds);
        write_integer(output,"frame_count",record.committed_frames.size());
        output << "frames_begin\n";
        output << std::setprecision(std::numeric_limits<double>::max_digits10);
        for (const auto& frame : record.committed_frames) {
            output << "frame " << frame.step << ' ' << frame.time_seconds << ' '
                   << frame.state.positions.size() << '\n';
            for (std::size_t node=0; node<frame.state.positions.size(); ++node) {
                output << "node " << node;
                for (const auto value : frame.state.positions[node]) output << ' ' << value;
                for (const auto value : frame.state.velocities[node]) output << ' ' << value;
                output << '\n';
            }
        }
        output << "frames_end\n";
        output.flush();
        output.close();
        std::filesystem::rename(temporary, path);
    } catch (...) {
        std::error_code cleanup;
        std::filesystem::remove(temporary, cleanup);
        throw;
    }
}

} // namespace cloth
