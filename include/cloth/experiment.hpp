#pragma once

#include "cloth/state.hpp"

#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace cloth {

// These names intentionally match the engineering-test protocol outcomes.
enum class RunStatus {
    committed,
    failed_numerical,
    failed_memory,
    failed_environment,
    timeout,
    interrupted,
    unsupported,
    not_run,
};

const char* run_status_name(RunStatus status) noexcept;

struct RunMetadata {
    // "<unset>" is an explicit placeholder, not evidence that a value was
    // discovered. A launcher should replace it before publishing a result.
    std::string run_id = "<unset>";
    std::string source_sha256 = "<unset>";
    std::string config_sha256 = "<unset>";
    std::string scene_sha256 = "<unset>";
    std::string model_identity = "<unset>";
    std::string mass_mode = "<unset>";
    std::string bending_model = "<unset>";
    std::string precision = "<unset>";
    std::string backend = "<unset>";
    std::string device = "<unset>";
    std::string toolchain = "<unset>";
    std::string container_digest = "<unset>";
    std::string random_state = "<unset>";
    double dt_seconds = 0.0;
    std::size_t requested_steps = 0;
    std::size_t committed_steps = 0;
    std::string solver_budget = "<unset>";
    std::string tolerances = "<unset>";
};

struct ExperimentFrame {
    std::size_t step = 0; // Committed physical step; starts at one.
    double time_seconds = 0.0;
    SimulationState state;
};

struct ExperimentRecord {
    RunMetadata metadata;
    RunStatus status = RunStatus::not_run;
    std::string failure_reason;
    double wall_clock_seconds = 0.0;
    double solver_seconds = 0.0;
    std::vector<ExperimentFrame> committed_frames;
};

// Adds one and only one committed frame. Candidate states from a failed step
// must never be passed here. The next required step is previous+1; positions
// and velocities must have equal, nonzero size and contain finite values.
void append_committed_frame(ExperimentRecord& record, ExperimentFrame frame);

// Changes the run outcome. A successful committed run must have exactly the
// requested number of frames and an empty reason. Every non-success status
// requires a reason; its already committed prefix remains valid evidence.
void finish_experiment(ExperimentRecord& record, RunStatus status,
                       std::string reason = {});

// Checks all invariants without touching the filesystem. Throws
// invalid_argument for malformed records and domain_error for non-finite
// timing or frame values.
void validate_experiment_record(const ExperimentRecord& record);

// Writes a deterministic UTF-8 text report. Keys and frame order are fixed;
// floating values use max_digits10. The function first creates
// <path>.tmp beside the destination, flushes/closes it, then renames it. To
// preserve atomic delivery on Windows, an existing destination is rejected;
// callers can archive or remove it explicitly before a new run.
void write_experiment_text(const std::filesystem::path& path,
                           const ExperimentRecord& record);

} // namespace cloth
