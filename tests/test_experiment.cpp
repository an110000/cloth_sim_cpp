#include "cloth/experiment.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
void check(bool condition, const char* message) { if (!condition) throw std::runtime_error(message); }
cloth::ExperimentRecord complete() {
    cloth::ExperimentRecord record;
    record.metadata.run_id="run-06b";
    record.metadata.model_identity="fixed-teaching-scene";
    record.metadata.requested_steps=2;
    record.metadata.dt_seconds=.01;
    for (std::size_t step=1;step<=2;++step) {
        cloth::ExperimentFrame frame;
        frame.step=step; frame.time_seconds=step*.01;
        frame.state.positions={cloth::Vec3{double(step),0.,0.}};
        frame.state.velocities={cloth::Vec3{1.,0.,0.}};
        cloth::append_committed_frame(record,std::move(frame));
    }
    record.wall_clock_seconds=.12; record.solver_seconds=.08;
    cloth::finish_experiment(record,cloth::RunStatus::committed);
    return record;
}
void complete_and_replay() {
    auto record=complete();
    const auto status=std::string(cloth::run_status_name(record.status));
    check(status=="committed","Committed status name");
    const auto path=std::filesystem::current_path()/"cloth_experiment_contract_test.txt";
    std::error_code error; std::filesystem::remove(path,error); std::filesystem::remove(path.string()+".tmp",error);
    cloth::write_experiment_text(path,record);
    check(std::filesystem::exists(path),"Report was not delivered");
    std::ifstream first(path,std::ios::binary); const std::string bytes((std::istreambuf_iterator<char>(first)),{});
    check(bytes.find("status=committed\n")!=std::string::npos,"Status is serialized");
    check(bytes.find("committed_steps=2\n")!=std::string::npos,"Committed count is serialized");
    check(bytes.find("frame 2 0.02 1\n")!=std::string::npos,"Frame order/time is serialized");
    bool rejected=false; try { cloth::write_experiment_text(path,record); } catch (const std::filesystem::filesystem_error&) { rejected=true; }
    check(rejected,"Existing report must not be replaced by a non-atomic write");
    std::ifstream second(path,std::ios::binary); const std::string replay((std::istreambuf_iterator<char>(second)),{});
    check(bytes==replay,"Existing report changed after rejected replay");
    check(!std::filesystem::exists(path.string()+".tmp"),"Temporary file leaked");
    std::filesystem::remove(path,error);
}
void failure_prefix_and_invalid_states() {
    cloth::ExperimentRecord record;
    record.metadata.run_id="failed-run"; record.metadata.model_identity="scene";
    record.metadata.requested_steps=3; record.metadata.dt_seconds=.01;
    cloth::ExperimentFrame frame; frame.step=1; frame.time_seconds=.01;
    frame.state.positions={cloth::Vec3{0.,0.,0.}}; frame.state.velocities={cloth::Vec3{0.,0.,0.}};
    cloth::append_committed_frame(record,std::move(frame));
    cloth::finish_experiment(record,cloth::RunStatus::failed_numerical,"Newton iteration limit");
    check(record.metadata.committed_steps==1 && record.committed_frames.size()==1,"Failure preserves committed prefix");
    bool append_rejected=false; try {
        cloth::ExperimentFrame candidate; candidate.step=2; candidate.time_seconds=.02;
        candidate.state.positions={cloth::Vec3{9.,0.,0.}}; candidate.state.velocities={cloth::Vec3{0.,0.,0.}};
        cloth::append_committed_frame(record,std::move(candidate));
    } catch (const std::logic_error&) { append_rejected=true; }
    check(append_rejected,"Failed run accepted a candidate frame");
    bool complete_rejected=false; try {
        cloth::ExperimentRecord incomplete; incomplete.metadata.run_id="incomplete";
        incomplete.metadata.model_identity="scene"; incomplete.metadata.dt_seconds=.01; incomplete.metadata.requested_steps=1;
        cloth::finish_experiment(incomplete,cloth::RunStatus::committed);
    } catch (const std::invalid_argument&) { complete_rejected=true; }
    check(complete_rejected,"Incomplete run marked committed");
    bool nan_rejected=false; try {
        auto bad=complete(); bad.wall_clock_seconds=std::numeric_limits<double>::quiet_NaN();
        cloth::validate_experiment_record(bad);
    } catch (const std::domain_error&) { nan_rejected=true; }
    check(nan_rejected,"NaN timing accepted");
}
}
int main() {
    try { complete_and_replay(); failure_prefix_and_invalid_states();
        std::cout<<"Experiment metadata, status, replay, failure prefix and atomic delivery passed\n"; return 0;
    } catch (const std::exception& error) { std::cerr<<error.what()<<'\n'; return 1; }
}
