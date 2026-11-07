#include "cloth/contact.hpp"
#include "cloth/implicit.hpp"
#include "cloth/xpbd.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
class Reader {
public:
    explicit Reader(const std::string& path) {
        std::ifstream input(path);
        check(static_cast<bool>(input), "Cannot open task-05B fixture");
        std::string line, all;
        while (std::getline(input, line)) if (!line.empty() && line.front() != '#') all += line+'\n';
        stream_.str(all);
    }
    template<class T> T get() {
        T value{};
        if (!(stream_ >> value)) throw std::runtime_error("Truncated task-05B fixture");
        return value;
    }
    void expect(const std::string& expected) {
        const auto actual = get<std::string>();
        check(actual == expected, "Expected "+expected+", got "+actual);
    }
    void end() { std::string extra; check(!(stream_ >> extra), "Trailing task-05B data"); }
private:
    std::istringstream stream_;
};
void near(double actual, double expected, double absolute, const std::string& label) {
    const double tolerance = absolute+1.e-7*std::abs(expected);
    if (!std::isfinite(actual) || !std::isfinite(expected) || std::abs(actual-expected) > tolerance) {
        std::ostringstream error;
        error.precision(17);
        error << label << ": actual=" << actual << " expected=" << expected
              << " error=" << std::abs(actual-expected) << " tolerance=" << tolerance;
        throw std::runtime_error(error.str());
    }
}
double scalar(Reader& reader, const std::string& tag) {
    reader.expect(tag); return reader.get<double>();
}
std::size_t count(Reader& reader, const std::string& tag, std::size_t maximum=1000) {
    reader.expect(tag);
    const auto value = reader.get<std::size_t>();
    check(value <= maximum, "Invalid count for "+tag);
    return value;
}
std::vector<double> scalars(Reader& reader, const std::string& tag, std::size_t size) {
    reader.expect(tag);
    std::vector<double> result(size);
    for (auto& value : result) value = reader.get<double>();
    return result;
}
std::vector<cloth::Vec3> points(Reader& reader, const std::string& tag, std::size_t size) {
    reader.expect(tag);
    std::vector<cloth::Vec3> result(size);
    for (auto& point : result) for (auto& value : point) value = reader.get<double>();
    return result;
}
std::vector<std::size_t> fixed(Reader& reader, std::size_t maximum) {
    const auto size = count(reader, "fixed", maximum);
    std::vector<std::size_t> result(size);
    for (auto& node : result) node = reader.get<std::size_t>();
    return result;
}
void compare(const std::vector<double>& actual, const std::vector<double>& expected,
             double absolute, const std::string& label) {
    check(actual.size() == expected.size(), label+" count");
    for (std::size_t i=0; i<actual.size(); ++i) near(actual[i], expected[i], absolute, label+"["+std::to_string(i)+"]");
}
void compare(const std::vector<cloth::Vec3>& actual, const std::vector<cloth::Vec3>& expected,
             double absolute, const std::string& label) {
    check(actual.size() == expected.size(), label+" count");
    for (std::size_t i=0; i<actual.size(); ++i) for (std::size_t axis=0; axis<3; ++axis)
        near(actual[i][axis], expected[i][axis], absolute, label+"["+std::to_string(i)+"]["+std::to_string(axis)+"]");
}
cloth::FrictionOptions friction_options(Reader& reader) {
    const auto data = scalars(reader, "friction_options", 2);
    return {data[0], data[1]};
}
struct SourceInfo {
    std::size_t active, sticking, sliding;
    double penetration, normal, tangent, dissipation;
};
SourceInfo source_info(Reader& reader, const std::string& prefix) {
    return {count(reader,prefix+"_active"), count(reader,prefix+"_sticking"), count(reader,prefix+"_sliding"),
            scalar(reader,prefix+"_max_penetration"), scalar(reader,prefix+"_normal_total"),
            scalar(reader,prefix+"_tangent_total"), scalar(reader,prefix+"_dissipation")};
}
void compare_info(const cloth::ContactInfo& actual, const SourceInfo& source, bool compare_counts,
                  double tolerance, const std::string& label) {
    check(actual.active_count == source.active, label+" active count");
    if (compare_counts) {
        check(actual.sticking_count == source.sticking, label+" sticking count");
        check(actual.sliding_count == source.sliding, label+" sliding count");
    }
    near(actual.max_penetration_before, source.penetration, tolerance, label+" penetration");
    near(actual.normal_impulse_total, source.normal, tolerance, label+" normal impulse");
    near(actual.tangential_impulse_total, source.tangent, tolerance, label+" tangent impulse");
    near(actual.friction_dissipation, source.dissipation, tolerance, label+" friction dissipation");
}

void friction_cases(Reader& reader) {
    check(count(reader,"friction_cases") == 4, "Expected four friction cases");
    std::set<std::string> names{"sticking","sliding","zero_normal","zero_tangent"};
    for (std::size_t index=0; index<4; ++index) {
        reader.expect("friction"); const auto name=reader.get<std::string>();
        check(names.erase(name)==1, "Unknown friction case");
        const auto n=count(reader,"vertices"); check(n==1, "Single-node friction fixture expected");
        const auto masses=scalars(reader,"masses",n), normal=scalars(reader,"normal_impulses",n);
        const auto options=friction_options(reader);
        const auto before=points(reader,"velocities",n), after=points(reader,"expected_velocities",n);
        const auto tangent=scalars(reader,"tangential_impulses",n), dissipation=scalars(reader,"dissipation_per_node",n);
        const auto sticking=count(reader,"source_sticking"), sliding=count(reader,"source_sliding");
        reader.expect("end_friction");
        const auto result=cloth::apply_coulomb_friction(before,normal,masses,{},options);
        compare(result.velocities,after,1.e-12,name+" friction velocity");
        compare(result.info.tangential_impulses,tangent,1.e-12,name+" friction impulse");
        compare(result.info.friction_dissipation_per_node,dissipation,1.e-12,name+" friction energy");
        if (name!="zero_tangent") {
            check(result.info.sticking_count==sticking && result.info.sliding_count==sliding,name+" friction state count");
        }
        near(result.info.friction_kinetic_change,-result.info.friction_dissipation,1.e-14,name+" actual kinetic loss");
    }
}

void ground_cases(Reader& reader) {
    check(count(reader,"ground_cases")==5,"Expected five ground cases");
    std::set<std::string> names{"incoming_sliding","incoming_sticking","plane_incoming","separated","zero_tangent"};
    for (std::size_t index=0; index<5; ++index) {
        reader.expect("ground"); const auto name=reader.get<std::string>();
        check(names.erase(name)==1,"Unknown ground case");
        const auto n=count(reader,"vertices"); check(n==1,"Single-node ground fixture expected");
        const auto masses=scalars(reader,"masses",n);
        const auto dt=scalar(reader,"dt");
        cloth::GroundContactOptions options; options.enabled=true;
        options.ground_y=scalar(reader,"ground_y"); options.friction=friction_options(reader);
        options.restitution=scalar(reader,"restitution"); options.stabilization_beta=scalar(reader,"beta");
        const cloth::SimulationState state{points(reader,"positions",n),points(reader,"velocities",n)};
        const cloth::SimulationState expected{points(reader,"expected_positions",n),points(reader,"expected_velocities",n)};
        const auto compare_counts=count(reader,"compare_counts",1)!=0;
        const auto source=source_info(reader,"source"); reader.expect("end_ground");
        const auto result=cloth::apply_ground_contact(state,masses,{},dt,options);
        compare(result.state.positions,expected.positions,1.e-12,name+" projected position");
        compare(result.state.velocities,expected.velocities,1.e-12,name+" projected velocity");
        compare_info(result.info,source,compare_counts,1.e-12,name);
        near(result.info.normal_impulse_total,masses[0]*(result.state.velocities[0][1]-state.velocities[0][1]),
             1.e-14,name+" actual normal momentum");
        near(result.info.max_penetration_after,0.,1.e-14,name+" no final penetration");
    }
}

void penalty_cases(Reader& reader) {
    check(count(reader,"penalty_cases")==2,"Expected two penalty cases");
    for (std::size_t index=0; index<2; ++index) {
        reader.expect("penalty"); const auto name=reader.get<std::string>();
        check(name==(index==0?"mixed":"separated"),"Penalty fixture order");
        const auto n=count(reader,"vertices"); check(n==3,"Three-node penalty fixture expected");
        const auto masses=scalars(reader,"masses",n); const auto anchors=fixed(reader,n);
        const auto dt=scalar(reader,"dt");
        cloth::GroundContactOptions options; options.enabled=true;
        options.ground_y=scalar(reader,"ground_y"); options.penalty_eta=scalar(reader,"eta");
        const auto x=points(reader,"positions",n); const auto energy=scalar(reader,"energy");
        const auto gradient=scalars(reader,"gradient_y",n), diagonal=scalars(reader,"diagonal_y",n);
        const auto impulses=scalars(reader,"normal_impulses",n), pen=scalars(reader,"penetration",n);
        reader.expect("end_penalty");
        const auto actual=cloth::evaluate_ground_penalty(x,masses,anchors,dt,options);
        near(actual.energy,energy,1.e-11,name+" penalty energy");
        compare(actual.gradient_y,gradient,1.e-10,name+" penalty gradient");
        compare(actual.diagonal_y,diagonal,1.e-8,name+" penalty diagonal");
        compare(actual.normal_impulses,impulses,1.e-11,name+" penalty Jn");
        compare(actual.penetration,pen,1.e-14,name+" penalty penetration");
    }
}

void legacy_audits(Reader& reader) {
    check(count(reader,"legacy_audits")==3,"Expected three legacy audits");
    std::set<std::string> names{"fixed_below_ground","separating_velocity","beta_reported_impulse"};
    for (std::size_t index=0; index<3; ++index) {
        reader.expect("audit"); const auto name=reader.get<std::string>();
        check(names.erase(name)==1,"Unknown legacy audit");
        const auto n=count(reader,"vertices"); check(n==1,"Single-node audit expected");
        const auto masses=scalars(reader,"masses",n); const auto anchors=fixed(reader,n);
        const auto dt=scalar(reader,"dt");
        cloth::GroundContactOptions options; options.enabled=true;
        options.ground_y=scalar(reader,"ground_y"); options.friction=friction_options(reader);
        options.restitution=scalar(reader,"restitution"); options.stabilization_beta=scalar(reader,"beta");
        const cloth::SimulationState state{points(reader,"positions",n),points(reader,"velocities",n)};
        const auto old_x=points(reader,"legacy_positions",n), old_v=points(reader,"legacy_velocities",n);
        const auto legacy=source_info(reader,"legacy");
        const auto old_momentum_change=scalar(reader,"actual_normal_momentum_change");
        reader.expect("end_audit");
        const auto port=cloth::apply_ground_contact(state,masses,anchors,dt,options);
        // These are explicit corrected contracts, never equality targets for the old bugs.
        if (name=="fixed_below_ground") {
            check(port.state.positions==state.positions && old_x!=state.positions,"Fixed projection bug audit");
            check(port.info.fixed_penetration_count==1 && port.info.active_count==0,"Fixed penetration report");
            near(port.info.normal_impulse_total,0.,0.,"Fixed vertex has no contact impulse");
        } else if (name=="separating_velocity") {
            check(port.state.velocities[0][1]==state.velocities[0][1] && old_v[0][1]==0.,"Outgoing normal velocity audit");
            near(port.info.normal_impulse_total,0.,0.,"Outgoing vertex gets no impulse");
        } else {
            check(std::abs(legacy.normal-old_momentum_change)>.1,"Legacy beta mismatch must be visible");
            near(port.info.normal_impulse_total,masses[0]*(port.state.velocities[0][1]-state.velocities[0][1]),
                 1.e-13,"Corrected beta impulse equals actual momentum change");
        }
    }
}

cloth::Mesh read_mesh(Reader& reader) {
    const auto n=count(reader,"vertices"), elements=count(reader,"triangles");
    check(n==4 && elements==2,"Expected one-cell reference mesh");
    cloth::Mesh mesh;
    mesh.reference_positions.resize(n); mesh.triangles.resize(elements);
    reader.expect("reference");
    for (auto& point:mesh.reference_positions) for (auto& value:point) value=reader.get<double>();
    mesh.positions=points(reader,"rest_positions",n);
    reader.expect("topology");
    for (auto& triangle:mesh.triangles) for (auto& node:triangle) node=reader.get<std::size_t>();
    return mesh;
}
std::vector<double> flatten(const std::vector<cloth::Vec3>& points) {
    std::vector<double> result;
    for (const auto& point:points) for (const auto value:point) result.push_back(value);
    return result;
}
double kinetic_change(const cloth::SparseMatrix& mass,const std::vector<cloth::Vec3>& before,
                      const std::vector<cloth::Vec3>& after) {
    const auto a=flatten(before), b=flatten(after);
    std::vector<double> delta(a.size()), average(a.size());
    for (std::size_t i=0;i<a.size();++i) { delta[i]=b[i]-a[i]; average[i]=.5*(a[i]+b[i]); }
    const auto momentum=mass.multiply(average);
    long double change=0.;
    for (std::size_t i=0;i<a.size();++i) change+=static_cast<long double>(delta[i])*momentum[i];
    return static_cast<double>(change);
}
double independent_implicit_residual(const cloth::ImplicitEulerSolver& solver,const cloth::Mesh& mesh,
    const cloth::Material& material,const cloth::SimulationState& previous,const cloth::SimulationState& current,
    const std::vector<std::size_t>& anchors,double dt,const cloth::GroundPenaltyEvaluation& penalty) {
    auto dx=flatten(current.positions); const auto x0=flatten(previous.positions),v0=flatten(previous.velocities);
    for (std::size_t i=0;i<dx.size();++i) dx[i]-=x0[i]+dt*v0[i];
    const auto inertia=solver.mass_matrix().multiply(dx);
    std::vector<double> gradient(dx.size(),0.);
    for (const auto& triangle:mesh.triangles) {
        std::array<cloth::Vec2,3> X{}; std::array<cloth::Vec3,3> x{};
        for (std::size_t local=0;local<3;++local) { X[local]=mesh.reference_positions[triangle[local]];x[local]=current.positions[triangle[local]]; }
        const auto element=cloth::evaluate_triangle(X,x,material);
        for (std::size_t local=0;local<3;++local) for (std::size_t axis=0;axis<3;++axis)
            gradient[3*triangle[local]+axis]+=element.gradient[local][axis];
    }
    long double squared=0.;
    for (std::size_t i=0;i<dx.size();++i) {
        if (std::find(anchors.begin(),anchors.end(),i/3)!=anchors.end()) continue;
        const double gy=i%3==1?penalty.gradient_y[i/3]:0.;
        const double residual=inertia[i]+dt*dt*(gradient[i]+gy-solver.external_force()[i]);
        squared+=static_cast<long double>(residual)*residual;
    }
    return std::sqrt(static_cast<double>(squared));
}

void trajectories(Reader& reader) {
    check(count(reader,"trajectories")==4,"Expected four contact trajectories");
    std::set<std::string> names{"implicit_lumped","implicit_consistent","distance","energy"};
    std::size_t total_steps=0;
    for (std::size_t trajectory=0;trajectory<4;++trajectory) {
        reader.expect("trajectory"); const auto name=reader.get<std::string>();
        check(names.erase(name)==1,"Unknown contact trajectory");
        const auto mesh=read_mesh(reader); const auto material_data=scalars(reader,"material",4);
        const cloth::Material material{material_data[0],material_data[1],material_data[2],material_data[3]};
        const auto density=scalar(reader,"density");
        const auto physical=scalars(reader,"physical_masses",4), original=scalars(reader,"original_node_masses",4);
        compare(cloth::lumped_masses(mesh,density,material.thickness),physical,1.e-13,name+" physical row sums");
        const auto original_eta=scalar(reader,"original_eta"),port_eta=scalar(reader,"port_eta");
        cloth::GroundContactOptions ground; ground.enabled=true;ground.penalty_eta=port_eta;
        ground.friction=friction_options(reader);ground.ground_y=scalar(reader,"ground_y");
        ground.restitution=scalar(reader,"restitution");ground.stabilization_beta=scalar(reader,"beta");
        const auto dt=scalar(reader,"dt"); const auto iterations=count(reader,"iterations");
        const auto newton_tolerance=scalar(reader,"newton_tolerance");
        const auto compliances=scalars(reader,"distance_compliances",3),energy_options=scalars(reader,"energy_options",2);
        const auto gravity_data=scalars(reader,"gravity",3);
        const cloth::Vec3 gravity{gravity_data[0],gravity_data[1],gravity_data[2]};
        const auto anchors=fixed(reader,4);
        cloth::SimulationState state{points(reader,"positions",4),points(reader,"velocities",4)};
        const bool implicit=name.rfind("implicit_",0)==0,consistent=name=="implicit_consistent";
        if (consistent) {
            check(ground.friction.kinetic==0.,"Consistent mapped comparison requires zero friction");
            near(port_eta,original_eta/2.,0.,"Explicit consistent eta mapping");
            for (std::size_t node=0;node<4;++node) {
                near(physical[node],2.*original[node],1.e-14,"Original consistent diagonal differs from physical row sum");
                near(port_eta*physical[node],original_eta*original[node],1.e-13,"Mapped penalty stiffness");
            }
        } else { near(port_eta,original_eta,0.,name+" same penalty eta");compare(original,physical,1.e-14,name+" original physical masses"); }
        cloth::NewtonOptions newton;newton.max_iterations=iterations;newton.residual_tolerance=newton_tolerance;
        cloth::XpbdOptions xpbd;xpbd.iterations=iterations;
        std::vector<cloth::Constraint> constraints;
        if (!implicit) constraints=name=="energy"?cloth::build_neohookean_constraints(mesh,material,energy_options[0],energy_options[1]):
            cloth::build_rectangular_constraints(mesh,1,1,compliances[0],compliances[1],compliances[2]);
        std::unique_ptr<cloth::ImplicitEulerSolver> implicit_solver;
        std::unique_ptr<cloth::XpbdSolver> xpbd_solver,free_xpbd;
        if (implicit) implicit_solver=std::make_unique<cloth::ImplicitEulerSolver>(mesh,material,density,
            consistent?cloth::MassMode::consistent:cloth::MassMode::lumped,gravity,newton,ground);
        else {
            xpbd_solver=std::make_unique<cloth::XpbdSolver>(physical,constraints,gravity,xpbd,ground);
            free_xpbd=std::make_unique<cloth::XpbdSolver>(physical,constraints,gravity,xpbd);
        }
        const auto mass=cloth::assemble_mass_matrix(mesh,density,material.thickness,
            consistent?cloth::MassMode::consistent:cloth::MassMode::lumped);
        check(count(reader,"steps")==8,"Expected eight contact steps");
        std::size_t contacted=0;
        for (std::size_t step=1;step<=8;++step) {
            reader.expect("step");check(reader.get<std::size_t>()==step,"Wrong step order");
            const auto pre_x=points(reader,"pre_contact_positions",4),pre_v=points(reader,"pre_contact_velocities",4);
            const auto expected_x=points(reader,"positions",4),expected_v=points(reader,"velocities",4);
            const auto pen=scalars(reader,"penetration",4),impulses=scalars(reader,"normal_impulses",4);
            const auto penalty_energy=scalar(reader,"penalty_energy"),source_residual=scalar(reader,"residual_norm");
            const auto source_threshold=scalar(reader,"residual_threshold"),physical_change=scalar(reader,"physical_mass_kinetic_change");
            const auto n_multipliers=count(reader,"multipliers");
            std::vector<double> multipliers(n_multipliers);for(auto& value:multipliers)value=reader.get<double>();
            const auto raw_before=scalars(reader,"raw_before",4),compliant_before=scalars(reader,"compliant_before",4);
            const auto raw_after=scalars(reader,"raw_after",4),compliant_after=scalars(reader,"compliant_after",4);
            const auto source=source_info(reader,"source");reader.expect("end_step");
            const auto previous=state; const auto label=name+" step "+std::to_string(step);
            cloth::ContactInfo contact;
            std::vector<cloth::Vec3> actual_pre_v;
            if (implicit) {
                const auto result=implicit_solver->step(state,anchors,dt);
                check(result.info.committed && result.info.solver_converged,label+" implicit solve failed");
                check(source_residual<=source_threshold*(1.+1.e-8),label+" original independently checked residual");
                check(result.info.residual_norm<=result.info.residual_threshold,label+" port residual threshold");
                state=result.state;contact=result.info.contact;
                compare(state.positions,pre_x,3.e-8,label+" no extra implicit projection");
                actual_pre_v=state.velocities;
                for(std::size_t node=0;node<4;++node)for(std::size_t axis=0;axis<3;++axis)
                    actual_pre_v[node][axis]=(state.positions[node][axis]-previous.positions[node][axis])/dt;
                compare(actual_pre_v,pre_v,3.e-6,label+" reconstructed pre-friction velocity");
                const auto p=cloth::evaluate_ground_penalty(state.positions,physical,anchors,dt,ground);
                const auto independent=independent_implicit_residual(*implicit_solver,mesh,material,previous,state,anchors,dt,p);
                near(result.info.residual_norm,independent,1.e-12,label+" independent port stationarity");
                near(contact.penalty_energy,penalty_energy,2.e-8,label+" penalty energy");
                near(contact.friction_kinetic_change,kinetic_change(mass,actual_pre_v,state.velocities),1.e-12,
                     label+" actual full-M friction kinetic change");
                near(contact.friction_kinetic_change,physical_change,2.e-8,label+" reference full-M kinetic change");
            } else {
                const auto free=free_xpbd->step(state,anchors,dt);
                check(free.info.committed,label+" valid pre-contact XPBD");
                compare(free.state.positions,pre_x,3.e-8,label+" pre-contact positions");
                compare(free.state.velocities,pre_v,3.e-6,label+" pre-contact velocities");
                const auto result=xpbd_solver->step(state,anchors,dt);
                check(result.info.committed,label+" budgeted contact XPBD failed");
                actual_pre_v=free.state.velocities;state=result.state;contact=result.info.contact;
                compare(result.info.multipliers,multipliers,2.e-8,label+" XPBD multipliers");
                bool pre_converged=true,post_converged=true;
                for(std::size_t group=0;group<4;++group) {
                    near(result.info.pre_contact_groups[group].max_violation,raw_before[group],2.e-8,label+" pre raw");
                    near(result.info.pre_contact_groups[group].max_compliant_residual,compliant_before[group],2.e-8,label+" pre R");
                    near(result.info.groups[group].max_violation,raw_after[group],2.e-8,label+" post raw");
                    near(result.info.groups[group].max_compliant_residual,compliant_after[group],2.e-8,label+" post R");
                    const auto tolerance=group==3?xpbd.energy_tolerance:group==2?xpbd.angle_tolerance:xpbd.distance_tolerance;
                    pre_converged=pre_converged && result.info.pre_contact_groups[group].max_compliant_residual<=tolerance;
                    post_converged=post_converged && result.info.groups[group].max_compliant_residual<=tolerance;
                }
                check(result.info.pre_contact_projection_converged==pre_converged,label+" pre-contact convergence label");
                check(result.info.projection_converged==post_converged,label+" post-contact convergence label");
                near(contact.normal_kinetic_change+contact.friction_kinetic_change,physical_change,2.e-8,label+" total contact kinetic change");
            }
            compare(state.positions,expected_x,3.e-8,label+" positions");
            compare(state.velocities,expected_v,3.e-6,label+" velocities");
            compare(contact.penetration,pen,3.e-8,label+" penetration vector");
            compare(contact.normal_impulses,impulses,2.e-8,label+" normal impulse vector");
            bool moving_tangent=true;
            for(std::size_t node=0;node<4;++node)if(impulses[node]>1.e-14 && std::hypot(pre_v[node][0],pre_v[node][2])<=1.e-14)
                moving_tangent=false;
            compare_info(contact,source,moving_tangent,2.e-8,label);
            for(const auto node:anchors) {
                check(state.positions[node]==previous.positions[node],label+" fixed position drift");
                check(state.velocities[node]==cloth::Vec3{0.,0.,0.},label+" fixed velocity");
                check(state.positions[node][1]>ground.ground_y,label+" fixture fixed-node domain");
            }
            contacted+=contact.active_count!=0;++total_steps;
        }
        check(contacted==8,name+" must exercise actual contact at every reference step");
        reader.expect("end_trajectory");
    }
    check(names.empty() && total_steps==32,"Incomplete contact trajectory coverage");
}
} // namespace
int main(int argc,char** argv) {
    try {
        check(argc==2,"Expected task-05B fixture path");Reader reader(argv[1]);
        reader.expect("CLOTH_TASK05B_REFERENCE");check(reader.get<int>()==1,"Unsupported task-05B fixture version");
        friction_cases(reader);ground_cases(reader);penalty_cases(reader);legacy_audits(reader);trajectories(reader);
        reader.expect("end");reader.end();
        std::cout<<"4 friction / 5 ground / 2 penalty / 3 legacy audits / 32 contact trajectory steps passed\n";
        return 0;
    } catch(const std::exception& error) { std::cerr<<error.what()<<'\n';return 1; }
}
