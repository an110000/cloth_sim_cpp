#include "cloth/xpbd.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void check(bool condition,const char* text) {if(!condition) throw std::runtime_error(text);}
void near(double a,double b,double tol,const char* text) {check(std::isfinite(a)&&std::abs(a-b)<=tol,text);}
void split_and_rollback() {
    const cloth::SimulationState state{{{0,.2,0},{.3,-.2,0}},{{0,0,0},{0,0,0}}};
    const std::vector<cloth::Constraint> constraints{
        cloth::DistanceConstraint{0,1,.5,0,0,cloth::ConstraintGroup::structural}};
    cloth::GroundContactOptions ground;ground.enabled=true;
    cloth::XpbdOptions options;options.iterations=1;
    const cloth::XpbdSolver budget({1,1},constraints,{0,0,0},options,ground);
    const auto result=budget.step(state,{0},.01);
    check(result.info.committed && result.info.pre_contact_projection_converged &&
          !result.info.projection_converged,"Final residual must be measured after contact");
    near(result.info.pre_contact_groups[0].max_compliant_residual,0,1e-14,"Pre-contact R");
    near(result.info.groups[0].max_compliant_residual,.5-std::sqrt(.13),1e-14,"Independent post-contact R");
    check(result.state.positions[0]==state.positions[0],"Fixed anchor cannot be projected");
    near(result.state.positions[1][1],0,0,"Free penetration removed");
    options.require_convergence=true;
    const cloth::XpbdSolver strict({1,1},constraints,{0,0,0},options,ground);
    const auto rejected=strict.step(state,{0},.01);
    check(!rejected.info.committed && rejected.info.failure_reason==cloth::XpbdFailure::projection_not_converged &&
          rejected.state.positions==state.positions && rejected.state.velocities==state.velocities,
          "Strict mode rejects post-contact material residual and rolls back");
    check(rejected.info.contact.max_penetration_before>.19,"Rejected candidate diagnostics retained");
    ground.enabled=false;
    const cloth::XpbdSolver legacy({1,1},constraints,{0,0,0},options);
    const cloth::XpbdSolver disabled({1,1},constraints,{0,0,0},options,ground);
    check(legacy.step(state,{0},.01).state.positions==disabled.step(state,{0},.01).state.positions,
          "Disabled contact preserves trajectory");
}
void both_models() {
    const auto mesh=cloth::make_rectangular_mesh(1,1,1,1);
    const auto masses=cloth::lumped_masses(mesh,800,.001);
    const cloth::Material material{20000,.3,.001,.05};
    cloth::GroundContactOptions ground;ground.enabled=true;ground.friction.kinetic=.3;
    cloth::SimulationState state{mesh.positions,std::vector<cloth::Vec3>(4,{1,-1,0})};
    for(auto& x:state.positions) x[1]=-.01;
    for(bool energy:{false,true}) {
        const auto constraints=energy?cloth::build_neohookean_constraints(mesh,material,1e-8,1)
                                     :cloth::build_rectangular_constraints(mesh,1,1,1e-8,1e-8,1e-5,true);
        const cloth::XpbdSolver solver(masses,constraints,{0,0,0},{},ground);
        const auto result=solver.step(state,{},.01);
        check(result.info.committed && result.info.projection_converged,"Uniform plane impact in both material models");
        check(result.info.contact.active_count==4 && result.info.contact.sliding_count==4,"Every node impacts and slides");
        near(result.info.contact.normal_impulse_total,.8,2e-12,"Actual nodal momentum impulse sum");
        near(result.info.contact.friction_dissipation,.5*.8*(1-.7*.7),2e-12,"Tangential exact kinetic loss");
        near(result.info.contact.friction_kinetic_change,-result.info.contact.friction_dissipation,1e-14,
             "Lumped model friction kinetic change");
        for(std::size_t i=0;i<4;++i) {
            near(result.state.positions[i][1],0,0,"Plane projection");
            near(result.state.velocities[i][0],.7,2e-12,"Coulomb sliding speed");
            near(result.state.velocities[i][1],0,1e-14,"Inelastic normal impact");
        }
    }
}
void collapsed_projection() {
    cloth::GroundContactOptions ground;ground.enabled=true;
    const cloth::SimulationState state{{{0,-1,0},{0,-2,0}},{{0,0,0},{0,0,0}}};
    const cloth::XpbdSolver solver({1,1},{cloth::DistanceConstraint{0,1,1,0,0,
                                  cloth::ConstraintGroup::structural}},{0,0,0},{},ground);
    const auto result=solver.step(state,{},.01);
    check(!result.info.committed && result.info.failure_reason==cloth::XpbdFailure::invalid_geometry &&
          result.state.positions==state.positions,"Contact-induced geometry collapse rolls back");
}
}
int main() {
    try {split_and_rollback();both_models();collapsed_projection();
         std::cout<<"XPBD contact integration passed\n";return 0;}
    catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
