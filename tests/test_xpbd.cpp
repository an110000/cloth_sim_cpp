#include "cloth/xpbd.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
void check(bool condition, const char* message) { if (!condition) throw std::runtime_error(message); }
void near(double a, double b, double tolerance, const char* message) {
    check(std::isfinite(a) && std::isfinite(b) && std::abs(a-b)<=tolerance, message);
}
cloth::SimulationState pair(double length) {
    return {{{0,0,0},{length,0,0}}, {{0,0,0},{0,0,0}}};
}
cloth::Constraint distance(double rest, double compliance=0) {
    return cloth::DistanceConstraint{0,1,rest,compliance,0,cloth::ConstraintGroup::structural};
}
void check_pair() {
    cloth::XpbdOptions options;
    options.iterations=1;
    options.require_convergence=true;
    const cloth::XpbdSolver hard({2,3},{distance(1)}, {0,0,0},options);
    const auto input=pair(1.2);
    const auto result=hard.step(input,{0},.1);
    check(result.info.committed && result.info.projection_converged, "Hard one-fixed pair");
    near(result.state.positions[1][0],1,1e-14,"Full correction uses zero inverse fixed mass");
    near(result.state.velocities[1][0],-2,1e-13,"Velocity reconstruction");
    check(input.positions==pair(1.2).positions,"Input state was changed");
    const cloth::XpbdSolver soft({1,1},{distance(1,.001)}, {0,0,0},options);
    const auto compliant=soft.step(input,{},.1);
    check(compliant.info.committed && compliant.info.projection_converged,"Soft residual convergence");
    check(compliant.info.groups[0].max_violation>1e-3,"Soft solution intentionally has nonzero C");
    near(compliant.info.groups[0].max_compliant_residual,0,1e-14,"Compliant residual, not raw C");
    const auto again=soft.step(input,{},.1);
    check(again.state.positions==compliant.state.positions && again.info.multipliers==compliant.info.multipliers,
          "Multipliers reset between independent time steps");
}
void rollback_and_status() {
    auto a=distance(1), b=distance(2);
    std::get<cloth::DistanceConstraint>(b).group=cloth::ConstraintGroup::shear;
    cloth::XpbdOptions options; options.iterations=1;
    const auto input=pair(1.4);
    const cloth::XpbdSolver budget({1,1},{a,b},{0,0,0},options);
    const auto approximation=budget.step(input,{},.01);
    check(approximation.info.finite_state && approximation.info.committed &&
          !approximation.info.projection_converged,"Budgeted finite approximation is not a convergence claim");
    options.require_convergence=true;
    const cloth::XpbdSolver strict({1,1},{a,b},{0,0,0},options);
    const auto failure=strict.step(input,{},.01);
    check(!failure.info.committed && failure.info.failure_reason==cloth::XpbdFailure::projection_not_converged,
          "Strict rejection");
    check(failure.state.positions==input.positions && failure.state.velocities==input.velocities,"Strict rollback");
    const auto invalid=strict.step(pair(0),{},.01);
    check(!invalid.info.committed && invalid.info.finite_state &&
          invalid.info.failure_reason==cloth::XpbdFailure::invalid_geometry,"Finite collapsed pair is not valid");
    const cloth::XpbdSolver fixed({1,1},{distance(1)},{0,-9.81,0});
    const auto frozen=fixed.step(pair(1),{1,0},.01);
    check(frozen.info.committed && frozen.info.projection_converged &&
          frozen.state.positions==pair(1).positions,"All fixed at rest");
    auto enormous=pair(1); for (auto& v:enormous.velocities) v[0]=std::numeric_limits<double>::max();
    const auto overflow=fixed.step(enormous,{},2);
    check(!overflow.info.committed && !overflow.info.projection_converged &&
          overflow.info.failure_reason==cloth::XpbdFailure::non_finite_arithmetic,
          "Prediction overflow rolls back");
}
void ballistic() {
    const auto mesh=cloth::make_rectangular_mesh(2,1,1,.7);
    const auto constraints=cloth::build_rectangular_constraints(mesh,2,1);
    const auto masses=cloth::lumped_masses(mesh,800,.001);
    const cloth::Vec3 g{.2,-9.81,.1}, v0{.1,.2,-.03};
    const cloth::XpbdSolver solver(masses,constraints,g);
    cloth::SimulationState state{mesh.positions,std::vector<cloth::Vec3>(6,v0)};
    for (int k=1;k<=8;++k) {
        const auto result=solver.step(state,{},.01);
        check(result.info.committed && result.info.projection_converged,"Free fall preserves rest constraints");
        state=result.state;
        for (std::size_t i=0;i<6;++i) for (std::size_t axis=0;axis<3;++axis) {
            near(state.positions[i][axis],mesh.positions[i][axis]+k*.01*v0[axis]+.5*k*(k+1)*.0001*g[axis],
                 3e-12,"Discrete free fall positions");
            near(state.velocities[i][axis],v0[axis]+k*.01*g[axis],3e-10,"Discrete free fall velocities");
        }
    }
}
template<class F> void rejects(F f) {
    bool rejected=false; try { f(); } catch(const std::invalid_argument&) { rejected=true; }
    check(rejected,"Invalid input accepted");
}
void invalid_input() {
    const cloth::XpbdSolver solver({1,1},{distance(1)});
    rejects([&]{(void)solver.step(pair(1),{},0);});
    rejects([&]{(void)solver.step(pair(1),{0,0},.01);});
    rejects([&]{(void)solver.step(pair(1),{2},.01);});
    rejects([&]{cloth::XpbdSolver bad({1,0},{distance(1)});});
    rejects([&]{cloth::XpbdOptions o; o.iterations=0; cloth::XpbdSolver bad({1,1},{distance(1)},{0,0,0},o);});
}
}
int main() {
    try { check_pair(); rollback_and_status(); ballistic(); invalid_input();
        std::cout<<"XPBD prediction, compliance, fixed nodes, resets, diagnostics and rollback passed\n"; return 0;
    } catch (const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
