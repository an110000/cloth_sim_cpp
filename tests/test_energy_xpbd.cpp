#include "cloth/xpbd.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

void check(bool value,const char* message) { if (!value) throw std::runtime_error(message); }
void near(double actual,double expected,double tolerance,const char* message) {
    check(std::isfinite(actual) && std::isfinite(expected) &&
          std::abs(actual-expected)<=tolerance*std::max(1.0,std::abs(expected)),message);
}
template<class Exception,class Function> void rejects(Function f,const char* message) {
    try { f(); } catch (const Exception&) { return; }
    throw std::runtime_error(message);
}
cloth::Material material() { return {2e4,.3,.001,.05}; }
cloth::Mesh deformed() {
    auto mesh=cloth::make_rectangular_mesh(1,1,1,.8);
    for (auto& point:mesh.positions) point={1.12*point[0]+.04*point[2], .03*point[0]-.02*point[2], .91*point[2]};
    return mesh;
}
cloth::SimulationState state_at(const cloth::Mesh& mesh) {
    return {mesh.positions,std::vector<cloth::Vec3>(mesh.positions.size(),{0,0,0})};
}

void value_gradient_and_builder() {
    const auto mesh=deformed();
    const auto constraints=cloth::build_neohookean_constraints(mesh,material(),1e-5,.7);
    check(constraints.size()==mesh.triangles.size(),"One energy constraint per triangle");
    check(static_cast<int>(cloth::ConstraintGroup::structural)==0 &&
          static_cast<int>(cloth::ConstraintGroup::shear)==1 &&
          static_cast<int>(cloth::ConstraintGroup::bending)==2 &&
          static_cast<int>(cloth::ConstraintGroup::neohookean_energy)==3,"Old group indices remain stable");
    for (std::size_t index=0;index<constraints.size();++index) {
        const auto& energy=std::get<cloth::NeoHookeanEnergyConstraint>(constraints[index]);
        check(energy.nodes==mesh.triangles[index],"Energy builder preserves triangle order");
        check(cloth::constraint_group(constraints[index])==cloth::ConstraintGroup::neohookean_energy,
              "Energy group is separate from geometry groups");
        std::array<cloth::Vec3,3> local{};
        for (std::size_t vertex=0;vertex<3;++vertex) {
            check(energy.reference[vertex]==mesh.reference_positions[energy.nodes[vertex]],"Reference data copied by local node order");
            local[vertex]=mesh.positions[energy.nodes[vertex]];
        }
        const auto element=cloth::evaluate_triangle(energy.reference,local,energy.material);
        const auto evaluated=cloth::evaluate_constraint(constraints[index],mesh.positions);
        check(evaluated.valid_geometry && evaluated.count==3,"Valid energy constraint has three nodes");
        near(evaluated.value,element.energy/.7,2e-14,"C is E divided by energy scale");
        for (std::size_t vertex=0;vertex<3;++vertex) for (std::size_t axis=0;axis<3;++axis) {
            near(evaluated.gradients[vertex][axis],element.gradient[vertex][axis]/.7,2e-14,"Constraint gradient uses the same scale");
            auto plus=mesh.positions,minus=mesh.positions;
            const auto node=energy.nodes[vertex];
            const double epsilon=1e-6;
            plus[node][axis]+=epsilon; minus[node][axis]-=epsilon;
            const double fd=(cloth::evaluate_constraint(constraints[index],plus).value-
                             cloth::evaluate_constraint(constraints[index],minus).value)/(2*epsilon);
            near(evaluated.gradients[vertex][axis],fd,3e-8,"Energy constraint centered finite difference");
        }
    }
    auto reversed=mesh;
    for (auto& triangle:reversed.triangles) std::swap(triangle[1],triangle[2]);
    const auto opposite=cloth::build_neohookean_constraints(reversed,material(),1e-5,.7);
    near(cloth::evaluate_constraint(opposite[0],reversed.positions).value,
         cloth::evaluate_constraint(constraints[0],mesh.positions).value,2e-13,"Energy supports reversed consistent reference winding");
    auto regularized=mesh;
    for (std::size_t i=0;i<regularized.positions.size();++i) {
        regularized.positions[i]={.025*mesh.reference_positions[i][0],0,mesh.reference_positions[i][1]};
    }
    const auto extension=cloth::build_neohookean_constraints(regularized,material());
    check(cloth::evaluate_constraint(extension[0],regularized.positions).valid_geometry,
          "Log continuation is valid and does not skip small-area elements");

    // Initial admissibility is relative det(C), not the bending area threshold.
    const auto tiny=cloth::make_rectangular_mesh(1,1,4e-7,4e-7);
    const auto tiny_constraints=cloth::build_neohookean_constraints(tiny,material());
    check(tiny_constraints.size()==2,"Small but valid reference elements use material-domain checks");
}

void scale_compensated_projection() {
    const auto mesh=deformed();
    auto original=cloth::build_neohookean_constraints(mesh,material(),1e-5,.7)[0];
    auto rescaled=cloth::build_neohookean_constraints(mesh,material(),1e-5/9.0,2.1)[0];
    auto x=mesh.positions,z=x;
    const std::vector<double> weights{1,.5,1.5,2};
    for (int iteration=0;iteration<8;++iteration) {
        const auto a=cloth::project_constraint(original,x,weights,.01);
        const auto b=cloth::project_constraint(rescaled,z,weights,.01);
        check(a.valid_geometry && b.valid_geometry && a.applied && b.applied,"Scale test stays above denominator guard");
        near(b.delta_lambda,3*a.delta_lambda,4e-12,"Multiplier scales inversely to scaled constraint");
        for (std::size_t vertex=0;vertex<x.size();++vertex) for (std::size_t axis=0;axis<3;++axis) {
            near(x[vertex][axis],z[vertex][axis],4e-12,"Scale/compliance compensation preserves positions");
        }
        near(cloth::constraint_multiplier(rescaled),3*cloth::constraint_multiplier(original),4e-12,"Accumulated multiplier compensation");
        near(cloth::evaluate_constraint(rescaled,z).value,
             cloth::evaluate_constraint(original,x).value/3,4e-12,"Scaled C remains dimensionless");
    }
}

void rest_fixed_reset_and_modes() {
    const auto mesh=cloth::make_rectangular_mesh(1,1,1,.8);
    auto constraints=cloth::build_neohookean_constraints(mesh,material(),0,1);
    auto positions=mesh.positions;
    const auto noop=cloth::project_constraint(constraints[0],positions,{0,0,0,0},.01);
    check(noop.valid_geometry && !noop.applied && positions==mesh.positions,
          "Zero-gradient/all-fixed hard energy constraint is valid, not degenerate");
    cloth::XpbdOptions legacy{1,1e-6,1e-6,true};
    check(legacy.require_convergence && legacy.energy_tolerance==1e-6,"Old option aggregate initializers keep their meaning");
    const cloth::XpbdSolver rest(cloth::lumped_masses(mesh,800,.001),constraints,{0,0,0},legacy);
    const auto initial=state_at(mesh);
    const auto step=rest.step(initial,{},.01);
    check(step.info.committed && step.info.projection_converged && step.info.groups[3].count==2,
          "Rest energy constraints pass the energy group residual");
    check(step.state.positions==initial.positions,"Rest state remains unchanged");
    const auto frozen=rest.step(initial,{3,1,0,2},.01);
    check(frozen.info.committed && frozen.info.projection_converged,"All-fixed rest energy constraints are valid");

    const auto distorted=deformed();
    auto input=state_at(distorted);
    cloth::XpbdOptions options;
    options.iterations=1;
    options.energy_tolerance=1e-30;
    const auto energy=cloth::build_neohookean_constraints(distorted,material(),1e-5,1);
    const auto masses=cloth::lumped_masses(distorted,800,.001);
    const cloth::XpbdSolver budget(masses,energy,{0,0,0},options);
    const auto approximation=budget.step(input,{2,3},.01);
    check(approximation.info.committed && approximation.info.finite_state &&
          !approximation.info.projection_converged,"Finite budget energy approximation is not residual convergence");
    check(approximation.info.groups[3].max_compliant_residual>options.energy_tolerance,
          "Energy convergence does not accidentally use distance or angle group");
    check(approximation.state.positions[2]==input.positions[2] && approximation.state.positions[3]==input.positions[3],
          "Energy projection preserves fixed nodes");
    const auto repeated=budget.step(input,{2,3},.01);
    check(repeated.state.positions==approximation.state.positions &&
          repeated.info.multipliers==approximation.info.multipliers,"Each step resets energy multipliers");
    for (const auto& constraint:budget.constraints()) near(cloth::constraint_multiplier(constraint),0,0,"Solver owns unchanged reset constraints");
    options.require_convergence=true;
    const cloth::XpbdSolver strict(masses,energy,{0,0,0},options);
    const auto rejected=strict.step(input,{2,3},.01);
    check(!rejected.info.committed && !rejected.info.projection_converged &&
          rejected.info.failure_reason==cloth::XpbdFailure::projection_not_converged,"Strict energy residual rejection");
    check(rejected.state.positions==input.positions && rejected.state.velocities==input.velocities,
          "Strict failure rolls back energy state");
}

void domain_overflow_and_invalid_parameters() {
    const auto mesh=cloth::make_rectangular_mesh(1,1,1,1);
    auto constraints=cloth::build_neohookean_constraints(mesh,material());
    const auto masses=cloth::lumped_masses(mesh,800,.001);
    const cloth::XpbdSolver solver(masses,constraints,{0,0,0});
    auto collapsed=state_at(mesh);
    collapsed.positions[1]=collapsed.positions[0];
    const auto failure=solver.step(collapsed,{},.01);
    check(!failure.info.committed && failure.info.finite_state &&
          failure.info.failure_reason==cloth::XpbdFailure::invalid_geometry,"Energy domain failure refuses the whole step");
    check(failure.state.positions==collapsed.positions,"Invalid energy state remains unchanged");
    std::get<cloth::NeoHookeanEnergyConstraint>(constraints[0]).lambda_value=2;
    auto old_positions=collapsed.positions;
    const auto invalid=cloth::project_constraint(constraints[0],collapsed.positions,{1,1,1,1},.01);
    check(!invalid.valid_geometry && !invalid.applied && std::isinf(invalid.value_before),"Collapsed energy is not a zero constraint");
    check(collapsed.positions==old_positions && cloth::constraint_multiplier(constraints[0])==2,
          "Invalid local energy projection preserves positions and multiplier");

    const double scale=std::numeric_limits<double>::denorm_min();
    auto tiny_scale_hard=cloth::build_neohookean_constraints(mesh,material(),0,scale);
    auto exact_rest=mesh.positions;
    const auto rest_evaluation=cloth::evaluate_constraint(tiny_scale_hard[0],exact_rest);
    check(rest_evaluation.valid_geometry && rest_evaluation.value==0,"Tiny positive scale preserves exact unit reference energy");
    const auto rest_projection=cloth::project_constraint(tiny_scale_hard[0],exact_rest,{1,1,1,1},.01);
    check(rest_projection.valid_geometry && !rest_projection.applied && exact_rest==mesh.positions,
          "Direct g/scale keeps zero gradient finite without reciprocal overflow");
    auto tiny_scale=cloth::build_neohookean_constraints(mesh,material(),1e-8,scale);
    const cloth::XpbdSolver overflow_solver(masses,tiny_scale,{0,0,0});
    const auto input=state_at(deformed());
    const auto overflow=overflow_solver.step(input,{},.01);
    check(!overflow.info.committed && overflow.info.failure_reason==cloth::XpbdFailure::non_finite_arithmetic,
          "Overflowed E/scale is a numeric failure with rollback");
    check(overflow.state.positions==input.positions && overflow.state.velocities==input.velocities,
          "Scaled energy overflow keeps original state");
    auto local=input.positions;
    const auto local_before=local;
    rejects<std::overflow_error>([&]{ (void)cloth::project_constraint(tiny_scale[0],local,{1,1,1,1},.01); },
                                 "Scaled energy overflow propagates before local writes");
    check(local==local_before && cloth::constraint_multiplier(tiny_scale[0])==0,
          "Local numeric failure does not partly update state");
    rejects<std::invalid_argument>([&]{ (void)cloth::build_neohookean_constraints(mesh,material(),-1,1); },"Negative energy compliance rejected");
    rejects<std::invalid_argument>([&]{ (void)cloth::build_neohookean_constraints(mesh,material(),0,0); },"Zero energy scale rejected");
    auto bad_material=material(); bad_material.thickness=0;
    rejects<std::invalid_argument>([&]{ (void)cloth::build_neohookean_constraints(mesh,bad_material); },"Invalid energy material rejected");
    auto bad_reference=mesh; bad_reference.reference_positions[1]=bad_reference.reference_positions[0];
    rejects<std::invalid_argument>([&]{ (void)cloth::build_neohookean_constraints(bad_reference,material()); },"Bad energy reference triangle rejected");
    auto bad_indices=mesh; bad_indices.triangles[0][0]=mesh.positions.size();
    rejects<std::invalid_argument>([&]{ (void)cloth::build_neohookean_constraints(bad_indices,material()); },"Bad energy node index rejected");
    cloth::XpbdOptions options; options.energy_tolerance=std::numeric_limits<double>::quiet_NaN();
    rejects<std::invalid_argument>([&]{ cloth::XpbdSolver invalid_options(masses,tiny_scale,{0,0,0},options); },"Non-finite energy tolerance rejected");
}

} // namespace

int main() {
    try {
        value_gradient_and_builder();
        scale_compensated_projection();
        rest_fixed_reset_and_modes();
        domain_overflow_and_invalid_parameters();
        std::cout<<"Energy XPBD values, gradients, scaling, modes, domains and numeric failures passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
