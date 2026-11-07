#include "cloth/constraints.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

constexpr double pi=3.14159265358979323846;
void check(bool value,const std::string& message) {
    if (!value) throw std::runtime_error(message);
}
void near(double actual,double expected,double tolerance,const std::string& message) {
    check(std::isfinite(actual) && std::abs(actual-expected)<=tolerance,message);
}
template<class Exception,class Function> void rejects(Function function,const std::string& message) {
    try { function(); } catch (const Exception&) { return; }
    throw std::runtime_error(message);
}
cloth::Vec3 cross(const cloth::Vec3& a,const cloth::Vec3& b) {
    return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
}

void distance_projection() {
    cloth::Constraint constraint=cloth::DistanceConstraint{0,1,1.0,0.02};
    std::vector<cloth::Vec3> x{{0,0,0},{2,0,0}};
    const auto evaluation=cloth::evaluate_constraint(constraint,x);
    near(evaluation.value,1.0,1e-14,"Distance residual");
    near(evaluation.gradients[0][0],-1.0,1e-14,"Distance direction sign");
    const auto old=x;
    const auto result=cloth::project_constraint(constraint,x,{1.0,0.5},0.1);
    near(result.delta_lambda,-1.0/3.5,1e-14,"Hand-computed compliant multiplier");
    near(x[0][0],1.0/3.5,1e-14,"First weighted correction");
    near(x[1][0],2.0-0.5/3.5,1e-14,"Second weighted correction");
    near((x[0][0]-old[0][0])+2.0*(x[1][0]-old[1][0]),0.0,1e-14,"Pair preserves center of mass");
    const auto updated=cloth::evaluate_constraint(constraint,x);
    near(updated.value+2.0*cloth::constraint_multiplier(constraint),0.0,1e-14,"Soft C need not be zero; compliant equation is");
    const auto repeated=cloth::project_constraint(constraint,x,{1.0,0.5},0.1);
    near(repeated.delta_lambda,0.0,1e-14,"Accumulated multiplier prevents oversolving soft pair");
    cloth::reset_constraint(constraint);
    near(cloth::constraint_multiplier(constraint),0.0,0.0,"Multiplier reset");

    constraint=cloth::DistanceConstraint{0,1,1.0,0.0};
    x=old;
    (void)cloth::project_constraint(constraint,x,{0.0,0.5},0.1);
    check(x[0]==old[0],"Zero inverse mass keeps fixed point unchanged");
    near(x[1][0],1.0,1e-14,"Fixed particle omitted from projection denominator");
    x=old;
    cloth::reset_constraint(constraint);
    const auto blocked=cloth::project_constraint(constraint,x,{0.0,0.0},0.1);
    check(!blocked.applied && blocked.valid_geometry && x==old,"All-fixed hard constraint is not projected");
    near(blocked.value_before,1.0,1e-14,"All-fixed violation remains visible");
}

void dihedral_gradients() {
    const double angle=0.4;
    std::vector<cloth::Vec3> x{{0,0,0},{1,0,0},{0.1,1,0},{0.2,-std::cos(angle),-std::sin(angle)}};
    cloth::Constraint constraint=cloth::DihedralConstraint{{0,1,2,3},0.0,1e-5};
    auto evaluation=cloth::evaluate_constraint(constraint,x);
    near(evaluation.angle,angle,2e-14,"Signed hinge angle convention");
    const double epsilon=1e-6;
    for (std::size_t vertex=0;vertex<4;++vertex) for (std::size_t axis=0;axis<3;++axis) {
        auto plus=x,minus=x;
        plus[vertex][axis]+=epsilon;
        minus[vertex][axis]-=epsilon;
        const double fd=cloth::wrap_angle(cloth::evaluate_constraint(constraint,plus).value-
                                         cloth::evaluate_constraint(constraint,minus).value)/(2.0*epsilon);
        near(evaluation.gradients[vertex][axis],fd,8e-8,"Analytic dihedral gradient versus centered difference");
    }
    cloth::Vec3 sum{},torque{};
    for (std::size_t vertex=0;vertex<4;++vertex) {
        const auto moment=cross(x[vertex],evaluation.gradients[vertex]);
        for (std::size_t axis=0;axis<3;++axis) {
            sum[axis]+=evaluation.gradients[vertex][axis];
            torque[axis]+=moment[axis];
        }
    }
    for (std::size_t axis=0;axis<3;++axis) {
        near(sum[axis],0.0,2e-14,"Hinge translation gradient null mode");
        near(torque[axis],0.0,2e-14,"Hinge rotational gradient null mode");
    }
    auto transformed=x;
    for (auto& point:transformed) point={-point[1]+0.2,point[0]-0.3,point[2]+0.7};
    const auto rotated=cloth::evaluate_constraint(constraint,transformed);
    near(rotated.value,evaluation.value,2e-14,"Hinge rigid-motion invariance");
    for (std::size_t vertex=0;vertex<4;++vertex) {
        const auto& g=evaluation.gradients[vertex];
        const cloth::Vec3 expected{-g[1],g[0],g[2]};
        for (std::size_t axis=0;axis<3;++axis) near(rotated.gradients[vertex][axis],expected[axis],2e-14,"Hinge gradient covariance");
    }
    auto doubled=x;
    for (auto& point:doubled) for (auto& coordinate:point) coordinate*=2.0;
    const auto scaled=cloth::evaluate_constraint(constraint,doubled);
    near(scaled.value,evaluation.value,2e-14,"Angle scale invariance");
    for (std::size_t vertex=0;vertex<4;++vertex) for (std::size_t axis=0;axis<3;++axis) {
        near(scaled.gradients[vertex][axis],0.5*evaluation.gradients[vertex][axis],2e-14,"Angle gradient inverse-length units");
    }
    auto reversed=constraint;
    std::get<cloth::DihedralConstraint>(reversed).nodes={1,0,3,2};
    near(cloth::evaluate_constraint(reversed,x).angle,angle,2e-14,"Reversing edge and faces preserves signed angle");
    near(cloth::wrap_angle(pi),-pi,1e-14,"Positive pi wraps to negative pi");
    near(cloth::wrap_angle(-pi),-pi,1e-14,"Negative pi is included");
    near(cloth::wrap_angle(3.0*pi),-pi,1e-14,"Multiple rotations wrap");
    std::get<cloth::DihedralConstraint>(constraint).rest_angle=angle-3.3;
    near(cloth::evaluate_constraint(constraint,x).value,3.3-2.0*pi,1e-14,"Residual uses shortest wrapped difference");
}

void builder_order_and_rest_state() {
    auto mesh=cloth::make_rectangular_mesh(2,1,1.0,0.7);
    const auto constraints=cloth::build_rectangular_constraints(mesh,2,1);
    check(constraints.size()==14,"2x1 mesh has seven structural, four shear, three hinges");
    const std::array<std::array<std::size_t,2>,11> pairs{{{0,1},{0,3},{1,2},{1,4},{2,5},{3,4},{4,5},
                                                        {0,4},{1,3},{1,5},{2,4}}};
    for (std::size_t i=0;i<pairs.size();++i) {
        const auto& value=std::get<cloth::DistanceConstraint>(constraints[i]);
        check(value.i==pairs[i][0] && value.j==pairs[i][1],"Original structural/shear insertion order");
        check(cloth::constraint_group(constraints[i])==(i<7?cloth::ConstraintGroup::structural:cloth::ConstraintGroup::shear),"Distance group tag");
    }
    const std::array<std::array<std::size_t,4>,3> hinges{{{1,4,0,5},{4,0,1,3},{5,1,2,4}}};
    for (std::size_t i=0;i<3;++i) {
        check(std::get<cloth::DihedralConstraint>(constraints[11+i]).nodes==hinges[i],"Original first-seen interior-edge order");
    }
    for (const auto& constraint:constraints) near(cloth::evaluate_constraint(constraint,mesh.positions).value,0.0,1e-14,"Rest constraints vanish");
    check(cloth::build_rectangular_constraints(mesh,2,1,1e-8,1e-8,1e-5,false).size()==11,"Bending can be disabled independently");
    mesh.positions[4][1]=0.2;
    const auto bent=cloth::build_bending_constraints(mesh);
    bool nonzero_rest=false;
    for (const auto& constraint:bent) {
        nonzero_rest=nonzero_rest || std::abs(std::get<cloth::DihedralConstraint>(constraint).rest_angle)>0.01;
        near(cloth::evaluate_constraint(constraint,mesh.positions).value,0.0,1e-14,"Rest angle uses actual initial embedding");
    }
    check(nonzero_rest,"Curved initial configuration does not get flat rest angles");
}

void invalid_geometry_and_parameters() {
    rejects<std::invalid_argument>([]{ cloth::validate_constraint(cloth::Constraint{},2); },"Default constraint is safely invalid");
    cloth::Constraint distance=cloth::DistanceConstraint{0,1,1.0,0.0};
    std::vector<cloth::Vec3> collapsed{{0,0,0},{0,0,0}};
    const auto before=collapsed;
    const auto bad=cloth::project_constraint(distance,collapsed,{1,1},0.01);
    check(!bad.valid_geometry && !bad.applied && std::isinf(bad.value_before),"Collapsed distance cannot report satisfied zero");
    check(collapsed==before && cloth::constraint_multiplier(distance)==0.0,"Collapsed distance does not mutate input");
    cloth::Constraint hinge=cloth::DihedralConstraint{{0,1,2,3},0.0};
    std::vector<cloth::Vec3> flat_line{{0,0,0},{1,0,0},{2,0,0},{0,-1,0}};
    const auto invalid=cloth::evaluate_constraint(hinge,flat_line);
    check(!invalid.valid_geometry && std::isinf(invalid.value),"Degenerate hinge not treated as zero angle");
    std::vector<cloth::Vec3> x{{0,0,0},{2,0,0}};
    rejects<std::invalid_argument>([&]{ (void)cloth::project_constraint(distance,x,{1,-1},0.01); },"Negative inverse mass rejected");
    rejects<std::invalid_argument>([&]{ (void)cloth::project_constraint(distance,x,{1,1},0.0); },"Zero dt rejected");
    std::get<cloth::DistanceConstraint>(distance).compliance=std::numeric_limits<double>::max();
    const auto original=x;
    rejects<std::overflow_error>([&]{ (void)cloth::project_constraint(distance,x,{1,1},1e-100); },"Overflowed compliance scaling rejected");
    check(x==original && cloth::constraint_multiplier(distance)==0.0,"Overflow projection rolls back all local values");
    std::get<cloth::DistanceConstraint>(distance).compliance=-1.0;
    rejects<std::invalid_argument>([&]{ cloth::validate_constraint(distance,2); },"Negative compliance rejected");
    distance=cloth::DistanceConstraint{0,0,1.0};
    rejects<std::invalid_argument>([&]{ cloth::validate_constraint(distance,2); },"Repeated distance nodes rejected");

    auto mesh=cloth::make_rectangular_mesh(1,1,1,1);
    auto wrong_winding=mesh;
    std::swap(wrong_winding.triangles[1][1],wrong_winding.triangles[1][2]);
    rejects<std::invalid_argument>([&]{ (void)cloth::build_bending_constraints(wrong_winding); },"Inconsistent shared-edge winding rejected");
    auto reference_collapse=mesh;
    reference_collapse.reference_positions[1]=reference_collapse.reference_positions[0];
    rejects<std::invalid_argument>([&]{ (void)cloth::build_bending_constraints(reference_collapse); },"Degenerate reference triangle rejected");
    auto overflow_normal=mesh;
    overflow_normal.positions[0]={0,0,0};
    overflow_normal.positions[1]={1e154,0,-1e154};
    overflow_normal.positions[3]={0,1.7e154,0};
    overflow_normal.positions[2]={-1,0,1};
    rejects<std::invalid_argument>([&]{ (void)cloth::build_bending_constraints(overflow_normal); },"Finite normal components with overflowing norm rejected");
    rejects<std::invalid_argument>([&]{ (void)cloth::build_rectangular_constraints(mesh,2,1); },"Grid shape mismatch rejected");
    cloth::Mesh nonmanifold;
    nonmanifold.reference_positions={{0,0},{1,0},{0,1},{0,-1},{0.5,2}};
    for (const auto& point:nonmanifold.reference_positions) nonmanifold.positions.push_back({point[0],0,point[1]});
    nonmanifold.triangles={{0,1,2},{1,0,3},{0,1,4}};
    rejects<std::invalid_argument>([&]{ (void)cloth::build_bending_constraints(nonmanifold); },"Three faces on edge rejected rather than silently skipped");
}

} // namespace

int main() {
    try {
        distance_projection();
        dihedral_gradients();
        builder_order_and_rest_state();
        invalid_geometry_and_parameters();
        std::cout<<"Distance/dihedral projections, derivatives, invariants, topology and failures passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
