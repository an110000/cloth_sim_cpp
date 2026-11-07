#include "cloth/constraints.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace cloth {
namespace {

constexpr double geometry_floor = 1.0e-12;
constexpr double pi = 3.141592653589793238462643383279502884;

Vec3 add(const Vec3& a, const Vec3& b) {
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}
Vec3 sub(const Vec3& a, const Vec3& b) {
    return {a[0]-b[0], a[1]-b[1], a[2]-b[2]};
}
Vec3 scaled(const Vec3& a, double s) { return {s*a[0], s*a[1], s*a[2]}; }
double dot(const Vec3& a, const Vec3& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
}
double length(const Vec3& a) { return std::hypot(std::hypot(a[0], a[1]), a[2]); }
bool finite(const Vec3& a) {
    return std::all_of(a.begin(), a.end(), [](double v) { return std::isfinite(v); });
}
void checked_compliance(double value) {
    if (!std::isfinite(value) || value < 0.0) throw std::invalid_argument("Compliance must be finite and nonnegative");
}
Vec3 normalized_derivative(const Vec3& unit, const Vec3& derivative, double norm) {
    return scaled(sub(derivative, scaled(unit, dot(unit, derivative))), 1.0/norm);
}

void validate_mesh(const Mesh& mesh, bool check_initial_face_areas = true) {
    const auto count = mesh.reference_positions.size();
    if (count == 0 || mesh.positions.size() != count || mesh.triangles.empty()) {
        throw std::invalid_argument("Constraint builder requires matching nonempty mesh arrays");
    }
    for (std::size_t i=0; i<count; ++i) {
        if (!finite(mesh.positions[i]) || !std::isfinite(mesh.reference_positions[i][0]) ||
            !std::isfinite(mesh.reference_positions[i][1])) {
            throw std::invalid_argument("Constraint builder coordinates must be finite");
        }
    }
    std::set<Triangle> faces;
    std::vector<bool> used(count, false);
    for (const auto& ids : mesh.triangles) {
        if (ids[0]>=count || ids[1]>=count || ids[2]>=count ||
            ids[0]==ids[1] || ids[1]==ids[2] || ids[2]==ids[0]) {
            throw std::invalid_argument("Invalid constraint-builder triangle indices");
        }
        auto key=ids;
        std::sort(key.begin(), key.end());
        if (!faces.insert(key).second) throw std::invalid_argument("Duplicate constraint-builder triangle");
        for (const auto id : ids) used[id]=true;
        const auto& a=mesh.reference_positions[ids[0]];
        const auto& b=mesh.reference_positions[ids[1]];
        const auto& c=mesh.reference_positions[ids[2]];
        const double det=(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);
        if (!std::isfinite(det) || std::abs(det)<1.0e-14) {
            throw std::invalid_argument("Degenerate reference triangle in constraint builder");
        }
        if (check_initial_face_areas) {
            const auto normal=cross(sub(mesh.positions[ids[1]], mesh.positions[ids[0]]),
                                    sub(mesh.positions[ids[2]], mesh.positions[ids[0]]));
            const double normal_length=length(normal);
            if (!finite(normal) || !std::isfinite(normal_length) || normal_length<geometry_floor) {
                throw std::invalid_argument("Degenerate initial triangle in constraint builder");
            }
        }
    }
    if (std::find(used.begin(), used.end(), false)!=used.end()) {
        throw std::invalid_argument("Constraint builder mesh contains an isolated vertex");
    }
}

ConstraintEvaluation distance_evaluation(const DistanceConstraint& constraint,
                                        const std::vector<Vec3>& x) {
    ConstraintEvaluation result;
    result.count=2;
    result.nodes={constraint.i, constraint.j, 0, 0};
    const auto difference=sub(x[constraint.i], x[constraint.j]);
    const double norm=length(difference);
    if (!finite(difference) || !std::isfinite(norm)) throw std::overflow_error("Non-finite distance geometry");
    if (norm<geometry_floor) {
        result.value=std::numeric_limits<double>::infinity();
        result.valid_geometry=false;
        return result;
    }
    result.value=norm-constraint.rest_length;
    result.gradients[0]=scaled(difference, 1.0/norm);
    result.gradients[1]=scaled(result.gradients[0], -1.0);
    return result;
}

ConstraintEvaluation dihedral_evaluation(const DihedralConstraint& constraint,
                                        const std::vector<Vec3>& x) {
    ConstraintEvaluation result;
    result.count=4;
    result.nodes=constraint.nodes;
    const auto& ids=constraint.nodes;
    const auto e=sub(x[ids[1]], x[ids[0]]);
    const auto r1=sub(x[ids[2]], x[ids[0]]);
    const auto r2=sub(x[ids[3]], x[ids[0]]);
    const auto n1=cross(e,r1), n2=cross(r2,e);
    const double le=length(e), ln1=length(n1), ln2=length(n2);
    if (!finite(e) || !finite(r1) || !finite(r2) || !finite(n1) || !finite(n2) ||
        !std::isfinite(le) || !std::isfinite(ln1) || !std::isfinite(ln2)) {
        throw std::overflow_error("Non-finite dihedral geometry");
    }
    if (le<geometry_floor || ln1<geometry_floor || ln2<geometry_floor) {
        result.value=result.angle=std::numeric_limits<double>::infinity();
        result.valid_geometry=false;
        return result;
    }
    const auto eu=scaled(e,1.0/le), u1=scaled(n1,1.0/ln1), u2=scaled(n2,1.0/ln2);
    const double s=dot(eu,cross(u1,u2));
    const double c=std::clamp(dot(u1,u2),-1.0,1.0);
    result.angle=std::atan2(s,c);
    result.value=wrap_angle(result.angle-constraint.rest_angle);
    // Twelve analytic directional derivatives of the normalized geometry.
    // No finite differences enter production evaluation.
    for (std::size_t vertex=0; vertex<4; ++vertex) {
        for (std::size_t axis=0; axis<3; ++axis) {
            Vec3 basis{};
            basis[axis]=1.0;
            const auto de=scaled(basis, vertex==0 ? -1.0 : (vertex==1 ? 1.0 : 0.0));
            const auto dr1=scaled(basis, vertex==0 ? -1.0 : (vertex==2 ? 1.0 : 0.0));
            const auto dr2=scaled(basis, vertex==0 ? -1.0 : (vertex==3 ? 1.0 : 0.0));
            const auto dn1=add(cross(de,r1),cross(e,dr1));
            const auto dn2=add(cross(dr2,e),cross(r2,de));
            const auto deu=normalized_derivative(eu,de,le);
            const auto du1=normalized_derivative(u1,dn1,ln1);
            const auto du2=normalized_derivative(u2,dn2,ln2);
            const double ds=dot(deu,cross(u1,u2))+dot(eu,add(cross(du1,u2),cross(u1,du2)));
            const double dc=dot(du1,u2)+dot(u1,du2);
            result.gradients[vertex][axis]=c*ds-s*dc;
        }
        if (!finite(result.gradients[vertex])) throw std::overflow_error("Non-finite dihedral gradient");
    }
    return result;
}

void validate_energy_parameters(const NeoHookeanEnergyConstraint& value) {
    const auto& material=value.material;
    if (!std::isfinite(value.energy_scale) || value.energy_scale<=0.0 ||
        !std::isfinite(material.young) || material.young<=0.0 ||
        !std::isfinite(material.poisson) || material.poisson<=-1.0 || material.poisson>=0.5 ||
        !std::isfinite(material.thickness) || material.thickness<=0.0 ||
        !std::isfinite(material.j_min) || material.j_min<0.0) {
        throw std::invalid_argument("Invalid Neo-Hookean energy scale or material");
    }
    for (const auto& point:value.reference) {
        for (const auto coordinate:point) {
            if (!std::isfinite(coordinate)) throw std::invalid_argument("Energy reference positions must be finite");
        }
    }
    const auto& a=value.reference[0];
    const auto& b=value.reference[1];
    const auto& c=value.reference[2];
    const double determinant=(b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0]);
    if (!std::isfinite(determinant) || std::abs(determinant)<1.0e-14) {
        throw std::invalid_argument("Degenerate energy reference triangle");
    }
}

ConstraintEvaluation energy_evaluation(const NeoHookeanEnergyConstraint& constraint,
                                      const std::vector<Vec3>& positions) {
    ConstraintEvaluation result;
    result.count=3;
    std::array<Vec3,3> local{};
    for (std::size_t vertex=0;vertex<3;++vertex) {
        result.nodes[vertex]=constraint.nodes[vertex];
        local[vertex]=positions[constraint.nodes[vertex]];
    }
    ElementResult element{};
    try {
        element=evaluate_triangle(constraint.reference,local,constraint.material);
    } catch (const std::domain_error&) {
        // Material-domain failures are not zero-gradient, satisfied constraints.
        result.value=std::numeric_limits<double>::infinity();
        result.valid_geometry=false;
        return result;
    }
    result.value=element.energy/constraint.energy_scale;
    if (!std::isfinite(result.value)) throw std::overflow_error("Non-finite scaled energy constraint");
    for (std::size_t vertex=0;vertex<3;++vertex) {
        for (std::size_t axis=0;axis<3;++axis) {
            result.gradients[vertex][axis]=element.gradient[vertex][axis]/constraint.energy_scale;
        }
        if (!finite(result.gradients[vertex])) throw std::overflow_error("Non-finite scaled energy gradient");
    }
    return result;
}

} // namespace

double wrap_angle(double angle) {
    if (!std::isfinite(angle)) throw std::invalid_argument("Angle must be finite");
    double value=std::fmod(angle+pi,2.0*pi);
    if (value<0.0) value+=2.0*pi;
    return value-pi;
}
ConstraintGroup constraint_group(const Constraint& constraint) noexcept {
    return std::visit([](const auto& value) {
        if constexpr (std::is_same_v<std::decay_t<decltype(value)>, DistanceConstraint>) return value.group;
        else if constexpr (std::is_same_v<std::decay_t<decltype(value)>, DihedralConstraint>) return ConstraintGroup::bending;
        else return ConstraintGroup::neohookean_energy;
    }, constraint);
}
const char* group_name(ConstraintGroup group) noexcept {
    switch (group) {
    case ConstraintGroup::structural: return "structural";
    case ConstraintGroup::shear: return "shear";
    case ConstraintGroup::bending: return "bending";
    case ConstraintGroup::neohookean_energy: return "neohookean_energy";
    }
    return "unknown";
}
double constraint_compliance(const Constraint& c) noexcept {
    return std::visit([](const auto& value) { return value.compliance; },c);
}
double constraint_multiplier(const Constraint& c) noexcept {
    return std::visit([](const auto& value) { return value.lambda_value; },c);
}
void reset_constraint(Constraint& c) noexcept {
    std::visit([](auto& value) { value.lambda_value=0.0; },c);
}
void validate_constraint(const Constraint& constraint, std::size_t vertex_count) {
    std::visit([&](const auto& value) {
        checked_compliance(value.compliance);
        if (!std::isfinite(value.lambda_value)) throw std::invalid_argument("Constraint multiplier must be finite");
        if constexpr (std::is_same_v<std::decay_t<decltype(value)>, DistanceConstraint>) {
            if (value.i>=vertex_count || value.j>=vertex_count || value.i==value.j ||
                !std::isfinite(value.rest_length) || value.rest_length<=0.0 ||
                (value.group!=ConstraintGroup::structural && value.group!=ConstraintGroup::shear)) {
                throw std::invalid_argument("Invalid distance constraint");
            }
        } else if constexpr (std::is_same_v<std::decay_t<decltype(value)>, DihedralConstraint>) {
            bool distinct=true;
            for (std::size_t i=0;i<4;++i) for (std::size_t j=i+1;j<4;++j) {
                if (value.nodes[i]==value.nodes[j]) distinct=false;
            }
            if (!distinct || *std::max_element(value.nodes.begin(),value.nodes.end())>=vertex_count ||
                !std::isfinite(value.rest_angle)) {
                throw std::invalid_argument("Invalid dihedral constraint");
            }
        } else {
            if (value.nodes[0]>=vertex_count || value.nodes[1]>=vertex_count || value.nodes[2]>=vertex_count ||
                value.nodes[0]==value.nodes[1] || value.nodes[1]==value.nodes[2] || value.nodes[2]==value.nodes[0]) {
                throw std::invalid_argument("Invalid energy constraint nodes");
            }
            validate_energy_parameters(value);
        }
    },constraint);
}
ConstraintEvaluation evaluate_constraint(const Constraint& constraint, const std::vector<Vec3>& positions) {
    validate_constraint(constraint,positions.size());
    return std::visit([&](const auto& value) {
        if constexpr (std::is_same_v<std::decay_t<decltype(value)>, DistanceConstraint>) {
            if (!finite(positions[value.i]) || !finite(positions[value.j])) throw std::invalid_argument("Constraint positions must be finite");
            return distance_evaluation(value,positions);
        } else if constexpr (std::is_same_v<std::decay_t<decltype(value)>, DihedralConstraint>) {
            for (const auto id : value.nodes) if (!finite(positions[id])) throw std::invalid_argument("Constraint positions must be finite");
            return dihedral_evaluation(value,positions);
        } else {
            for (const auto id : value.nodes) if (!finite(positions[id])) throw std::invalid_argument("Constraint positions must be finite");
            return energy_evaluation(value,positions);
        }
    },constraint);
}
ProjectionResult project_constraint(Constraint& constraint, std::vector<Vec3>& positions,
                                     const std::vector<double>& inverse_masses, double dt) {
    if (inverse_masses.size()!=positions.size() || !std::isfinite(dt) || dt<=0.0 ||
        !std::isfinite(dt*dt) || dt*dt==0.0) throw std::invalid_argument("Invalid projection mass count or time step");
    const auto evaluation=evaluate_constraint(constraint,positions);
    for (std::size_t local=0; local<evaluation.count; ++local) {
        const auto weight=inverse_masses[evaluation.nodes[local]];
        if (!std::isfinite(weight) || weight<0.0) throw std::invalid_argument("Inverse masses must be finite and nonnegative");
    }
    if (!evaluation.valid_geometry) return {0.0,evaluation.value,false,false};
    const double alpha=constraint_compliance(constraint)/(dt*dt);
    double denominator=alpha;
    for (std::size_t local=0; local<evaluation.count; ++local) {
        const auto weight=inverse_masses[evaluation.nodes[local]];
        if (weight>0.0) denominator+=weight*dot(evaluation.gradients[local],evaluation.gradients[local]);
    }
    if (!std::isfinite(alpha) || !std::isfinite(denominator)) throw std::overflow_error("Non-finite XPBD denominator");
    if (denominator<1.0e-14) return {0.0,evaluation.value,true,false};
    const double old_multiplier=constraint_multiplier(constraint);
    const double delta=(-evaluation.value-alpha*old_multiplier)/denominator;
    const double multiplier=old_multiplier+delta;
    if (!std::isfinite(delta) || !std::isfinite(multiplier)) throw std::overflow_error("Non-finite XPBD multiplier");
    std::array<Vec3,4> candidates{};
    for (std::size_t local=0; local<evaluation.count; ++local) {
        const auto node=evaluation.nodes[local];
        candidates[local]=add(positions[node],scaled(evaluation.gradients[local],inverse_masses[node]*delta));
        if (!finite(candidates[local])) throw std::overflow_error("Non-finite XPBD position correction");
    }
    std::visit([&](auto& value) { value.lambda_value=multiplier; },constraint);
    for (std::size_t local=0; local<evaluation.count; ++local) positions[evaluation.nodes[local]]=candidates[local];
    return {delta,evaluation.value,true,true};
}

std::vector<Constraint> build_bending_constraints(const Mesh& mesh, double compliance) {
    checked_compliance(compliance);
    validate_mesh(mesh);
    struct DirectedFace { std::size_t i,j,opposite; };
    struct Edge { std::vector<DirectedFace> faces; };
    std::map<std::pair<std::size_t,std::size_t>,std::size_t> lookup;
    std::vector<Edge> edges;
    for (const auto& triangle : mesh.triangles) {
        for (std::size_t edge=0; edge<3; ++edge) {
            const auto i=triangle[edge], j=triangle[(edge+1)%3], opposite=triangle[(edge+2)%3];
            const auto key=std::make_pair(std::min(i,j),std::max(i,j));
            const auto [iterator,inserted]=lookup.emplace(key,edges.size());
            if (inserted) edges.emplace_back();
            auto& adjacent=edges[iterator->second].faces;
            adjacent.push_back({i,j,opposite});
            if (adjacent.size()>2) throw std::invalid_argument("Nonmanifold shared edge");
        }
    }
    std::vector<Constraint> result;
    for (const auto& edge : edges) {
        if (edge.faces.size()==1) continue;
        const auto& a=edge.faces[0];
        const auto& b=edge.faces[1];
        if (b.i!=a.j || b.j!=a.i) throw std::invalid_argument("Inconsistent triangle winding on shared edge");
        Constraint constraint=DihedralConstraint{{a.i,a.j,a.opposite,b.opposite},0.0,compliance};
        const auto evaluation=evaluate_constraint(constraint,mesh.positions);
        if (!evaluation.valid_geometry) throw std::invalid_argument("Degenerate initial hinge");
        std::get<DihedralConstraint>(constraint).rest_angle=evaluation.angle;
        result.push_back(std::move(constraint));
    }
    return result;
}

std::vector<Constraint> build_neohookean_constraints(const Mesh& mesh, const Material& material,
                                                    double compliance, double energy_scale) {
    checked_compliance(compliance);
    // Energy admissibility uses the relative det(C) floor from material.cpp,
    // not the absolute normal-length threshold of dihedral bending.
    validate_mesh(mesh,false);
    std::vector<Constraint> result;
    result.reserve(mesh.triangles.size());
    for (const auto& triangle:mesh.triangles) {
        NeoHookeanEnergyConstraint energy;
        energy.nodes=triangle;
        energy.material=material;
        energy.compliance=compliance;
        energy.energy_scale=energy_scale;
        for (std::size_t vertex=0;vertex<3;++vertex) energy.reference[vertex]=mesh.reference_positions[triangle[vertex]];
        Constraint constraint=energy;
        validate_constraint(constraint,mesh.positions.size());
        if (!evaluate_constraint(constraint,mesh.positions).valid_geometry) {
            throw std::invalid_argument("Initial energy constraint is outside the material domain");
        }
        result.push_back(std::move(constraint));
    }
    return result;
}

std::vector<Constraint> build_rectangular_constraints(const Mesh& mesh, std::size_t nx, std::size_t ny,
    double structural_compliance, double shear_compliance, double bending_compliance, bool include_bending) {
    checked_compliance(structural_compliance);
    checked_compliance(shear_compliance);
    checked_compliance(bending_compliance);
    const auto maximum=std::numeric_limits<std::size_t>::max();
    if (nx==0 || ny==0 || nx==maximum || ny==maximum ||
        nx+1>maximum/(ny+1) || nx>maximum/ny/2) {
        throw std::invalid_argument("Invalid rectangular constraint dimensions");
    }
    if (mesh.reference_positions.size()!=(nx+1)*(ny+1) || mesh.triangles.size()!=2*nx*ny) {
        throw std::invalid_argument("Rectangular dimensions do not match mesh topology");
    }
    validate_mesh(mesh);
    std::set<Triangle> expected,actual;
    const auto row_width=nx+1;
    for (std::size_t j=0; j<ny; ++j) for (std::size_t i=0; i<nx; ++i) {
        const auto a=j*row_width+i,b=a+1,c=a+row_width,d=c+1;
        Triangle first{a,b,d},second{a,d,c};
        std::sort(first.begin(),first.end()); std::sort(second.begin(),second.end());
        expected.insert(first); expected.insert(second);
    }
    for (auto triangle : mesh.triangles) { std::sort(triangle.begin(),triangle.end()); actual.insert(triangle); }
    if (actual!=expected) throw std::invalid_argument("Triangles do not match rectangular cell connectivity");
    std::vector<Constraint> result;
    auto append_distance=[&](std::size_t i,std::size_t j,double compliance,ConstraintGroup group) {
        const auto& a=mesh.reference_positions[i]; const auto& b=mesh.reference_positions[j];
        const double rest=std::hypot(a[0]-b[0],a[1]-b[1]);
        Constraint c=DistanceConstraint{i,j,rest,compliance,0.0,group};
        validate_constraint(c,mesh.positions.size());
        result.push_back(std::move(c));
    };
    for (std::size_t j=0; j<=ny; ++j) for (std::size_t i=0; i<=nx; ++i) {
        const auto node=j*row_width+i;
        if (i<nx) append_distance(node,node+1,structural_compliance,ConstraintGroup::structural);
        if (j<ny) append_distance(node,node+row_width,structural_compliance,ConstraintGroup::structural);
    }
    for (std::size_t j=0; j<ny; ++j) for (std::size_t i=0; i<nx; ++i) {
        const auto a=j*row_width+i,b=a+1,c=a+row_width,d=c+1;
        append_distance(a,d,shear_compliance,ConstraintGroup::shear);
        append_distance(b,c,shear_compliance,ConstraintGroup::shear);
    }
    // Validate winding/nonmanifold topology even when bending is disabled.
    auto bending=build_bending_constraints(mesh,bending_compliance);
    if (include_bending) result.insert(result.end(),bending.begin(),bending.end());
    return result;
}

} // namespace cloth
