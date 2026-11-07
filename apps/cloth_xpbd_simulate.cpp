#include "cloth/xpbd.hpp"
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
#include <utility>

namespace {
struct Arguments {
    std::size_t steps=8;
    double dt=.01;
    cloth::XpbdOptions options;
    bool bending=true, help=false;
    bool energy_model=false;
    double energy_scale=1.0, energy_compliance=1e-8;
    cloth::GroundContactOptions ground;
    std::filesystem::path output="output/task04";
};
double positive(const std::string& text, bool allow_zero=false) {
    if (text.empty() || text.find_first_not_of("0123456789.eE+-")!=std::string::npos)
        throw std::invalid_argument("Expected a finite positive number");
    std::size_t used=0; const auto value=std::stod(text,&used);
    if (used!=text.size() || !std::isfinite(value) || value<0 || (!allow_zero && value==0))
        throw std::invalid_argument("Number is outside its allowed finite range");
    return value;
}
std::size_t count(const std::string& text) {
    std::size_t value=0;
    const auto parsed=std::from_chars(text.data(),text.data()+text.size(),value);
    if (text.empty() || parsed.ec!=std::errc{} || parsed.ptr!=text.data()+text.size() || value==0)
        throw std::invalid_argument("Expected a positive integer");
    return value;
}
Arguments parse(int argc,char** argv) {
    Arguments args; std::set<std::string> seen;
    for(int i=1;i<argc;++i) {
        const std::string option=argv[i];
        if (!seen.insert(option).second) throw std::invalid_argument("Duplicate option");
        if(option=="--help") {args.help=true;continue;}
        if(option=="--strict") {args.options.require_convergence=true;continue;}
        if(option=="--no-bending") {args.bending=false;continue;}
        if(option!="--steps" && option!="--dt" && option!="--iterations" && option!="--output" &&
           option!="--distance-tolerance" && option!="--angle-tolerance" && option!="--model" &&
           option!="--energy-scale" && option!="--energy-compliance" && option!="--energy-tolerance" &&
           !cloth_cli::contact_option(option,true))
            throw std::invalid_argument("Unknown option: "+option);
        if(++i>=argc) throw std::invalid_argument("Missing option value");
        const std::string value=argv[i];
        if(cloth_cli::contact_option(option,true)) cloth_cli::set_contact(args.ground,option,value);
        else if(option=="--steps") args.steps=count(value);
        else if(option=="--dt") args.dt=positive(value);
        else if(option=="--iterations") args.options.iterations=count(value);
        else if(option=="--distance-tolerance") args.options.distance_tolerance=positive(value);
        else if(option=="--angle-tolerance") args.options.angle_tolerance=positive(value);
        else if(option=="--energy-tolerance") args.options.energy_tolerance=positive(value);
        else if(option=="--energy-scale") args.energy_scale=positive(value);
        else if(option=="--energy-compliance") args.energy_compliance=positive(value,true);
        else if(option=="--model") {
            if(value!="distance" && value!="energy") throw std::invalid_argument("Model must be distance or energy");
            args.energy_model=value=="energy";
        }
        else { if(value.empty() || value.rfind("--",0)==0) throw std::invalid_argument("Invalid output directory"); args.output=value; }
    }
    if (!std::isfinite(args.dt*args.dt) || args.dt*args.dt==0 ||
        static_cast<long double>(args.steps)*args.dt>std::numeric_limits<double>::max())
        throw std::invalid_argument("Time range is outside finite double precision");
    if(!args.energy_model && (seen.count("--energy-scale") || seen.count("--energy-compliance") ||
                             seen.count("--energy-tolerance")))
        throw std::invalid_argument("Energy parameters require --model energy");
    if(args.energy_model && seen.count("--no-bending"))
        throw std::invalid_argument("--no-bending belongs to --model distance; energy model has no hinges");
    cloth::validate_ground_contact_options(args.ground);
    for(const auto& option:seen)
        if(cloth_cli::contact_option(option,true) && !args.ground.enabled)
            throw std::invalid_argument("Contact parameters require --ground-y");
    return args;
}
std::ofstream file(const std::filesystem::path& path) {
    std::ofstream stream(path); stream.exceptions(std::ios::badbit|std::ios::failbit);
    stream<<std::setprecision(17)<<std::boolalpha; return stream;
}
template<class T> void points(std::ostream& out,const T& values) {
    out<<'[';
    for(std::size_t i=0;i<values.size();++i) {
        if(i)out<<',';
        out<<'[';
        for(std::size_t k=0;k<values[i].size();++k) {if(k)out<<',';out<<values[i][k];}
        out<<']';
    }
    out<<']';
}
void number(std::ostream& out,double value) { if(std::isfinite(value)) out<<value; else out<<"null"; }
void frame(std::ostream& out,std::size_t step,double dt,const cloth::SimulationState& state) {
    for(std::size_t i=0;i<state.positions.size();++i) {
        out<<step<<','<<static_cast<double>(step)*dt<<','<<i;
        for(auto v:state.positions[i])out<<','<<v;
        for(auto v:state.velocities[i])out<<','<<v;
        out<<'\n';
    }
}
void config(const Arguments& args,const cloth::Mesh& mesh,const cloth::SimulationState& initial,
            const cloth::Material& material) {
    auto out=file(args.output/"config.json");
    out<<"{\n\"format_version\":2,\"solver\":\""<<(args.energy_model?"xpbd_neohookean_energy":"xpbd_distance")
       <<"\",\"scene\":\"pinned_cloth\",\n"
         "\"units\":{\"position\":\"m\",\"velocity\":\"m/s\",\"time\":\"s\",\"mass\":\"kg\",\"angle\":\"rad\",\n"
         "\"distance_compliance\":\"s^2/kg\",\"bending_compliance\":\"s^2/(kg*m^2)\","
         "\"energy_scale\":\"J\",\"energy_constraint\":\"1\",\"energy_compliance\":\"1/J\"},\n"
         "\"nx_cells\":2,\"ny_cells\":1,\"width\":1,\"height\":0.7,\"density\":800,\"thickness\":0.001,\n"
         "\"structural_compliance\":1e-8,\"shear_compliance\":1e-8,\"bending_compliance\":1e-5,\n"
         "\"gravity\":[0,-9.81,0],\"fixed_vertices\":[3,4,5],";
    cloth_cli::config(out,args.ground);
    out<<",\n\"reference_positions\":"; points(out,mesh.reference_positions);
    out<<",\n\"triangles\":";points(out,mesh.triangles);
    out<<",\n\"rest_positions\":";points(out,mesh.positions);
    out<<",\n\"initial_positions\":";points(out,initial.positions);
    out<<",\n\"initial_velocities\":";points(out,initial.velocities);
    out<<",\n\"dt\":"<<args.dt<<",\"steps\":"<<args.steps<<",\"iterations\":"<<args.options.iterations
       <<",\"bending_enabled\":"<<(!args.energy_model && args.bending)
       <<",\"require_convergence\":"<<args.options.require_convergence
       <<",\"distance_tolerance_m\":"<<args.options.distance_tolerance
       <<",\"angle_tolerance_rad\":"<<args.options.angle_tolerance
       <<",\"energy_tolerance_scaled\":"<<args.options.energy_tolerance
       <<",\"energy_scale_J\":"<<args.energy_scale<<",\"energy_compliance_inverse_J\":"<<args.energy_compliance
       <<",\"material\":{\"young_Pa\":"<<material.young<<",\"poisson\":"<<material.poisson
       <<",\"thickness_m\":"<<material.thickness<<",\"j_min\":"<<material.j_min<<"}"
       <<",\"diagnostic_membrane_energy\":\"sum of NH element energies; not the selected XPBD model potential\"\n}\n";
    out.close();
}
void summary(const Arguments& args,std::size_t completed,const cloth::XpbdStepInfo& info,
             bool diagnostic_energies_valid) {
    auto out=file(args.output/"summary.json");
    out<<"{\n\"format_version\":2,\"model\":\""<<(args.energy_model?"energy":"distance")
       <<"\",\"steps_requested\":"<<args.steps<<",\"steps_completed\":"<<completed
       <<",\"finite_state\":"<<info.finite_state<<",\"committed\":"<<info.committed
       <<",\"projection_converged\":"<<info.projection_converged<<",\"failure_reason\":\""
       <<cloth::failure_name(info.failure_reason)<<"\",\n\"require_convergence\":"<<args.options.require_convergence
       <<",\"last_attempt_sweeps\":"<<info.iterations<<",\"groups\":{";
    for(std::size_t i=0;i<info.groups.size();++i) {
        if(i)out<<',';
        out<<'"'<<cloth::group_name(static_cast<cloth::ConstraintGroup>(i))<<"\":{\"count\":"<<info.groups[i].count
           <<",\"unit\":\""<<(i==3?"1":(i==2?"rad":"m"))<<"\",\"max_violation\":";
        number(out,info.groups[i].max_violation);out<<",\"max_compliant_residual\":";
        number(out,info.groups[i].max_compliant_residual);out<<'}';
    }
    out<<"},\"pre_contact_projection_converged\":"<<info.pre_contact_projection_converged
       <<",\"pre_contact_residuals\":[";
    for(std::size_t i=0;i<info.pre_contact_groups.size();++i) {
        if(i) out<<',';
        number(out,info.pre_contact_groups[i].max_compliant_residual);
    }
    out<<"],\"diagnostic_energies_valid\":"<<diagnostic_energies_valid
       <<",\"last_attempt_diagnostics\":true,";
    cloth_cli::summary(out,info.contact);
    out<<"\n}\n";out.close();
}
std::pair<double,double> energies(const Arguments& args,const cloth::Mesh& mesh,
                                const cloth::SimulationState& state,const cloth::Material& material,
                                bool& diagnostic_energies_valid) {
    long double membrane=0, effective=0;
    try {
    for(const auto& triangle:mesh.triangles) {
        std::array<cloth::Vec2,3> reference;
        std::array<cloth::Vec3,3> positions;
        for(std::size_t i=0;i<3;++i) {
            reference[i]=mesh.reference_positions[triangle[i]];
            positions[i]=state.positions[triangle[i]];
        }
        const auto energy=cloth::evaluate_triangle(reference,positions,material).energy;
        membrane+=energy;
        if(args.energy_model && args.energy_compliance>0) {
            const long double c=static_cast<long double>(energy)/args.energy_scale;
            effective+=c*c/(2*args.energy_compliance);
        }
    }
    } catch(const std::domain_error&) {
        diagnostic_energies_valid=false;
        return {std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN()};
    } catch(const std::overflow_error&) {
        diagnostic_energies_valid=false;
        return {std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN()};
    }
    const auto finite_double=[](long double value) {
        return std::isfinite(value) && std::abs(value)<=std::numeric_limits<double>::max()
            ? static_cast<double>(value) : std::numeric_limits<double>::infinity();
    };
    const auto membrane_value=finite_double(membrane);
    const auto effective_value=finite_double(effective);
    if(!std::isfinite(membrane_value) || (args.energy_model && args.energy_compliance>0 && !std::isfinite(effective_value)))
        diagnostic_energies_valid=false;
    return {membrane_value,args.energy_model && args.energy_compliance>0
        ?effective_value:std::numeric_limits<double>::quiet_NaN()};
}
}
int main(int argc,char** argv) {
    try {
        const auto args=parse(argc,argv);
        if(args.help) {std::cout<<"cloth_xpbd_simulate [--steps N] [--dt H] [--iterations N] [--no-bending]\n"
            " [--strict] [--distance-tolerance M] [--angle-tolerance RAD] [--output DIR]\n"
            " [--model distance|energy] [--energy-scale J] [--energy-compliance INV_J]\n"
            " [--energy-tolerance SCALED]\n"
            " [--ground-y M] [--friction MU] [--static-friction MU_S]\n"
            " [--restitution E] [--contact-beta BETA]\n"
            "Default: finite fixed-budget approximation; strict mode requires compliant residuals.\n";return 0;}
        const auto mesh=cloth::make_rectangular_mesh(2,1,1,.7);
        const cloth::Material material{20000,.3,.001,.05};
        const auto constraints=args.energy_model
            ?cloth::build_neohookean_constraints(mesh,material,args.energy_compliance,args.energy_scale)
            :cloth::build_rectangular_constraints(mesh,2,1,1e-8,1e-8,1e-5,args.bending);
        const cloth::XpbdSolver solver(cloth::lumped_masses(mesh,800,.001),constraints,{0,-9.81,0},args.options,args.ground);
        cloth::SimulationState state{mesh.positions,std::vector<cloth::Vec3>(6,{0,0,0})};
        for(std::size_t i=0;i<3;++i) {
            state.positions[i][0]=1.06*state.positions[i][0]+.025;
            state.positions[i][1]=-.025*(1+mesh.reference_positions[i][0]);
            state.positions[i][2]-=.025; state.velocities[i]={.015,-.02,.01};
        }
        std::filesystem::create_directories(args.output);config(args,mesh,state,material);
        auto trajectory=file(args.output/"trajectory.csv");trajectory<<"step,time,vertex,x,y,z,vx,vy,vz\n";
        frame(trajectory,0,args.dt,state);
        auto contact_csv=file(args.output/"contact.csv");cloth_cli::header(contact_csv);
        auto split_csv=file(args.output/"contact_split.csv");
        split_csv<<"step,group,pre_contact_R,post_contact_R\n";
        auto diagnostics=file(args.output/"diagnostics.csv");
        auto multipliers=file(args.output/"multipliers.csv");multipliers<<"step,constraint,lambda\n";
        diagnostics<<"step,projection_converged,structural_C_m,shear_C_m,bending_C_rad,energy_C_scaled,"
            "structural_R_m,shear_R_m,bending_R_rad,energy_R_scaled,membrane_energy_J,energy_constraint_potential_J\n";
        std::size_t completed=0;cloth::XpbdStepInfo info;
        bool diagnostic_energies_valid=true;
        while(completed<args.steps) {
            auto result=solver.step(state,{3,4,5},args.dt);info=std::move(result.info);
            if(!info.committed)break;
            state=std::move(result.state);++completed;frame(trajectory,completed,args.dt,state);
            diagnostics<<completed<<','<<info.projection_converged;
            for(const auto& group:info.groups)diagnostics<<','<<group.max_violation;
            for(const auto& group:info.groups)diagnostics<<','<<group.max_compliant_residual;
            const auto energy=energies(args,mesh,state,material,diagnostic_energies_valid);
            diagnostics<<',';number(diagnostics,energy.first);
            diagnostics<<',';number(diagnostics,energy.second);
            diagnostics<<'\n';
            for(std::size_t i=0;i<info.multipliers.size();++i)
                multipliers<<completed<<','<<i<<','<<info.multipliers[i]<<'\n';
            cloth_cli::row(contact_csv,completed,info.contact);
            for(std::size_t i=0;i<info.groups.size();++i)
                split_csv<<completed<<','<<cloth::group_name(static_cast<cloth::ConstraintGroup>(i))<<','
                         <<info.pre_contact_groups[i].max_compliant_residual<<','
                         <<info.groups[i].max_compliant_residual<<'\n';
        }
        trajectory.close();diagnostics.close();multipliers.close();contact_csv.close();split_csv.close();
        summary(args,completed,info,diagnostic_energies_valid);
        std::cout<<"steps_completed="<<completed<<'/'<<args.steps<<" committed="<<std::boolalpha<<info.committed
                 <<" projection_converged="<<info.projection_converged<<" failure_reason="<<cloth::failure_name(info.failure_reason)<<'\n';
        return completed==args.steps?0:2;
    } catch(const std::exception& e) {std::cerr<<"cloth_xpbd_simulate: "<<e.what()<<'\n';return 1;}
}
