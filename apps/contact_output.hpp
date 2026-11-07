#pragma once

#include "cloth/contact.hpp"

#include <cmath>
#include <ostream>
#include <stdexcept>
#include <string>

namespace cloth_cli {
inline bool contact_option(const std::string& key, bool projected) {
    return key=="--ground-y" || key=="--friction" || key=="--static-friction" ||
        (projected ? key=="--restitution" || key=="--contact-beta" : key=="--contact-eta");
}
inline void set_contact(cloth::GroundContactOptions& ground, const std::string& key,
                        const std::string& text) {
    if(text.empty() || text.find_first_not_of("0123456789.eE+-")!=std::string::npos)
        throw std::invalid_argument("Contact option needs a finite number");
    std::size_t used=0; const double value=std::stod(text,&used);
    if(used!=text.size() || !std::isfinite(value))
        throw std::invalid_argument("Contact option needs a finite number");
    if(key=="--ground-y") {ground.enabled=true;ground.ground_y=value;}
    else if(key=="--friction") ground.friction.kinetic=value;
    else if(key=="--static-friction") ground.friction.static_coefficient=value;
    else if(key=="--contact-eta") ground.penalty_eta=value;
    else if(key=="--restitution") ground.restitution=value;
    else if(key=="--contact-beta") ground.stabilization_beta=value;
    else throw std::invalid_argument("Unknown contact option");
}
inline void config(std::ostream& out,const cloth::GroundContactOptions& g) {
    out<<"\"contact\":"<<g.enabled<<",\"contact_parameters\":{\"ground_y_m\":"<<g.ground_y
       <<",\"normal\":[0,1,0],\"kinetic_friction\":"<<g.friction.kinetic
       <<",\"static_friction\":"<<g.friction.static_coefficient.value_or(g.friction.kinetic)
       <<",\"restitution\":"<<g.restitution<<",\"stabilization_beta\":"<<g.stabilization_beta
       <<",\"penalty_eta\":"<<g.penalty_eta<<",\"nodal_mass\":\"physical_lumped\"}";
}
inline void summary(std::ostream& out,const cloth::ContactInfo& c) {
    out<<"\"contact\":{\"active_count\":"<<c.active_count<<",\"sticking_count\":"<<c.sticking_count
       <<",\"sliding_count\":"<<c.sliding_count<<",\"fixed_penetration_count\":"<<c.fixed_penetration_count
       <<",\"max_penetration_before_m\":"<<c.max_penetration_before
       <<",\"max_penetration_after_m\":"<<c.max_penetration_after
       <<",\"max_fixed_penetration_m\":"<<c.max_fixed_penetration
       <<",\"normal_impulse_Ns\":"<<c.normal_impulse_total
       <<",\"tangential_impulse_Ns\":"<<c.tangential_impulse_total
       <<",\"lumped_friction_dissipation_J\":"<<c.friction_dissipation
       <<",\"actual_friction_kinetic_change_J\":"<<c.friction_kinetic_change
       <<",\"normal_kinetic_change_J\":"<<c.normal_kinetic_change
       <<",\"penalty_energy_J\":"<<c.penalty_energy<<'}';
}
inline void header(std::ostream& out) {
    out<<"step,active,sticking,sliding,fixed_penetration_count,penetration_before_m,penetration_after_m,"
          "normal_impulse_Ns,tangential_impulse_Ns,lumped_friction_dissipation_J,"
          "actual_friction_kinetic_change_J,normal_kinetic_change_J,penalty_energy_J\n";
}
inline void row(std::ostream& out,std::size_t step,const cloth::ContactInfo& c) {
    out<<step<<','<<c.active_count<<','<<c.sticking_count<<','<<c.sliding_count<<','
       <<c.fixed_penetration_count<<','<<c.max_penetration_before<<','<<c.max_penetration_after<<','
       <<c.normal_impulse_total<<','<<c.tangential_impulse_total<<','<<c.friction_dissipation<<','
       <<c.friction_kinetic_change<<','<<c.normal_kinetic_change<<','<<c.penalty_energy<<'\n';
}
} // namespace cloth_cli
