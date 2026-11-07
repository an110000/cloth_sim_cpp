#pragma once

#include "cloth/mesh.hpp"

#include <vector>

namespace cloth {

struct SimulationState {
    std::vector<Vec3> positions;
    std::vector<Vec3> velocities;
};

} // namespace cloth
