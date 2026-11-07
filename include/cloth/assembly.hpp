#pragma once

#include "cloth/material.hpp"
#include "cloth/sparse.hpp"

#include <vector>

namespace cloth {

struct ElasticAssembly {
    double energy;
    std::vector<double> gradient; // Vertex-interleaved: x0,y0,z0,x1,y1,z1,...
    SparseMatrix hessian;
};

// Assemble membrane elasticity only. No mass, external force, or boundary
// elimination is hidden here. The mesh and its positions are never modified.
// All vertices must be referenced by valid, nondegenerate triangles.
ElasticAssembly assemble_elastic(const Mesh& mesh, const Material& material);

} // namespace cloth
