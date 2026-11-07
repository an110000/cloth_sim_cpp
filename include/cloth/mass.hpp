#pragma once

#include "cloth/mesh.hpp"
#include "cloth/sparse.hpp"

namespace cloth {

enum class MassMode { lumped, consistent };

// Both modes return a (3N) x (3N) matrix, with interleaved vertex coordinates.
// Consistent element mass is rho*t*A/12 * [[2,1,1],[1,2,1],[1,1,2]] tensor I3.
// Row sums equal the lumped diagonal. Each scalar vertex mass appears 3 times.
SparseMatrix assemble_mass_matrix(const Mesh& mesh, double density,
                                  double thickness, MassMode mode);

} // namespace cloth
