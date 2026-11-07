#pragma once

#include "cloth/mesh.hpp"

#include <array>

namespace cloth {

// SI units: young [Pa], thickness [m], element energy [J].
struct Material {
    double young = 1.0e6;
    double poisson = 0.3;
    double thickness = 0.001;
    double j_min = 0.05; // Zero selects the unregularized logarithm.
};

struct ElementResult {
    double energy;
    std::array<Vec3, 3> gradient; // dE/dx; force is its negative.
    double area_ratio;
};

// Rows and columns follow [x0,y0,z0,x1,y1,z1,x2,y2,z2]. Units: N/m.
using Mat9 = std::array<std::array<double, 9>, 9>;

struct ElementDerivatives {
    ElementResult value;
    Mat9 hessian;
};

// 2D reference triangle embedded in 3D; matches core/energy.py.
// Throws invalid_argument for invalid input and domain_error for collapse.
ElementResult evaluate_triangle(const std::array<Vec2, 3>& reference,
                                const std::array<Vec3, 3>& positions,
                                const Material& material);

// Analytic second derivatives, including the J < j_min continuation.
// The returned Hessian is not projected to be positive definite or symmetrized.
// evaluate_triangle remains the cheaper energy/gradient-only path.
ElementDerivatives evaluate_triangle_with_hessian(
    const std::array<Vec2, 3>& reference,
    const std::array<Vec3, 3>& positions,
    const Material& material);

} // namespace cloth
