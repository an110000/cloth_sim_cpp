#include "cloth/material.hpp"

#include <cmath>
#include <stdexcept>

namespace cloth {
namespace {

using Mat2 = std::array<Vec2, 2>;
using Mat32 = std::array<Vec2, 3>;

// Shared first-order state; no 9x9 work or allocation on the gradient-only path.
struct ElementState {
    Mat2 inverse_dm{};
    Mat2 inverse_c{};
    Mat32 f{};
    Mat32 fa{};
    double mu = 0.0;
    double lambda = 0.0;
    double ell = 0.0;
    double q = 1.0;
    double q_prime = 0.0;
    double volume = 0.0;
    ElementResult value{};
};

ElementState prepare_element(const std::array<Vec2, 3>& reference,
                             const std::array<Vec3, 3>& positions,
                             const Material& material) {
    if (!std::isfinite(material.young) || material.young <= 0.0 ||
        !std::isfinite(material.poisson) || material.poisson <= -1.0 ||
        material.poisson >= 0.5 || !std::isfinite(material.thickness) ||
        material.thickness <= 0.0 || !std::isfinite(material.j_min) ||
        material.j_min < 0.0) {
        throw std::invalid_argument("Invalid material parameters");
    }
    for (const auto& vertex : reference) {
        for (const auto value : vertex) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument("Non-finite reference position");
            }
        }
    }
    for (const auto& vertex : positions) {
        for (const auto value : vertex) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument("Non-finite current position");
            }
        }
    }

    const double a = reference[1][0] - reference[0][0];
    const double b = reference[2][0] - reference[0][0];
    const double c = reference[1][1] - reference[0][1];
    const double d = reference[2][1] - reference[0][1];
    const double det_dm = a * d - b * c;
    if (!std::isfinite(det_dm) || std::abs(det_dm) < 1.0e-14) {
        throw std::domain_error("Degenerate reference triangle");
    }
    ElementState state;
    state.inverse_dm = {{
        {d / det_dm, -b / det_dm}, {-c / det_dm, a / det_dm}}};
    const auto& inverse_dm = state.inverse_dm;
    auto& f = state.f;
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const double edge1 = positions[1][axis] - positions[0][axis];
        const double edge2 = positions[2][axis] - positions[0][axis];
        for (std::size_t col = 0; col < 2; ++col) {
            f[axis][col] = edge1 * inverse_dm[0][col] +
                           edge2 * inverse_dm[1][col];
        }
    }
    double c00 = 0.0;
    double c01 = 0.0;
    double c11 = 0.0;
    for (const auto& row : f) {
        c00 += row[0] * row[0];
        c01 += row[0] * row[1];
        c11 += row[1] * row[1];
    }
    const double det_c = c00 * c11 - c01 * c01;
    if (!std::isfinite(det_c) || det_c <= 1.0e-14) {
        throw std::domain_error("Collapsed or non-finite triangle");
    }
    const double j = std::sqrt(det_c);
    state.inverse_c = {{{c11 / det_c, -c01 / det_c},
                        {-c01 / det_c, c00 / det_c}}};
    double ell = std::log(j);
    double g = 1.0;
    if (material.j_min > 0.0 && j < material.j_min) {
        const double delta = j - material.j_min;
        const double scaled = delta / material.j_min;
        ell = std::log(material.j_min) + scaled - 0.5 * scaled * scaled;
        g = (j / material.j_min) * (1.0 - scaled);
        state.q_prime = (2.0 / material.j_min) * (1.0 - j / material.j_min);
    }
    const double mu = material.young / (2.0 * (1.0 + material.poisson));
    const double lambda = material.young * material.poisson /
                          ((1.0 + material.poisson) * (1.0 - material.poisson));
    const double volume = 0.5 * std::abs(det_dm) * material.thickness;
    state.mu = mu;
    state.lambda = lambda;
    state.ell = ell;
    state.q = g;
    state.volume = volume;
    const double density = 0.5 * mu * (c00 + c11 - 2.0 - 2.0 * ell) +
                           0.5 * lambda * ell * ell;
    state.value = ElementResult{volume * density, {}, j};
    auto& result = state.value;
    const double coefficient = (lambda * ell - mu) * g;
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const double fa0 = (f[axis][0] * c11 - f[axis][1] * c01) / det_c;
        const double fa1 = (f[axis][1] * c00 - f[axis][0] * c01) / det_c;
        state.fa[axis] = {fa0, fa1};
        const double p0 = mu * f[axis][0] + coefficient * fa0;
        const double p1 = mu * f[axis][1] + coefficient * fa1;
        result.gradient[1][axis] = volume *
            (p0 * inverse_dm[0][0] + p1 * inverse_dm[0][1]);
        result.gradient[2][axis] = volume *
            (p0 * inverse_dm[1][0] + p1 * inverse_dm[1][1]);
        result.gradient[0][axis] = -result.gradient[1][axis] -
                                   result.gradient[2][axis];
    }
    if (!std::isfinite(result.energy)) {
        throw std::domain_error("Non-finite element energy");
    }
    for (const auto& vertex : result.gradient) {
        for (const auto value : vertex) {
            if (!std::isfinite(value)) {
                throw std::domain_error("Non-finite element gradient");
            }
        }
    }
    return state;
}

} // namespace

ElementResult evaluate_triangle(const std::array<Vec2, 3>& reference,
                                const std::array<Vec3, 3>& positions,
                                const Material& material) {
    return prepare_element(reference, positions, material).value;
}

ElementDerivatives evaluate_triangle_with_hessian(
    const std::array<Vec2, 3>& reference,
    const std::array<Vec3, 3>& positions,
    const Material& material) {
    const auto state = prepare_element(reference, positions, material);
    ElementDerivatives result{state.value, {}};
    const double a = state.lambda * state.ell - state.mu;
    const double aq = a * state.q;
    const double chain = state.lambda * state.q * state.q +
                         a * state.value.area_ratio * state.q_prime;

    // Each analytic directional derivative becomes one Hessian column.
    for (std::size_t vertex = 0; vertex < 3; ++vertex) {
        for (std::size_t axis = 0; axis < 3; ++axis) {
            Mat32 df{};
            for (std::size_t col = 0; col < 2; ++col) {
                if (vertex == 0)
                    df[axis][col] = -state.inverse_dm[0][col] - state.inverse_dm[1][col];
                else
                    df[axis][col] = state.inverse_dm[vertex - 1][col];
            }

            Mat2 dc{};
            double dlog_j = 0.0;
            for (std::size_t row = 0; row < 3; ++row) {
                for (std::size_t i = 0; i < 2; ++i) {
                    dlog_j += state.fa[row][i] * df[row][i];
                    for (std::size_t k = 0; k < 2; ++k)
                        dc[i][k] += df[row][i] * state.f[row][k] +
                                    state.f[row][i] * df[row][k];
                }
            }
            // d(C^-1) = -C^-1 dC C^-1; fixed 2x2 loops, not a solver.
            Mat2 da{};
            for (std::size_t i = 0; i < 2; ++i)
                for (std::size_t k = 0; k < 2; ++k)
                    for (std::size_t u = 0; u < 2; ++u)
                        for (std::size_t v = 0; v < 2; ++v)
                            da[i][k] -= state.inverse_c[i][u] * dc[u][v] * state.inverse_c[v][k];

            const std::size_t column = 3 * vertex + axis;
            for (std::size_t row = 0; row < 3; ++row) {
                Vec2 dp{};
                for (std::size_t i = 0; i < 2; ++i) {
                    double dfa = 0.0;
                    for (std::size_t k = 0; k < 2; ++k)
                        dfa += df[row][k] * state.inverse_c[k][i] + state.f[row][k] * da[k][i];
                    dp[i] = state.mu * df[row][i] + chain * dlog_j * state.fa[row][i] + aq * dfa;
                }
                const double dg1 = state.volume *
                    (dp[0] * state.inverse_dm[0][0] + dp[1] * state.inverse_dm[0][1]);
                const double dg2 = state.volume *
                    (dp[0] * state.inverse_dm[1][0] + dp[1] * state.inverse_dm[1][1]);
                result.hessian[row][column] = -dg1 - dg2;
                result.hessian[3 + row][column] = dg1;
                result.hessian[6 + row][column] = dg2;
            }
        }
    }
    for (const auto& row : result.hessian)
        for (const auto value : row)
            if (!std::isfinite(value)) throw std::domain_error("Non-finite element Hessian");
    return result;
}

} // namespace cloth
