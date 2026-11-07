"""Freeze original distance/shear/dihedral XPBD values and ordered trajectories.

This script never edits the original Python modules or earlier fixtures. Python's
XPBDSolver info['converged'] means finite state, not a residual convergence claim.
Fixed masses follow run_experiments.py exactly: physical mass is replaced by inf.
"""
import argparse
import hashlib
from pathlib import Path
import sys

import numpy as np


def values(array):
    return ' '.join(format(float(x), '.17g') for x in np.asarray(array).ravel())


def constraint_record(constraint):
    if hasattr(constraint, 'rest_length'):
        group = 'shear' if constraint.__class__.__name__ == 'ShearConstraint' else 'structural'
        return ('distance ' + group + ' ' + str(constraint.i) + ' ' + str(constraint.j) +
                ' ' + values([constraint.rest_length, constraint.compliance]))
    return ('bending bending ' + ' '.join(str(i) for i in constraint.nodes) + ' ' +
            values([constraint.rest_angle, constraint.compliance]))


def local_value_gradient(constraint, x):
    if hasattr(constraint, 'rest_length'):
        value, direction = constraint.evaluate(x)
        return value, np.array([direction, -direction])
    return constraint.evaluate(x), constraint.compute_gradients(x)


def check_gradient(constraint, x, wrap_angle):
    _, analytic = local_value_gradient(constraint, x)
    numerical = np.zeros_like(x)
    eps = 1.e-6
    for node in range(len(x)):
        for axis in range(3):
            plus, minus = x.copy(), x.copy()
            plus[node, axis] += eps
            minus[node, axis] -= eps
            vp = local_value_gradient(constraint, plus)[0]
            vm = local_value_gradient(constraint, minus)[0]
            difference = vp - vm
            if not hasattr(constraint, 'rest_length'):
                difference = wrap_angle(difference)
            numerical[node, axis] = difference / (2 * eps)
    np.testing.assert_allclose(analytic, numerical, atol=2.e-7, rtol=2.e-7)
    np.testing.assert_allclose(np.sum(analytic, axis=0), 0., atol=2.e-13)


def emit_local(lines, kind, name, constraint, positions, inverse, dt, initial_lambda, wrap_angle):
    x = np.array(positions, dtype=float)
    inverse = np.asarray(inverse, dtype=float)
    check_gradient(constraint, x, wrap_angle)
    constraint.lambda_value = initial_lambda
    value, gradient = local_value_gradient(constraint, x)
    lines.extend(['local ' + name, 'type ' + kind, 'vertices ' + str(len(x)),
                  'constraint ' + constraint_record(constraint),
                  'dt ' + values(dt), 'initial_lambda ' + values(initial_lambda),
                  'inverse_masses ' + values(inverse), 'positions ' + values(x),
                  'value ' + values(value), 'gradient ' + values(gradient), 'projections 2'])
    for index in range(1, 3):
        old = x.copy()
        delta, before = constraint.solve(x, inverse, dt)
        after, after_gradient = local_value_gradient(constraint, x)
        if not np.isfinite(x).all() or not np.isfinite(after) or not np.isfinite(delta):
            raise RuntimeError('Non-finite local reference ' + name)
        np.testing.assert_array_equal(x[inverse == 0], old[inverse == 0])
        lines.extend(['projection ' + str(index), 'delta_lambda ' + values(delta),
                      'value_before ' + values(before), 'lambda ' + values(constraint.lambda_value),
                      'positions ' + values(x), 'value_after ' + values(after),
                      'gradient_after ' + values(after_gradient), 'end_projection'])
    constraint.reset()
    assert constraint.lambda_value == 0.
    lines.extend(['reset_lambda 0', 'end_local'])


def emit_mesh(lines, mesh):
    lines.extend(['grid ' + str(mesh.nx) + ' ' + str(mesh.ny),
                  'vertices ' + str(len(mesh.X)), 'triangles ' + str(len(mesh.triangles)),
                  'reference ' + values(mesh.X), 'rest_positions ' + values(mesh.x),
                  'topology ' + ' '.join(str(int(i)) for i in mesh.triangles.ravel())])


def emit_constraints(lines, constraints):
    lines.append('constraints ' + str(len(constraints)))
    lines.extend('constraint ' + constraint_record(c) for c in constraints)


def group_diagnostics(groups, x, dt):
    raw, compliant = [], []
    for group in groups:
        rows = [local_value_gradient(c, x)[0] for c in group]
        corrected = [value + c.compliance / (dt * dt) * c.lambda_value
                     for c, value in zip(group, rows)]
        raw.append(max((abs(v) for v in rows), default=0.))
        compliant.append(max((abs(v) for v in corrected), default=0.))
    return np.array(raw), np.array(compliant)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--output', type=Path, default=root / 'tests/data/task04_reference.txt')
    args = parser.parse_args()
    source = args.source.resolve()
    modules = ['core/mesh.py', 'core/constraints.py', 'core/bending.py', 'core/mass.py',
               'core/energy.py', 'core/material.py', 'solvers/xpbd.py', 'run_experiments.py']
    for module in modules:
        if not (source / module).is_file():
            parser.error('Missing original source: ' + module)
    sys.path.insert(0, str(source))
    from core.constraints import (DistanceConstraint, build_structural_constraints,
                                  build_shear_constraints)
    from core.bending import DihedralBendingConstraint, build_bending_constraints, wrap_angle
    from core.mass import assemble_mass
    from core.mesh import ClothMesh
    from solvers.xpbd import XPBDSolver

    lines = ['CLOTH_TASK04_REFERENCE 1', '# NumPy ' + np.__version__,
             '# Original code only; finite Python converged flag is not residual convergence.',
             '# All local reference geometries are valid; degenerate policy is tested separately.',
             '# Constraint groups in diagnostics are structural, shear, bending in that order.',
             '# Fixed inverse masses are zero, as in run_experiments.py with mass[fixed]=inf.']
    for module in modules:
        lines.append('# sha256 ' + module + ' ' + hashlib.sha256((source / module).read_bytes()).hexdigest())
    lines.append('locals 7')
    distance_cases = [
        ('distance_soft', [[0., .1, .2], [1.2, -.15, .4]], [2., 3.], 1., 1.e-4, .02, 0.),
        ('distance_fixed_multiplier', [[0., .1, .2], [1.2, -.15, .4]], [0., 2.], 1., 2.e-5, .03, .004),
        ('distance_hard_multiplier', [[-.1, 0., .1], [1.05, .2, -.1]], [2., 1.], 1., 0., .015, -.002)]
    for name, x, inverse, length, compliance, dt, initial_lambda in distance_cases:
        emit_local(lines, 'distance', name, DistanceConstraint(0, 1, length, compliance),
                   x, inverse, dt, initial_lambda, wrap_angle)
    bending_cases = [('bending_flat', 0., 0., 0.), ('bending_ordinary', .7, 0., 0.),
                     ('bending_nonzero_rest', -.4, .2, .0002),
                     ('bending_wrapped_rest', -np.pi + .03, np.pi - .04, -.0001)]
    for name, theta, rest, initial_lambda in bending_cases:
        x = [[0., 0., 0.], [1., 0., 0.], [.2, 1., 0.],
             [.3, -np.cos(theta), -np.sin(theta)]]
        constraint = DihedralBendingConstraint((0, 1), (2, 3), rest, 1.e-5)
        emit_local(lines, 'bending', name, constraint, x, [0., 2., 3., 4.], .02,
                   initial_lambda, wrap_angle)

    def build(mesh, bending=True):
        structural = build_structural_constraints(mesh, compliance=1.e-8)
        shear = build_shear_constraints(mesh, compliance=1.e-8)
        bends = build_bending_constraints(mesh, compliance=1.e-5) if bending else []
        return structural, shear, bends

    lines.append('builders 2')
    for nx, ny in [(2, 1), (2, 2)]:
        mesh = ClothMesh(nx, ny, 1., .7)
        mesh.x = np.column_stack((mesh.X[:, 0], .04 * mesh.X[:, 0] * mesh.X[:, 1], mesh.X[:, 1]))
        groups = build(mesh)
        constraints = [c for group in groups for c in group]
        lines.append('builder grid_' + str(nx) + '_' + str(ny))
        emit_mesh(lines, mesh)
        lines.append('compliances ' + values([1.e-8, 1.e-8, 1.e-5]))
        emit_constraints(lines, constraints)
        lines.append('end_builder')

    lines.append('trajectories 4')
    largest_compliant = np.zeros(3)
    for name in ['rest', 'freefall', 'pinned_deformed', 'bending_folded']:
        mesh = ClothMesh(2, 1, 1., .7)
        mesh.x = np.column_stack((mesh.X[:, 0], np.zeros(len(mesh.X)), mesh.X[:, 1]))
        groups = build(mesh)
        constraints = [c for group in groups for c in group]
        x = mesh.x.copy()
        v = np.zeros_like(x)
        fixed = [3, 4, 5] if name in ('pinned_deformed', 'bending_folded') else []
        gravity = np.array([0., 0., 0.] if name == 'rest' else [0., -9.81, 0.])
        iterations = 10 if name in ('rest', 'freefall') else 30
        if name == 'freefall':
            v[:] = [.12, .07, -.03]
        elif name == 'pinned_deformed':
            x[:3, 0] = 1.06 * x[:3, 0] + .025
            x[:3, 1] = -.025 * (1 + mesh.X[:3, 0])
            x[:3, 2] -= .025
            v[:3] = [.015, -.02, .01]
        elif name == 'bending_folded':
            x[1, 1] = -.18
            x[2, 1] = -.09
            v[:3] = [.01, -.03, 0.]
        density, thickness, dt, steps = 800., .001, .01, 8
        physical_mass = assemble_mass(mesh.X, mesh.triangles, density, thickness, mode='lumped')[::3]
        mass = physical_mass.copy()
        mass[fixed] = np.inf
        solver = XPBDSolver(mass, constraints, gravity, iterations,
                            dict(zip(['structural', 'shear', 'bending'], groups)))
        initial_x, initial_v = x.copy(), v.copy()
        lines.append('trajectory ' + name)
        emit_mesh(lines, mesh)
        lines.extend(['density_thickness ' + values([density, thickness]),
                      'physical_masses ' + values(physical_mass),
                      'compliances ' + values([1.e-8, 1.e-8, 1.e-5]),
                      'iterations ' + str(iterations), 'gravity ' + values(gravity),
                      'dt ' + values(dt), 'fixed ' + str(len(fixed)) + ''.join(' ' + str(i) for i in fixed),
                      'positions ' + values(x), 'velocities ' + values(v)])
        emit_constraints(lines, constraints)
        lines.append('steps ' + str(steps))
        for index in range(1, steps + 1):
            old_x = x.copy()
            x, v, info = solver.step(x, v, fixed, dt)
            if not info['converged'] or not np.isfinite(x).all() or not np.isfinite(v).all():
                raise RuntimeError(name + ': original solver returned non-finite state')
            raw, compliant = group_diagnostics(groups, x, dt)
            if not np.isfinite(raw).all() or not np.isfinite(compliant).all():
                raise RuntimeError(name + ': invalid final diagnostics')
            np.testing.assert_allclose(raw, [info['constraint_errors'][g]
                                           for g in ['structural', 'shear', 'bending']], atol=1.e-14)
            np.testing.assert_array_equal(x[fixed], old_x[fixed])
            np.testing.assert_array_equal(v[fixed], np.zeros((len(fixed), 3)))
            if name == 'rest':
                np.testing.assert_allclose(x, initial_x, atol=2.e-13, rtol=0.)
                np.testing.assert_allclose(v, initial_v, atol=2.e-11, rtol=0.)
            elif name == 'freefall':
                expected_x = (initial_x + index * dt * initial_v +
                              .5 * index * (index + 1) * dt * dt * gravity)
                np.testing.assert_allclose(x, expected_x, atol=2.e-12, rtol=0.)
                np.testing.assert_allclose(v, initial_v + index * dt * gravity, atol=2.e-10, rtol=0.)
            lines.extend(['step ' + str(index), 'positions ' + values(x), 'velocities ' + values(v),
                          'multipliers ' + values([c.lambda_value for c in constraints]),
                          'raw_errors ' + values(raw), 'compliant_errors ' + values(compliant), 'end_step'])
            largest_compliant = np.maximum(largest_compliant, compliant)
        lines.append('end_trajectory')
    lines.append('end')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print('Exported 7 local cases with 2 projections each, 2 ordered builders, 4 trajectories / 32 steps')
    print('Largest final compliant residuals [structural m, shear m, bending rad]:', largest_compliant)
    print(args.output.resolve())


if __name__ == '__main__':
    main()
