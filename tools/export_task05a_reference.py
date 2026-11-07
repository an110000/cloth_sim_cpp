"""Freeze original Neo-Hookean energy-XPBD local projections and trajectories.

No original source is edited. All sampled geometries are valid. Python's finite
`converged` flag is not interpreted as compliant residual convergence. Collapsed
geometry policy is intentionally tested separately by the C++ core tests.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np


def values(array):
    return ' '.join(format(float(value), '.17g') for value in np.asarray(array).ravel())


def material_record(material):
    return values([material.young, material.poisson, material.t, material.j_min or 0.0])


def emit_local(lines, name, constraint, positions, inverse, dt, initial_lambda):
    x = np.asarray(positions, dtype=float).copy()
    inverse = np.asarray(inverse, dtype=float)
    constraint.lambda_value = initial_lambda
    gradient = constraint.compute_gradients(x)
    numerical = np.zeros_like(gradient)
    epsilon = 1.e-7 if name == 'log_continuation' else 1.e-6
    for local, node in enumerate(constraint.nodes):
        for axis in range(3):
            plus, minus = x.copy(), x.copy()
            plus[node, axis] += epsilon
            minus[node, axis] -= epsilon
            numerical[local, axis] = (constraint.evaluate(plus)-constraint.evaluate(minus))/(2*epsilon)
    gradient_error = np.linalg.norm(gradient-numerical)/max(1., np.linalg.norm(gradient))
    if gradient_error > 2.e-6:
        raise RuntimeError(name + ': source energy constraint failed independent finite difference')
    np.testing.assert_allclose(np.sum(gradient, axis=0), 0., atol=2.e-12, rtol=0.)
    value = constraint.evaluate(x)
    lines.extend(['local ' + name,
                  'nodes ' + ' '.join(str(int(i)) for i in constraint.nodes),
                  'reference ' + values(constraint.X),
                  'material ' + material_record(constraint.material),
                  'compliance ' + values(constraint.compliance),
                  'energy_scale ' + values(constraint.energy_scale),
                  'dt ' + values(dt), 'initial_lambda ' + values(initial_lambda),
                  'inverse_masses ' + values(inverse), 'positions ' + values(x),
                  'value ' + values(value), 'gradient ' + values(gradient), 'projections 2'])
    for iteration in range(1, 3):
        previous = x.copy()
        delta, before = constraint.solve(x, inverse, dt)
        after, gradient_after = constraint.evaluate(x), constraint.compute_gradients(x)
        if not np.isfinite(x).all() or not np.isfinite(gradient_after).all() or not np.isfinite([delta, before, after]).all():
            raise RuntimeError(name + ': non-finite source projection')
        if constraint.degenerate_skips:
            raise RuntimeError(name + ': reference crossed invalid geometry')
        np.testing.assert_array_equal(x[inverse == 0], previous[inverse == 0])
        lines.extend(['projection ' + str(iteration), 'delta_lambda ' + values(delta),
                      'value_before ' + values(before), 'lambda ' + values(constraint.lambda_value),
                      'positions ' + values(x), 'value_after ' + values(after),
                      'gradient_after ' + values(gradient_after), 'end_projection'])
    constraint.reset()
    assert constraint.lambda_value == 0. and constraint.degenerate_skips == 0
    lines.extend(['reset_lambda 0', 'end_local'])
    return float(gradient_error)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--output', type=Path, default=root/'tests/data/task05a_reference.txt')
    args = parser.parse_args()
    source = args.source.resolve()
    modules = ['core/mesh.py', 'core/constraints.py', 'core/energy.py', 'core/material.py',
               'core/mass.py', 'solvers/xpbd.py', 'solvers/xpbd_neohookean.py']
    for module in modules:
        if not (source/module).is_file():
            parser.error('Missing original source: ' + module)
    sys.path.insert(0, str(source))
    from core.constraints import NeoHookeanEnergyConstraint
    from core.material import NeoHookeanMaterial
    from core.mesh import ClothMesh
    from core.mass import assemble_mass
    from solvers.xpbd_neohookean import NeoHookeanXPBDSolver

    hashes = {name: hashlib.sha256((source/name).read_bytes()).hexdigest() for name in modules}
    lines = ['CLOTH_TASK05A_REFERENCE 1', '# NumPy ' + np.__version__,
             '# Original energy-XPBD only; no contact or friction.',
             '# C=E/energy_scale, R=C+compliance/dt^2*lambda; both dimensionless.',
             '# Python converged reports finite states, not residual convergence.',
             '# Fixed physical masses are encoded as infinity only in the original solver.']
    lines.extend('# sha256 ' + name + ' ' + digest for name, digest in hashes.items())
    lines.append('locals 6')
    reference = np.array([[0., 0.], [.7, .1], [.2, .9]])
    rest = np.column_stack((reference[:, 0], np.zeros(3), reference[:, 1]))
    deformed = rest + np.array([[.02, .01, -.01], [.09, .06, -.02], [-.03, -.02, .05]])
    material = NeoHookeanMaterial(young=2000., poisson=.27, thickness=.001)
    errors = {}
    # Nodes intentionally use a nontrivial permutation to exercise local/global ordering.
    nodes = np.array([2, 0, 1])
    def local(name, X, local_x, mat, compliance, scale, initial_lambda=0.):
        global_x = np.empty_like(local_x)
        global_x[nodes] = local_x
        constraint = NeoHookeanEnergyConstraint(nodes, X, mat, compliance, scale)
        errors[name] = emit_local(lines, name, constraint, global_x, [0., 2., 3.], .02, initial_lambda)
    local('rest', reference, rest, material, 1.e-6, 1.)
    local('nonlinear', reference, deformed, material, 1.e-6, 1.)
    local('scaled_constraint', reference, deformed, material, 64.e-6, .125)
    local('nonzero_multiplier', reference, deformed, material, 2.e-5, .7, -.002)
    local('hard_energy', reference, deformed, material, 0., 1., -.001)
    continuation_X = np.array([[0., 0.], [1., 0.], [0., 1.]])
    continuation_x = np.array([[0., .02, 0.], [1.03, .02, 0.], [.01, .0201, .018]])
    local('log_continuation', continuation_X, continuation_x, material, 1.e-6, 1.)

    lines.append('trajectories 4')
    largest_raw, largest_compliant = 0., 0.
    trajectory_results = {}
    for name in ['rest', 'freefall', 'pinned_deformed', 'pinned_scaled_equivalent']:
        mesh = ClothMesh(2, 1, 1., .7)
        mesh.x = np.column_stack((mesh.X[:, 0], np.zeros(len(mesh.X)), mesh.X[:, 1]))
        material = NeoHookeanMaterial(young=20000., poisson=.3, thickness=.001)
        physical_mass = assemble_mass(mesh.X, mesh.triangles, 800., material.t, mode='lumped')[::3]
        fixed = np.array([3, 4, 5] if name.startswith('pinned') else [], dtype=int)
        gravity = np.array([0., 0., 0.] if name == 'rest' else [.1, -9.81, .05])
        x, v = mesh.x.copy(), np.zeros_like(mesh.x)
        compliance, scale = 1.e-6, 1.
        if name == 'freefall':
            v[:] = [.12, .07, -.03]
            compliance, scale = 2.e-6, .5
        elif name.startswith('pinned'):
            x[:3, 0] = 1.04*x[:3, 0]+.015
            x[:3, 1] = -.04*(1+mesh.X[:3, 0])
            x[:3, 2] -= .015
            v[:3] = [.015, -.02, .01]
            compliance, scale = 1.e-6, .2
            if name == 'pinned_scaled_equivalent':
                compliance, scale = compliance/9., scale*3.
        scene = SimpleNamespace(mesh=mesh, material=material, lumped_mass=np.repeat(physical_mass, 3),
                                n_vertices=len(mesh.X), fixed_nodes=fixed, gravity=gravity)
        iterations, dt, steps = (10 if name in ('rest', 'freefall') else 30), .005, 8
        solver = NeoHookeanXPBDSolver(scene, compliance, scale, iterations)
        initial_x, initial_v = x.copy(), v.copy()
        lines.extend(['trajectory ' + name, 'vertices ' + str(len(mesh.X)),
                      'triangles ' + str(len(mesh.triangles)), 'reference ' + values(mesh.X),
                      'rest_positions ' + values(mesh.x),
                      'topology ' + ' '.join(str(int(i)) for i in mesh.triangles.ravel()),
                      'material ' + material_record(material), 'density ' + values(800.),
                      'physical_masses ' + values(physical_mass),
                      'compliance ' + values(compliance), 'energy_scale ' + values(scale),
                      'iterations ' + str(iterations), 'gravity ' + values(gravity), 'dt ' + values(dt),
                      'fixed ' + str(len(fixed)) + ''.join(' ' + str(i) for i in fixed),
                      'positions ' + values(x), 'velocities ' + values(v), 'steps ' + str(steps)])
        states = []
        for index in range(1, steps+1):
            old_x = x.copy()
            x, v, info = solver.step(x, v, fixed, dt)
            if not info['converged'] or not np.isfinite(x).all() or not np.isfinite(v).all() or info['degenerate_skips']:
                raise RuntimeError(name + ': invalid original trajectory')
            multipliers = np.array([c.lambda_value for c in solver.constraints])
            raw_values = np.array([c.evaluate(x) for c in solver.constraints])
            raw = float(max(abs(raw_values)))
            compliant = float(max(abs(raw_values+compliance/(dt*dt)*multipliers)))
            if not np.isfinite([raw, compliant]).all():
                raise RuntimeError(name + ': non-finite diagnostics')
            np.testing.assert_allclose(raw, info['constraint_errors']['neohookean_energy'], atol=1.e-13, rtol=1.e-12)
            np.testing.assert_array_equal(x[fixed], old_x[fixed])
            np.testing.assert_array_equal(v[fixed], np.zeros((len(fixed), 3)))
            if name == 'rest':
                np.testing.assert_allclose(x, initial_x, atol=2.e-13, rtol=0.)
                np.testing.assert_allclose(v, initial_v, atol=2.e-11, rtol=0.)
            elif name == 'freefall':
                np.testing.assert_allclose(x, initial_x+index*dt*initial_v+
                                          .5*index*(index+1)*dt*dt*gravity, atol=2.e-12, rtol=0.)
                np.testing.assert_allclose(v, initial_v+index*dt*gravity, atol=2.e-10, rtol=0.)
            lines.extend(['step ' + str(index), 'positions ' + values(x), 'velocities ' + values(v),
                          'multipliers ' + values(multipliers), 'raw_error ' + values(raw),
                          'compliant_error ' + values(compliant), 'end_step'])
            largest_raw, largest_compliant = max(largest_raw, raw), max(largest_compliant, compliant)
            states.append((x.copy(), v.copy(), multipliers.copy()))
        trajectory_results[name] = states
        lines.append('end_trajectory')
    scale_position_error = 0.
    for base, scaled in zip(trajectory_results['pinned_deformed'], trajectory_results['pinned_scaled_equivalent']):
        np.testing.assert_allclose(base[0], scaled[0], atol=2.e-12, rtol=2.e-12)
        np.testing.assert_allclose(base[1], scaled[1], atol=2.e-10, rtol=2.e-10)
        np.testing.assert_allclose(3.*base[2], scaled[2], atol=2.e-12, rtol=2.e-10)
        scale_position_error = max(scale_position_error, float(np.max(abs(base[0]-scaled[0]))))
    lines.append('end')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text('\n'.join(lines)+'\n', encoding='utf-8')
    manifest = {'task': '05A', 'fixture_version': 1, 'numpy': np.__version__,
                'source_sha256': hashes, 'fixture_sha256': hashlib.sha256(args.output.read_bytes()).hexdigest(),
                'counts': {'local_cases': 6, 'local_projections': 12, 'trajectories': 4, 'time_steps': 32},
                'gradient_relative_errors': errors, 'largest_raw_error': largest_raw,
                'largest_compliant_error': largest_compliant, 'scale_equivalence_max_position_difference': scale_position_error,
                'note': 'All source geometries valid. Finite Python converged flag is not residual convergence. No contact/friction.'}
    args.output.with_suffix('.sources.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({k: manifest[k] for k in ('counts', 'gradient_relative_errors', 'largest_raw_error',
                                              'largest_compliant_error', 'scale_equivalence_max_position_difference')}, indent=2))
    print(args.output.resolve())


if __name__ == '__main__':
    main()
