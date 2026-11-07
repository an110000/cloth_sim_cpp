"""Freeze contact/friction from original Python, with explicit legacy audits.

Only valid beta=0, restitution=0 incoming ground cases are equality references.
The consistent penalty trajectory maps eta_new=eta_old/2 and uses zero friction:
the original takes diag(M), while the port uses physical row-sum nodal masses.
This is a parameter-mapped comparison, not same-parameter equivalence.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
import scipy


def values(array):
    return ' '.join(format(float(v), '.17g') for v in np.asarray(array).ravel())


def fixed_record(fixed):
    return 'fixed ' + str(len(fixed)) + ''.join(' ' + str(int(i)) for i in fixed)


def penalty(x, masses, fixed, dt, ground_y, eta):
    """Isolate original implicit_sparse.step's nested quadratic penalty formula."""
    pen = np.maximum(ground_y-np.asarray(x)[:, 1], 0.)
    kappa = eta*np.asarray(masses)/(dt*dt)
    kappa = kappa.copy()
    kappa[np.asarray(fixed, dtype=int)] = 0.
    return .5*np.dot(kappa, pen*pen), -kappa*pen, kappa*(pen > 0.), dt*kappa*pen, pen


def emit_contact_info(lines, prefix, info):
    lines.extend([prefix+'_active ' + str(info['contact_count']),
                  prefix+'_sticking ' + str(info['sticking_count']),
                  prefix+'_sliding ' + str(info['sliding_count']),
                  prefix+'_max_penetration ' + values(info['max_penetration']),
                  prefix+'_normal_total ' + values(info['normal_impulse']),
                  prefix+'_tangent_total ' + values(info['tangential_impulse']),
                  prefix+'_dissipation ' + values(info['friction_dissipation'])])


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--output', type=Path, default=root/'tests/data/task05b_reference.txt')
    args = parser.parse_args()
    source = args.source.resolve()
    modules = ['core/contact.py', 'core/friction.py', 'core/simulation.py', 'core/mesh.py',
               'core/mass.py', 'core/mass_sparse.py', 'core/material.py', 'core/energy.py',
               'core/constraints.py', 'core/bending.py', 'solvers/implicit.py',
               'solvers/implicit_sparse.py', 'solvers/xpbd.py', 'solvers/xpbd_neohookean.py']
    for name in modules:
        if not (source/name).is_file():
            parser.error('Missing original source: '+name)
    sys.path.insert(0, str(source))
    from core.contact import apply_ground_contact, apply_friction_from_normal_impulse
    from core.friction import apply_coulomb_friction
    from core.mesh import ClothMesh
    from core.material import NeoHookeanMaterial
    from core.mass import assemble_mass
    from core.constraints import build_structural_constraints, build_shear_constraints
    from core.bending import build_bending_constraints
    from solvers.implicit_sparse import ImplicitEulerSolverSparse
    from solvers.xpbd import XPBDSolver
    from solvers.xpbd_neohookean import NeoHookeanXPBDSolver

    hashes = {name: hashlib.sha256((source/name).read_bytes()).hexdigest() for name in modules}
    lines = ['CLOTH_TASK05B_REFERENCE 1', '# NumPy '+np.__version__+'; SciPy '+scipy.__version__,
             '# Equality references use beta=0, restitution=0 and fixed nodes above ground.',
             '# Legacy audits document intentional changes; their old outputs are not port targets.',
             '# Original zero-vt sticking counts differ between friction.py and contact.py.',
             '# Consistent trajectory: zero friction, eta_port=eta_original/2 because original uses diag(M).']
    lines.extend('# sha256 '+name+' '+digest for name, digest in hashes.items())
    lines.append('friction_cases 4')
    for name, velocity, impulse in [('sticking', [.02, -.4, -.01], .5),
                                     ('sliding', [1.2, -.4, .7], .5),
                                     ('zero_normal', [1.2, -.4, .7], 0.),
                                     ('zero_tangent', [0., -.4, 0.], .5)]:
        v = np.array([velocity], dtype=float)
        mass, normal = np.array([2.]), np.array([impulse])
        result = apply_coulomb_friction(v, normal, mass, .3, static_friction_coeff=.6)
        after = result['velocities_new']
        exact_dissipation = .5*np.sum(mass[:, None]*(v[:, [0, 2]]**2-after[:, [0, 2]]**2))
        np.testing.assert_allclose(exact_dissipation, result['friction_dissipation_total'], atol=1.e-14, rtol=1.e-13)
        np.testing.assert_array_equal(after[:, 1], v[:, 1])
        assert result['friction_satisfied'].all() and exact_dissipation >= 0.
        lines.extend(['friction '+name, 'vertices 1', 'masses '+values(mass), 'normal_impulses '+values(normal),
                      'friction_options .3 .6', 'velocities '+values(v), 'expected_velocities '+values(after),
                      'tangential_impulses '+values(result['tangential_impulses']),
                      'dissipation_per_node '+values(result['friction_dissipation']),
                      'source_sticking '+str(result['sticking_count']), 'source_sliding '+str(result['sliding_count']),
                      'end_friction'])

    ground_cases = [('incoming_sliding', -.015, [1.2, -.4, .7]),
                    ('incoming_sticking', -.015, [.02, -.4, -.01]),
                    ('plane_incoming', 0., [.4, -.1, .2]),
                    ('separated', .1, [.4, -.3, .2]),
                    ('zero_tangent', -.015, [0., -.5, 0.])]
    lines.append('ground_cases '+str(len(ground_cases)))
    for name, height, velocity in ground_cases:
        x = np.array([[.2, height, .1]])
        v = np.array([velocity], dtype=float)
        initial_x, initial_v = x.copy(), v.copy()
        masses = np.array([2.])
        info = apply_ground_contact(x, v, 0., .3, dt=.02, masses=masses,
                                    static_friction=.6, restitution=0., contact_beta=0.)
        actual_jn = masses*(v[:, 1]-initial_v[:, 1])
        np.testing.assert_allclose(max(0., actual_jn.sum()), info['normal_impulse'], atol=1.e-14)
        dissipation = .5*np.sum(masses[:, None]*(initial_v[:, [0, 2]]**2-v[:, [0, 2]]**2))
        np.testing.assert_allclose(dissipation, info['friction_dissipation'], atol=1.e-14)
        assert np.min(x[:, 1]) >= 0. and dissipation >= 0.
        lines.extend(['ground '+name, 'vertices 1', 'masses '+values(masses), 'dt .02',
                      'ground_y 0', 'friction_options .3 .6', 'restitution 0', 'beta 0',
                      'positions '+values(initial_x), 'velocities '+values(initial_v),
                      'expected_positions '+values(x), 'expected_velocities '+values(v),
                      'compare_counts '+str(int(name != 'zero_tangent'))])
        emit_contact_info(lines, 'source', info)
        lines.append('end_ground')

    lines.append('penalty_cases 2')
    maximum_penalty_fd_error = 0.
    for name, y in [('mixed', [-.01, .02, -.004]), ('separated', [.01, .02, .004])]:
        x = np.array([[.1, y[0], .2], [.2, y[1], .3], [.3, y[2], .4]])
        masses, fixed, dt, eta = np.array([2., 3., 5.]), [1], .02, 30.
        energy, gradient, diagonal, impulse, pen = penalty(x, masses, fixed, dt, 0., eta)
        numerical = np.empty(3)
        numeric_diagonal = np.empty(3)
        for node in range(3):
            plus, minus = x.copy(), x.copy()
            plus[node, 1] += 1.e-6
            minus[node, 1] -= 1.e-6
            ep, gp, _, _, _ = penalty(plus, masses, fixed, dt, 0., eta)
            em, gm, _, _, _ = penalty(minus, masses, fixed, dt, 0., eta)
            numerical[node] = (ep-em)/2.e-6
            numeric_diagonal[node] = (gp[node]-gm[node])/2.e-6
        np.testing.assert_allclose(gradient, numerical, atol=1.e-7, rtol=1.e-10)
        np.testing.assert_allclose(diagonal, numeric_diagonal, atol=1.e-6, rtol=1.e-10)
        maximum_penalty_fd_error = max(maximum_penalty_fd_error, float(np.max(abs(gradient-numerical))))
        lines.extend(['penalty '+name, 'vertices 3', 'masses '+values(masses), fixed_record(fixed),
                      'dt '+values(dt), 'ground_y 0', 'eta '+values(eta), 'positions '+values(x),
                      'energy '+values(energy), 'gradient_y '+values(gradient),
                      'diagonal_y '+values(diagonal), 'normal_impulses '+values(impulse),
                      'penetration '+values(pen), 'end_penalty'])

    audit_summaries = []
    audits = [('fixed_below_ground', -.01, [0., -.2, 0.], [0], 0.),
              ('separating_velocity', 0., [.1, .5, .2], [], 0.),
              ('beta_reported_impulse', -.01, [.1, -.2, .2], [], .2)]
    lines.append('legacy_audits '+str(len(audits)))
    for name, height, velocity, fixed, beta in audits:
        x, v = np.array([[.2, height, .1]]), np.array([velocity], dtype=float)
        old_x, old_v = x.copy(), v.copy()
        info = apply_ground_contact(x, v, 0., .3, fixed, dt=.01, masses=np.array([2.]),
                                    static_friction=.6, restitution=0., contact_beta=beta)
        lines.extend(['audit '+name, 'vertices 1', 'masses 2', fixed_record(fixed), 'dt .01',
                      'ground_y 0', 'friction_options .3 .6', 'restitution 0', 'beta '+values(beta),
                      'positions '+values(old_x), 'velocities '+values(old_v),
                      'legacy_positions '+values(x), 'legacy_velocities '+values(v)])
        emit_contact_info(lines, 'legacy', info)
        actual_normal_change = float(2.*(v[0, 1]-old_v[0, 1]))
        lines.extend(['actual_normal_momentum_change '+values(actual_normal_change), 'end_audit'])
        audit_summaries.append({'name': name, 'reported_normal_impulse': info['normal_impulse'],
                                'actual_normal_momentum_change': actual_normal_change})

    lines.append('trajectories 4')
    maximum_final_residual = 0.
    trajectory_summaries = []
    for method in ['implicit_lumped', 'implicit_consistent', 'distance', 'energy']:
        mesh = ClothMesh(1, 1, 1., .7)
        theta = .2
        mesh.x = np.column_stack((mesh.X[:, 0], .001+np.sin(theta)*mesh.X[:, 1],
                                  np.cos(theta)*mesh.X[:, 1]))
        material = NeoHookeanMaterial(young=1000., poisson=.3, thickness=.001)
        physical = assemble_mass(mesh.X, mesh.triangles, 800., .001, mode='lumped')[::3]
        fixed, gravity = np.array([2, 3], dtype=int), np.array([0., -9.81, 0.])
        x, v = mesh.x.copy(), np.zeros_like(mesh.x)
        v[:2] = [.3, -.4, .1]
        dt, steps, eta_original, eta_port = .01, 8, 30., 30.
        mu, mu_static = .25, .4
        original_node_masses = physical.copy()
        if method.startswith('implicit'):
            mode = method.split('_')[1]
            solver = ImplicitEulerSolverSparse(mesh.X, mesh.triangles, material, 800., .001, mode,
                                               gravity, 100, 1.e-10, ground_y=0., contact_eta=eta_original)
            original_node_masses = solver.mass[::3].copy()
            if mode == 'consistent':
                np.testing.assert_allclose(physical, 2.*original_node_masses, atol=1.e-15, rtol=1.e-14)
                eta_port = eta_original/2.
                mu = mu_static = 0.
            constraints = []
            iterations = 100
        else:
            mass = physical.copy()
            mass[fixed] = np.inf
            iterations = 30
            if method == 'distance':
                constraints = (build_structural_constraints(mesh, 1.e-6)+
                               build_shear_constraints(mesh, 1.e-6)+build_bending_constraints(mesh, 1.e-5))
                solver = XPBDSolver(mass, constraints, gravity, iterations)
            else:
                scene = SimpleNamespace(mesh=mesh, material=material, lumped_mass=np.repeat(physical, 3),
                                        n_vertices=len(mesh.X), fixed_nodes=fixed, gravity=gravity)
                solver = NeoHookeanXPBDSolver(scene, 1.e-6, .2, iterations)
                constraints = solver.constraints
        lines.extend(['trajectory '+method, 'vertices 4', 'triangles 2', 'reference '+values(mesh.X),
                      'rest_positions '+values(mesh.x),
                      'topology '+' '.join(str(int(i)) for i in mesh.triangles.ravel()),
                      'material '+values([material.young, material.poisson, material.t, material.j_min]),
                      'density 800', 'physical_masses '+values(physical),
                      'original_node_masses '+values(original_node_masses),
                      'original_eta '+values(eta_original), 'port_eta '+values(eta_port),
                      'friction_options '+values([mu, mu_static]), 'ground_y 0', 'restitution 0', 'beta 0',
                      'dt '+values(dt), 'iterations '+str(iterations), 'newton_tolerance 1e-10',
                      'distance_compliances 1e-6 1e-6 1e-5', 'energy_options 1e-6 .2',
                      'gravity '+values(gravity), fixed_record(fixed), 'positions '+values(x),
                      'velocities '+values(v), 'steps '+str(steps)])
        contact_steps, total_dissipation = 0, 0.
        for step in range(1, steps+1):
            old_x, old_v = x.copy(), v.copy()
            x, v, info = solver.step(x, v, fixed, dt)
            if not info['converged'] or not np.isfinite(x).all() or not np.isfinite(v).all():
                raise RuntimeError(method+' step '+str(step)+': invalid source solve '+str(info))
            pre_x, pre_v = x.copy(), v.copy()
            if method.startswith('implicit'):
                energy, gy, _, impulses, pen = penalty(x, original_node_masses, fixed, dt, 0., eta_original)
                np.testing.assert_allclose(impulses, info['contact_normal_impulse'], atol=1.e-13, rtol=1.e-12)
                _, elastic_gradient = solver.energy_gradient(x)
                gradient = elastic_gradient.copy()
                gradient[1::3] += gy
                inertia = np.asarray(solver.M @ (x-old_x-dt*old_v).ravel()).ravel()
                scaled_force, scaled_gradient = dt*dt*solver.f_ext_flat, dt*dt*gradient
                residual = inertia-scaled_force+scaled_gradient
                free = np.ones(residual.size, dtype=bool)
                for node in fixed:
                    free[3*node:3*node+3] = False
                residual_norm = float(np.linalg.norm(residual[free]))
                threshold = 1.e-10*max(1., np.linalg.norm(inertia[free]), np.linalg.norm(scaled_force[free]),
                                      np.linalg.norm(scaled_gradient[free]))
                if residual_norm > threshold*(1+1.e-8):
                    raise RuntimeError(method+': final independent contact residual failed')
                maximum_final_residual = max(maximum_final_residual, residual_norm)
                contact = apply_friction_from_normal_impulse(v, impulses, pen, original_node_masses,
                                                             mu, mu_static, fixed, dt)
                kinetic_change = float(.5*(v.ravel() @ (solver.M @ v.ravel())-
                                             pre_v.ravel() @ (solver.M @ pre_v.ravel())))
                raw_before = raw_after = compliant_before = compliant_after = np.zeros(4)
                multipliers = np.empty(0)
            else:
                def diagnostics(points):
                    raw, corrected = np.zeros(4), np.zeros(4)
                    for c in constraints:
                        value = c.evaluate(points)
                        value = float(value if np.isscalar(value) else value[0])
                        group = 3 if method == 'energy' else (2 if hasattr(c, 'rest_angle') else
                                1 if c.__class__.__name__ == 'ShearConstraint' else 0)
                        raw[group] = max(raw[group], abs(value))
                        corrected[group] = max(corrected[group], abs(value+c.compliance/(dt*dt)*c.lambda_value))
                    return raw, corrected
                raw_before, compliant_before = diagnostics(x)
                active = (x[:, 1] < 0.) | np.isclose(x[:, 1], 0., atol=1.e-12)
                active[fixed] = False
                if np.any(v[active, 1] > 1.e-12):
                    raise RuntimeError(method+': trajectory entered legacy separating-velocity bug domain')
                pen = np.maximum(-x[:, 1], 0.)
                impulses = np.where(active, physical*np.maximum(-v[:, 1], 0.), 0.)
                contact = apply_ground_contact(x, v, 0., mu, fixed, dt, physical, mu_static,
                                               restitution=0., contact_beta=0.)
                raw_after, compliant_after = diagnostics(x)
                multipliers = np.array([c.lambda_value for c in constraints])
                residual_norm, threshold, energy = 0., 0., 0.
                kinetic_change = float(.5*np.sum(physical[:, None]*(v*v-pre_v*pre_v)))
            np.testing.assert_array_equal(x[fixed], old_x[fixed])
            np.testing.assert_array_equal(v[fixed], np.zeros((len(fixed), 3)))
            assert np.min(x[fixed, 1]) > 0. and contact['friction_dissipation'] >= -1.e-14
            if not method.startswith('implicit'):
                np.testing.assert_allclose(impulses.sum(), contact['normal_impulse'], atol=1.e-13, rtol=1.e-12)
                assert np.min(x[:, 1]) >= -1.e-14
            lines.extend(['step '+str(step), 'pre_contact_positions '+values(pre_x),
                          'pre_contact_velocities '+values(pre_v), 'positions '+values(x), 'velocities '+values(v),
                          'penetration '+values(pen), 'normal_impulses '+values(impulses),
                          'penalty_energy '+values(energy), 'residual_norm '+values(residual_norm),
                          'residual_threshold '+values(threshold),
                          'physical_mass_kinetic_change '+values(kinetic_change),
                          'multipliers '+str(len(multipliers))+' '+values(multipliers),
                          'raw_before '+values(raw_before), 'compliant_before '+values(compliant_before),
                          'raw_after '+values(raw_after), 'compliant_after '+values(compliant_after)])
            emit_contact_info(lines, 'source', contact)
            lines.append('end_step')
            contact_steps += int(contact['contact_count'] > 0)
            total_dissipation += contact['friction_dissipation']
        if contact_steps == 0:
            raise RuntimeError(method+': no actual contact in reference trajectory')
        trajectory_summaries.append({'name': method, 'steps': steps, 'contact_steps': contact_steps,
                                     'original_eta': eta_original, 'port_eta': eta_port,
                                     'friction': mu, 'total_friction_dissipation': total_dissipation})
        lines.append('end_trajectory')
    lines.append('end')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text('\n'.join(lines)+'\n', encoding='utf-8')
    manifest = {'task': '05B', 'fixture_version': 1, 'numpy': np.__version__, 'scipy': scipy.__version__,
                'source_sha256': hashes, 'fixture_sha256': hashlib.sha256(args.output.read_bytes()).hexdigest(),
                'counts': {'friction_cases': 4, 'ground_cases': 5, 'penalty_cases': 2,
                           'legacy_audits': 3, 'trajectories': 4, 'time_steps': 32},
                'maximum_penalty_gradient_fd_error': maximum_penalty_fd_error,
                'maximum_final_implicit_residual': maximum_final_residual,
                'legacy_audits': audit_summaries, 'trajectories': trajectory_summaries,
                'consistent_policy': 'physical row-sum=2*original diagonal masses checked; eta_port=eta_original/2; friction=0; parameter-mapped comparison only',
                'count_policy': 'zero-vt sticking counts from contact.py and friction.py intentionally differ; no cross-port equality assertion for those counts'}
    args.output.with_suffix('.sources.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({key: manifest[key] for key in ('counts', 'maximum_final_implicit_residual', 'trajectories')}, indent=2))
    print(args.output.resolve())


if __name__ == '__main__':
    main()
