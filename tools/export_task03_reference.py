"""Freeze no-contact implicit-Euler trajectories from the original Python solver.

Every recorded step must report convergence and pass a freshly recomputed final
residual check. The original modules and the task-01/02 fixtures are not edited.
Normal C++ builds consume the checked-in text and do not require Python.
"""
import argparse
import contextlib
import hashlib
import io
from pathlib import Path
import sys

import numpy as np
import scipy


def values(array):
    return " ".join(format(float(value), ".17g") for value in np.asarray(array).ravel())


def final_residual(solver, previous_x, previous_v, x, fixed, dt):
    """Independently evaluate the no-contact stationarity test at returned x."""
    y = (previous_x + dt * previous_v).ravel()
    _, gradient = solver.energy_gradient(x)
    inertia = np.asarray(solver.M @ (x.ravel() - y)).ravel()
    gravity = dt * dt * solver.f_ext_flat
    elastic = dt * dt * np.asarray(gradient).ravel()
    residual = inertia - gravity + elastic
    mask = np.ones(residual.size, dtype=bool)
    for vertex in fixed:
        mask[3 * vertex:3 * vertex + 3] = False
    scale = max(1.0, np.linalg.norm(inertia[mask]),
                np.linalg.norm(gravity[mask]), np.linalg.norm(elastic[mask]))
    return float(np.linalg.norm(residual[mask])), float(solver.newton_tol * scale)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, type=Path)
    parser.add_argument("--output", type=Path, default=root / "tests/data/task03_reference.txt")
    args = parser.parse_args()
    source = args.source.resolve()
    modules = ["core/mesh.py", "core/energy.py", "core/material.py", "core/mass.py",
               "core/mass_sparse.py", "solvers/implicit.py", "solvers/implicit_sparse.py"]
    for relative in modules:
        if not (source / relative).is_file():
            parser.error("Missing original source: " + relative)
    sys.path.insert(0, str(source))
    from core.material import NeoHookeanMaterial
    from core.mesh import ClothMesh
    from solvers.implicit import ImplicitEulerSolver as DenseSolver
    from solvers.implicit_sparse import ImplicitEulerSolverSparse as SparseSolver

    parameters = np.array([2.e4, .3, .001, .05])
    density, dt, tolerance, max_iter, steps = 800.0, .01, 1.e-10, 50, 8
    mesh = ClothMesh(2, 1, 1.0, .7)
    X, triangles = mesh.X.copy(), mesh.triangles.copy()
    base = np.column_stack((X[:, 0], np.zeros(len(X)), X[:, 1]))
    names = [(kind, mode) for kind in ("rest", "fall", "pinned")
             for mode in ("lumped", "consistent")]
    lines = ["CLOTH_TASK03_REFERENCE 1",
             "# Original Python sparse solver, no ground contact, no friction, no bending.",
             "# NumPy " + np.__version__ + "; SciPy " + scipy.__version__,
             "# residual_norm/threshold are recomputed at each returned state, not copied from info.",
             "# pinned_lumped additionally checked against the original dense solver at every step."]
    for relative in modules:
        lines.append("# sha256 " + relative + " " +
                     hashlib.sha256((source / relative).read_bytes()).hexdigest())
    lines.append("cases " + str(len(names)))
    total_steps, largest_residual, dense_checks = 0, 0.0, 0
    for kind, mode in names:
        name = kind + "_" + mode
        fixed = [3, 4, 5] if kind == "pinned" else []
        gravity = np.array([0., 0., 0.] if kind == "rest" else [0., -9.81, 0.])
        x, v = base.copy(), np.zeros_like(base)
        if kind == "fall":
            v[:] = [.12, .07, -.03]
        elif kind == "pinned":
            x[:3, 0] = 1.06 * x[:3, 0] + .025
            x[:3, 1] = -.025 * (1.0 + X[:3, 0])
            x[:3, 2] -= .025
            v[:3] = [.015, -.02, .01]
        material = NeoHookeanMaterial(parameters[0], parameters[1], parameters[3],
                                     thickness=parameters[2])
        kwargs = dict(X=X, triangles=triangles, material=material, density=density,
                      thickness=parameters[2], mass_mode=mode, gravity=gravity,
                      newton_max_iter=max_iter, newton_tol=tolerance)
        solver = SparseSolver(**kwargs)
        dense = None
        if name == "pinned_lumped":
            with contextlib.redirect_stdout(io.StringIO()):
                dense = DenseSolver(**kwargs)
            dense_x, dense_v = x.copy(), v.copy()
        lines.extend(["case " + name, "vertices " + str(len(X)),
                      "triangles " + str(len(triangles)), "reference " + values(X),
                      "topology " + " ".join(str(int(i)) for i in triangles.ravel()),
                      "material " + values(parameters), "density " + values(density),
                      "mass_mode " + mode, "gravity " + values(gravity),
                      "dt " + values(dt), "newton " + values(tolerance) + " " + str(max_iter),
                      "fixed " + str(len(fixed)) + "".join(" " + str(i) for i in fixed),
                      "positions " + values(x), "velocities " + values(v),
                      "steps " + str(steps)])
        for index in range(1, steps + 1):
            old_x, old_v = x.copy(), v.copy()
            x, v, info = solver.step(old_x, old_v, fixed, dt)
            residual, threshold = final_residual(solver, old_x, old_v, x, fixed, dt)
            if not info["converged"] or not np.isfinite(x).all() or not np.isfinite(v).all():
                raise RuntimeError(name + " step " + str(index) + " did not converge: " + str(info))
            if not np.isfinite(residual) or residual > threshold * (1.0 + 1.e-8):
                raise RuntimeError(name + " returned residual fails convergence: " +
                                   str((index, residual, threshold)))
            if fixed:
                np.testing.assert_array_equal(x[fixed], old_x[fixed])
                np.testing.assert_array_equal(v[fixed], np.zeros((len(fixed), 3)))
            if kind == "rest":
                np.testing.assert_allclose(x, base, rtol=0., atol=1.e-13)
                np.testing.assert_allclose(v, 0., rtol=0., atol=1.e-12)
            if kind == "fall":
                expected_v = np.tile([.12, .07, -.03], (len(X), 1)) + index * dt * gravity
                expected_x = (base + index * dt * np.tile([.12, .07, -.03], (len(X), 1)) +
                              dt * dt * index * (index + 1) / 2.0 * gravity)
                np.testing.assert_allclose(x, expected_x, rtol=0., atol=2.e-12)
                np.testing.assert_allclose(v, expected_v, rtol=0., atol=2.e-11)
            if dense is not None:
                dense_x, dense_v, dense_info = dense.step(dense_x, dense_v, fixed, dt)
                if not dense_info["converged"]:
                    raise RuntimeError("Dense reference did not converge at step " + str(index))
                np.testing.assert_allclose(x, dense_x, rtol=2.e-8, atol=2.e-9)
                np.testing.assert_allclose(v, dense_v, rtol=2.e-7, atol=2.e-7)
                dense_checks += 1
            lines.extend(["step " + str(index), "positions " + values(x),
                          "velocities " + values(v), "residual_norm " + values(residual),
                          "residual_threshold " + values(threshold),
                          "iterations " + str(info["iterations"]), "end_step"])
            total_steps += 1
            largest_residual = max(largest_residual, residual)
        lines.append("end_case")
    lines.append("end")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("Exported", len(names), "trajectories /", total_steps, "steps;",
          dense_checks, "dense comparisons; max final residual", largest_residual)
    print(args.output.resolve())


if __name__ == "__main__":
    main()
