"""Export analytic derivatives and assembled systems from the original Python code.

Task-01 fixtures are read, never regenerated or overwritten by this script.
The checked-in output lets CTest run without Python, NumPy, or SciPy.
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


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", required=True, type=Path)
    args = parser.parse_args()
    source = args.source.resolve()
    modules = ["core/mesh.py", "core/energy.py", "core/material.py", "core/mass.py",
               "core/mass_sparse.py", "solvers/implicit.py", "solvers/implicit_sparse.py"]
    for relative in modules:
        if not (source / relative).is_file():
            parser.error("Missing reference source: " + relative)
    sys.path.insert(0, str(source))
    from core.energy import deformation_gradient, triangle_energy, triangle_gradient, triangle_hessian
    from core.mass import assemble_mass
    from core.material import NeoHookeanMaterial
    from core.mesh import ClothMesh
    from solvers.implicit import ImplicitEulerSolver
    from solvers.implicit_sparse import ImplicitEulerSolverSparse

    lines = ["# Generated from the original Python implementation; review changes before accepting.",
             "# NumPy " + np.__version__ + "; SciPy " + scipy.__version__]
    for relative in modules:
        digest = hashlib.sha256((source / relative).read_bytes()).hexdigest()
        lines.append("# sha256 " + relative + " " + digest)
    task01 = root / "tests/data/python_reference.txt"
    lines.append("# sha256 task01_fixture " + hashlib.sha256(task01.read_bytes()).hexdigest())
    cases = []
    for line in task01.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        tokens = line.split()
        numbers = np.array([float(token) for token in tokens[1:]])
        cases.append((tokens[0], numbers[:4], numbers[4:10].reshape(3, 2),
                      numbers[10:19].reshape(3, 3)))
    if len(cases) != 11:
        raise RuntimeError("Expected the frozen 11 task-01 input cases")
    lines.extend(["CLOTH_TASK02_REFERENCE 1", "elements " + str(len(cases))])
    for name, parameters, X, x in cases:
        young, poisson, thickness, j_min = parameters
        material = NeoHookeanMaterial(young, poisson, j_min, thickness=thickness)
        F = deformation_gradient(X, x)
        lines.extend([
            "element " + name, "material " + values(parameters),
            "reference " + values(X), "positions " + values(x),
            "energy " + values(triangle_energy(X, x, material)),
            "gradient " + values(triangle_gradient(X, x, material)),
            "area_ratio " + values(np.sqrt(np.linalg.det(F.T @ F))),
            "hessian " + values(triangle_hessian(X, x, material)), "end",
        ])

    square = ClothMesh(1, 1, 1.0, 1.0)
    strip = ClothMesh(2, 1, 1.0, 0.8)
    irregular_X = np.array([[0., 0.], [1.2, 0.], [.1, .9], [1.1, 1.3], [2., .2]])
    irregular_triangles = np.array([[0, 1, 2], [1, 3, 2], [1, 4, 3]], dtype=int)
    systems = [("two_triangles", square.X, square.triangles, [0, 2], False),
               ("nonuniform", irregular_X, irregular_triangles, [1], False),
               ("regularized_strip", strip.X, strip.triangles,
                list(range(len(strip.X))), True)]
    lines.append("systems " + str(len(systems)))
    for name, X, triangles, fixed, regularized in systems:
        parameters = [8.e5, .27, .0012, .05]
        density = 720.0
        x = np.column_stack((1.07 * X[:, 0] + .08 * X[:, 1],
                             .03 * np.sin(X[:, 0]) - .05 * X[:, 1],
                             .96 * X[:, 1] + .02 * X[:, 0]))
        if regularized:
            x = np.column_stack((.025 * X[:, 0], np.zeros(len(X)), X[:, 1]))
        material = NeoHookeanMaterial(parameters[0], parameters[1], parameters[3],
                                     thickness=parameters[2])
        solver = ImplicitEulerSolverSparse(X, triangles, material, density=density,
                                          thickness=parameters[2], gravity=(0., 0., 0.))
        energy, gradient, hessian = solver.assemble_energy_gradient_hessian(x)
        hessian = hessian.toarray()
        # An independent dense assembly must agree before accepting the fixture.
        with contextlib.redirect_stdout(io.StringIO()):
            dense_solver = ImplicitEulerSolver(X, triangles, material, density=density,
                                              thickness=parameters[2], gravity=(0., 0., 0.))
            dense_energy, dense_gradient, dense_hessian = dense_solver.assemble_energy_gradient_hessian(x)
        np.testing.assert_allclose(energy, dense_energy, rtol=2.e-10, atol=1.e-9)
        np.testing.assert_allclose(np.asarray(gradient).ravel(), np.asarray(dense_gradient).ravel(),
                                   rtol=2.e-10, atol=1.e-9)
        np.testing.assert_allclose(hessian, dense_hessian, rtol=2.e-10, atol=2.e-8)
        lumped = np.diag(assemble_mass(X, triangles, density, parameters[2], "lumped"))
        consistent = assemble_mass(X, triangles, density, parameters[2], "consistent")
        free = np.array([3 * i + axis for i in range(len(X)) if i not in fixed
                         for axis in range(3)], dtype=int)
        lines.extend([
            "system " + name, "vertices " + str(len(X)),
            "reference " + values(X), "positions " + values(x),
            "triangles " + str(len(triangles)) + " " + " ".join(str(i) for i in triangles.ravel()),
            "material " + values(parameters), "density " + values(density),
            "fixed " + str(len(fixed)) + " " + " ".join(str(i) for i in fixed),
            "energy " + values(energy), "gradient " + values(gradient),
            "hessian " + values(hessian), "mass_lumped " + values(lumped),
            "mass_consistent " + values(consistent),
            "free_indices " + str(len(free)) + " " + " ".join(str(i) for i in free),
            "hessian_free " + values(hessian[np.ix_(free, free)]), "end",
        ])
    lines.append("done")
    output = root / "tests/data/task02_reference.txt"
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("Exported", len(cases), "element Hessians and", len(systems), "assembled systems to", output)


if __name__ == "__main__":
    main()
