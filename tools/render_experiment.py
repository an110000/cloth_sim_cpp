"""Render a headless C++ experiment report into PNG/GIF diagnostics.

This is an analysis tool, not part of the C++ solver. It consumes the stable
CLOTH_EXPERIMENT_TEXT report emitted by cloth_scene_run.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np


def parse_report(path: Path) -> tuple[str, list[dict]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    status = next((line.split("=", 1)[1] for line in lines if line.startswith("status=")), "unknown")
    frames: list[dict] = []
    index = 0
    while index < len(lines):
        match = re.match(r"^frame (\d+) ([^ ]+) (\d+)$", lines[index])
        if not match:
            index += 1
            continue
        step, time, count = int(match.group(1)), float(match.group(2)), int(match.group(3))
        positions = []; velocities = []
        for node in range(count):
            index += 1
            fields = lines[index].split()
            if len(fields) != 8 or fields[0] != "node" or int(fields[1]) != node:
                raise ValueError(f"Malformed node line near frame {step}: {lines[index]}")
            positions.append([float(value) for value in fields[2:5]])
            velocities.append([float(value) for value in fields[5:8]])
        frames.append({"step": step, "time": time,
                       "positions": np.asarray(positions),
                       "velocities": np.asarray(velocities)})
        index += 1
    if not frames:
        raise ValueError(f"No committed frames found in {path}")
    return status, frames


def parse_scene(path: Path | None, vertex_count: int) -> tuple[list[tuple[int, int, int]], list[int]]:
    if path is None:
        return [], []
    values = {}
    fixed: list[int] = []
    nx = ny = None
    base = "generated"
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.split("#", 1)[0].strip()
        if "=" not in line:
            continue
        key, value = (part.strip() for part in line.split("=", 1))
        if key == "base": base = value
        elif key == "fixed_vertices":
            fixed = [int(part.strip()) for part in value.split(",") if part.strip()]
        elif key == "nx": nx = int(value)
        elif key == "ny": ny = int(value)
        elif key.startswith("mesh_triangle."):
            values[int(key.split(".", 1)[1])] = tuple(int(part.strip()) for part in value.split(","))
    if base == "teaching" and not fixed:
        fixed = [3, 4, 5]
    if base == "teaching" and nx is None:
        nx, ny = 2, 1
    triangles = [values[key] for key in sorted(values)]
    if not triangles and nx and ny and vertex_count == (nx + 1) * (ny + 1):
        stride = nx + 1
        for row in range(ny):
            for col in range(nx):
                n0 = row * stride + col; n1 = n0 + 1
                n2 = n0 + stride; n3 = n2 + 1
                triangles.extend([(n0, n1, n3), (n0, n3, n2)])
    return triangles, fixed


def limits(frames: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    points = np.concatenate([frame["positions"] for frame in frames], axis=0)
    low = points.min(axis=0); high = points.max(axis=0)
    span = np.maximum(high - low, 1.0e-6)
    padding = 0.12 * span
    return low - padding, high + padding


def draw_frame(frame: dict, triangles: list[tuple[int, int, int]], fixed: list[int], low: np.ndarray,
               high: np.ndarray, path: Path, title: str) -> None:
    position = frame["positions"]; speed = np.linalg.norm(frame["velocities"], axis=1)
    fig = plt.figure(figsize=(8, 6), dpi=140)
    axis = fig.add_subplot(111, projection="3d")
    if triangles:
        faces = [[position[i], position[j], position[k]] for i, j, k in triangles]
        collection = Poly3DCollection(faces, alpha=0.58, linewidths=0.7,
                                      edgecolors="#264653", cmap="viridis")
        collection.set_array(np.asarray([speed[list(triangle)].mean() for triangle in triangles]))
        axis.add_collection3d(collection)
        fig.colorbar(collection, ax=axis, shrink=0.65, pad=0.08, label="speed")
    scatter = axis.scatter(position[:, 0], position[:, 1], position[:, 2], c=speed,
                           cmap="viridis", s=24, depthshade=True)
    if fixed:
        pinned = position[fixed]
        axis.scatter(pinned[:, 0], pinned[:, 1], pinned[:, 2],
                     c="#d62828", s=44, marker="s", label="fixed")
        axis.legend(loc="upper left")
    axis.set_xlim(low[0], high[0]); axis.set_ylim(low[1], high[1]); axis.set_zlim(low[2], high[2])
    axis.set_xlabel("x [m]"); axis.set_ylabel("y [m]"); axis.set_zlabel("z [m]")
    axis.set_title(title)
    fig.tight_layout(); fig.savefig(path); plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--scene", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--fps", type=int, default=12)
    args = parser.parse_args()
    status, frames = parse_report(args.report)
    if args.fps <= 0:
        raise SystemExit("--fps must be positive")
    args.output.mkdir(parents=True, exist_ok=True)
    triangles, fixed = parse_scene(args.scene, len(frames[0]["positions"]))
    low, high = limits(frames)
    for frame in (frames[0], frames[-1]):
        draw_frame(frame, triangles, fixed, low, high,
                   args.output / f"frame_{frame['step']:06d}.png",
                   f"step {frame['step']}  t={frame['time']:.4f}s")
    fig, axis = plt.subplots(2, 1, figsize=(8, 6), dpi=140, sharex=True)
    times = np.asarray([frame["time"] for frame in frames])
    speed = np.asarray([np.linalg.norm(frame["velocities"], axis=1).max() for frame in frames])
    origin = frames[0]["positions"]
    displacement = np.asarray([np.linalg.norm(frame["positions"] - origin, axis=1).max() for frame in frames])
    axis[0].plot(times, speed, "o-", label="max speed")
    axis[1].plot(times, displacement, "o-", label="max displacement")
    for current in axis:
        current.grid(alpha=0.25); current.legend(); current.set_ylabel("value")
    axis[-1].set_xlabel("time [s]")
    fig.suptitle("C++ experiment diagnostics"); fig.tight_layout()
    fig.savefig(args.output / "diagnostics.png"); plt.close(fig)
    animation_path = args.output / "animation.gif"
    fig = plt.figure(figsize=(8, 6), dpi=110)
    def update(frame: dict) -> None:
        fig.clear(); axis = fig.add_subplot(111, projection="3d")
        position = frame["positions"]; speed = np.linalg.norm(frame["velocities"], axis=1)
        if triangles:
            axis.add_collection3d(Poly3DCollection([[position[i], position[j], position[k]]
                for i, j, k in triangles], alpha=0.58, edgecolors="#264653"))
        axis.scatter(position[:, 0], position[:, 1], position[:, 2], c=speed, cmap="viridis", s=24)
        if fixed:
            pinned = position[fixed]; axis.scatter(pinned[:, 0], pinned[:, 1], pinned[:, 2], c="#d62828", marker="s")
        axis.set_xlim(low[0], high[0]); axis.set_ylim(low[1], high[1]); axis.set_zlim(low[2], high[2])
        axis.set_xlabel("x [m]"); axis.set_ylabel("y [m]"); axis.set_zlabel("z [m]")
        axis.set_title(f"step {frame['step']}  t={frame['time']:.4f}s")
    if len(frames) > 1:
        animation = FuncAnimation(fig, update, frames=frames, repeat=False)
        animation.save(animation_path, writer=PillowWriter(fps=args.fps))
    else:
        update(frames[0])
    if len(frames) == 1:
        fig.savefig(animation_path.with_suffix(".png"))
    plt.close(fig)
    summary = {"status": status, "frames": len(frames), "vertices": len(frames[0]["positions"]),
               "triangles": len(triangles), "fixed_vertices": fixed,
               "outputs": sorted(path.name for path in args.output.iterdir())}
    (args.output / "visualization_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
