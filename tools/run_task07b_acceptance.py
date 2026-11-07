"""Run the bounded task-07B arbitrary-topology acceptance."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

from run_task07_acceptance import check_ctest, check_fixture_sources, check_python, run

PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_ORIGINAL = PROJECT.parents[1] / "python" / "job" / "cloth_sim"
DEFAULT_BUILD = PROJECT / "build" / "cpu-release"


def check_custom_scene(build: Path, output_root: Path) -> dict:
    executable = build / ("cloth_scene_run.exe" if sys.platform.startswith("win") else "cloth_scene_run")
    scene = PROJECT / "tests" / "scenes" / "task07_custom_mesh.scene"
    report = output_root / "task07_custom_mesh.txt"
    result = run([str(executable), "--scene", str(scene), "--output", str(report)], PROJECT, 30.0)
    text = report.read_text(encoding="utf-8") if report.is_file() else ""
    scene_text = scene.read_text(encoding="utf-8")
    expected_vertices = len(re.findall(r"^mesh_position\.\d+=", scene_text, re.MULTILINE))
    expected_triangles = len(re.findall(r"^mesh_triangle\.\d+=", scene_text, re.MULTILINE))
    status = re.search(r"^status=(.+)$", text, re.MULTILINE)
    committed = re.search(r"^committed_steps=(\d+)$", text, re.MULTILINE)
    frames = len(re.findall(r"^frame \d+ ", text, re.MULTILINE))
    return {
        "scene": scene.name,
        "run": result,
        "status": status.group(1) if status else None,
        "committed_steps": int(committed.group(1)) if committed else None,
        "frame_count": frames,
        "arbitrary_topology": expected_vertices == 4 and expected_triangles == 2
        and bool(re.search(r"^frame 1 0\.01 4$", text, re.MULTILINE)),
        "passed": result["returncode"] == 0 and status and status.group(1) == "committed"
        and committed and frames == int(committed.group(1))
        and expected_vertices == 4 and expected_triangles == 2
        and bool(re.search(r"^frame 1 0\.01 4$", text, re.MULTILINE)),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original", type=Path, default=DEFAULT_ORIGINAL)
    parser.add_argument("--build", type=Path, default=DEFAULT_BUILD)
    parser.add_argument("--output", type=Path, default=PROJECT / "output" / "task07b_acceptance.json")
    args = parser.parse_args()
    original = args.original.resolve(); build = args.build.resolve(); output = args.output.resolve()
    if not original.is_dir() or not build.is_dir():
        raise SystemExit("original project and Release build directory are required")
    output.parent.mkdir(parents=True, exist_ok=True)
    temp_root = PROJECT / "tmp" / "task07b_scene_acceptance"
    if temp_root.exists():
        shutil.rmtree(temp_root)
    temp_root.mkdir(parents=True, exist_ok=True)
    custom = check_custom_scene(build, temp_root)
    report = {
        "task": "07B", "date": "2026-11-05", "scope": "arbitrary-topology scene migration",
        "original_project": str(original), "build": str(build),
        "python_validation": check_python(original),
        "cpp_release_ctest": check_ctest(build),
        "custom_mesh_scene": custom,
        "task05b_fixture_sources": check_fixture_sources(original),
        "limitations": ["custom topology is validated for triangle meshes and distance XPBD only",
                        "C01-C08 complex scenarios remain planned", "GPU/CUDA not executed",
                        "cloud execution not executed", "Linux/MSVC/sanitizer/Docker not executed"],
    }
    report["passed"] = all(report[key]["passed"] for key in
                            ("python_validation", "cpp_release_ctest", "custom_mesh_scene",
                             "task05b_fixture_sources"))
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"passed": report["passed"],
                      "python_tests": report["python_validation"]["passed_tests"],
                      "cpp_tests": report["cpp_release_ctest"]["tests"],
                      "custom_mesh": custom["passed"],
                      "fixture_sources": report["task05b_fixture_sources"]["count"],
                      "output": str(output)}, ensure_ascii=False, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
