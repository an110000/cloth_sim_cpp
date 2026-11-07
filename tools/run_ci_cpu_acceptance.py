"""Run the repository-only CPU acceptance used by Linux CI.

The original Python project is intentionally not required here: this report
checks the C++ repository, its CTest result, and the arbitrary-topology scene.
The cross-project Python comparison remains a separate local acceptance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

PROJECT = Path(__file__).resolve().parents[1]


def run(command: list[str], cwd: Path, timeout: float) -> dict:
    started = time.perf_counter()
    try:
        result = subprocess.run(command, cwd=cwd, text=True, capture_output=True,
                                timeout=timeout, check=False)
        return {"command": command, "returncode": result.returncode,
                "seconds": time.perf_counter() - started,
                "stdout": result.stdout[-8000:], "stderr": result.stderr[-8000:]}
    except subprocess.TimeoutExpired as error:
        return {"command": command, "returncode": None,
                "seconds": time.perf_counter() - started,
                "stdout": str(error.stdout or "")[-8000:],
                "stderr": (str(error.stderr or "") + " timeout").strip()[-8000:]}


def sha256(path: Path) -> str | None:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def probe(command: list[str], cwd: Path) -> str:
    result = run(command, cwd, 20.0)
    text = (result["stdout"] + result["stderr"]).strip()
    if result["returncode"] != 0:
        return "<unavailable>"
    return text[-2000:]


def ctest_summary(build: Path, junit: Path) -> dict:
    details = {"tests": None, "failures": None, "skipped": None}
    if junit.is_file():
        root = ET.parse(junit).getroot()
        details = {"tests": int(root.attrib.get("tests", "-1")),
                   "failures": int(root.attrib.get("failures", "-1")),
                   "skipped": int(root.attrib.get("skipped", "0"))}
    return details


def run_custom_scene(build: Path, output_root: Path) -> dict:
    executable = build / ("cloth_scene_run.exe" if os.name == "nt" else "cloth_scene_run")
    scene = PROJECT / "tests" / "scenes" / "task07_custom_mesh.scene"
    report = output_root / "task07_custom_mesh.txt"
    result = run([str(executable), "--scene", str(scene), "--output", str(report)], PROJECT, 30.0)
    text = report.read_text(encoding="utf-8") if report.is_file() else ""
    status = re.search(r"^status=(.+)$", text, re.MULTILINE)
    requested = re.search(r"^requested_steps=(\d+)$", text, re.MULTILINE)
    committed = re.search(r"^committed_steps=(\d+)$", text, re.MULTILINE)
    frame = re.search(r"^frame 1 0\.01 (\d+)$", text, re.MULTILINE)
    scene_text = scene.read_text(encoding="utf-8")
    vertices = len(re.findall(r"^mesh_position\.\d+=", scene_text, re.MULTILINE))
    triangles = len(re.findall(r"^mesh_triangle\.\d+=", scene_text, re.MULTILINE))
    return {
        "scene": scene.name,
        "scene_sha256": sha256(scene),
        "run": result,
        "status": status.group(1) if status else None,
        "requested_steps": int(requested.group(1)) if requested else None,
        "committed_steps": int(committed.group(1)) if committed else None,
        "vertex_count": int(frame.group(1)) if frame else None,
        "expected_vertices": vertices,
        "expected_triangles": triangles,
        "passed": result["returncode"] == 0 and status and status.group(1) == "committed"
        and requested and committed and requested.group(1) == committed.group(1)
        and frame and int(frame.group(1)) == vertices == 4 and triangles == 2,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build = args.build.resolve(); output = args.output.resolve()
    if not build.is_dir():
        raise SystemExit(f"build directory does not exist: {build}")
    if output.exists():
        shutil.rmtree(output)
    output.mkdir(parents=True)
    junit = build / "ci-ctest-results.xml"
    ctest = run(["ctest", "--output-on-failure", "--output-junit", str(junit)], build, 300.0)
    # CTest resolves a relative JUnit path relative to its working directory.
    if not junit.is_file() and (build / str(junit)).is_file():
        junit = build / str(junit)
    summary = ctest_summary(build, junit)
    custom = run_custom_scene(build, output)
    commit = probe(["git", "rev-parse", "HEAD"], PROJECT)
    status = probe(["git", "status", "--porcelain"], PROJECT)
    report = {
        "task": "07C", "date": "2026-11-06", "scope": "repository-only Linux/CPU CI acceptance",
        "project": str(PROJECT), "build": str(build),
        "source_identity": {"git_commit": commit, "git_status_porcelain": status},
        "environment": {
            "os": platform.platform(), "machine": platform.machine(),
            "python": sys.version, "cmake": probe(["cmake", "--version"], PROJECT),
            "compiler": probe(["c++", "--version"], PROJECT),
        },
        "cpp_ctest": {"run": ctest, **summary,
                      "junit_sha256": sha256(junit)},
        "custom_mesh_scene": custom,
        "limitations": [
            "original Python cross-project validation is a separate local acceptance",
            "GPU/CUDA not executed", "complex C01-C08 scenarios remain planned",
            "this job does not claim Linux CI passed until a remote runner uploads this report",
        ],
    }
    report["passed"] = ctest["returncode"] == 0 and summary["failures"] == 0 and custom["passed"]
    report_path = output / "task07c_cpu_acceptance.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"passed": report["passed"], "ctest_tests": summary["tests"],
                      "custom_mesh": custom["passed"], "report": str(report_path)}, ensure_ascii=False, indent=2))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
