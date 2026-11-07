"""Run the bounded task-07 migration acceptance without building or downloading."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_ORIGINAL = PROJECT.parents[1] / "python" / "job" / "cloth_sim"
DEFAULT_BUILD = PROJECT / "build" / "cpu-release"

def run(command: list[str], cwd: Path, timeout: float) -> dict:
    started = time.perf_counter()
    try:
        result = subprocess.run(command, cwd=cwd, text=True, capture_output=True,
                                timeout=timeout, check=False)
        return {"command": command, "returncode": result.returncode,
                "seconds": time.perf_counter() - started,
                "stdout": result.stdout[-4000:], "stderr": result.stderr[-4000:]}
    except subprocess.TimeoutExpired as error:
        return {"command": command, "returncode": None,
                "seconds": time.perf_counter() - started,
                "stdout": str(error.stdout or "")[-4000:],
                "stderr": (str(error.stderr or "") + " timeout").strip()[-4000:]}

def check_python(original: Path) -> dict:
    result = run([sys.executable, "-m", "pytest", "validation", "-q"], original, 300.0)
    match = re.search(r"(\d+) passed", result["stdout"] + result["stderr"])
    passed = int(match.group(1)) if match else None
    result.update({"passed_tests": passed, "passed": result["returncode"] == 0 and passed is not None})
    return result

def check_ctest(build: Path) -> dict:
    result = run(["ctest", "--output-junit", "task07-ctest-results.xml"], build, 180.0)
    junit = build / "task07-ctest-results.xml"
    details = {"tests": None, "failures": None, "skipped": None}
    if junit.is_file():
        root = ET.parse(junit).getroot()
        details = {"tests": int(root.attrib.get("tests", "-1")),
                   "failures": int(root.attrib.get("failures", "-1")),
                   "skipped": int(root.attrib.get("skipped", "0"))}
    result.update(details)
    result["passed"] = result["returncode"] == 0 and details["tests"] is not None and details["failures"] == 0
    return result

def check_scene_runs(build: Path, output_root: Path) -> dict:
    executable = build / ("cloth_scene_run.exe" if os.name == "nt" else "cloth_scene_run")
    scenes = [PROJECT / "tests" / "scenes" / name for name in
              ("task06_teaching.scene", "task06_implicit.scene", "task06_energy.scene")]
    records = []
    for scene in scenes:
        report = output_root / (scene.stem + ".txt")
        result = run([str(executable), "--scene", str(scene), "--output", str(report)], PROJECT, 30.0)
        text = report.read_text(encoding="utf-8") if report.is_file() else ""
        status = re.search(r"^status=(.+)$", text, re.MULTILINE)
        committed = re.search(r"^committed_steps=(\d+)$", text, re.MULTILINE)
        frames = len(re.findall(r"^frame \d+ ", text, re.MULTILINE))
        bad_number = bool(re.search(r"(?:^|[= ])(?:nan|inf)(?:$|\s)", text, re.IGNORECASE))
        item = {"scene": scene.name, "run": result, "status": status.group(1) if status else None,
                "committed_steps": int(committed.group(1)) if committed else None,
                "frame_count": frames, "finite_text": not bad_number,
                "passed": result["returncode"] == 0 and status and status.group(1) == "committed"
                and committed and frames == int(committed.group(1)) and not bad_number}
        records.append(item)
    return {"runs": records, "passed": all(item["passed"] for item in records)}

def check_fixture_sources(original: Path) -> dict:
    manifest_path = PROJECT / "tests" / "data" / "task05b_reference.sources.json"
    data = json.loads(manifest_path.read_text(encoding="utf-8"))
    mismatches = []
    for relative, expected in data["source_sha256"].items():
        source = original / relative
        actual = hashlib.sha256(source.read_bytes()).hexdigest() if source.is_file() else None
        if actual != expected:
            mismatches.append({"source": relative, "expected": expected, "actual": actual})
    return {"count": len(data["source_sha256"]), "mismatches": mismatches,
            "passed": not mismatches}

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--original", type=Path, default=DEFAULT_ORIGINAL)
    parser.add_argument("--build", type=Path, default=DEFAULT_BUILD)
    parser.add_argument("--output", type=Path, default=PROJECT / "output" / "task07_acceptance.json")
    args = parser.parse_args()
    original = args.original.resolve(); build = args.build.resolve(); output = args.output.resolve()
    if not original.is_dir() or not build.is_dir():
        raise SystemExit("original project and Release build directory are required")
    output.parent.mkdir(parents=True, exist_ok=True)
    temp_root = PROJECT / "tmp" / "task07_scene_acceptance"
    if temp_root.exists(): shutil.rmtree(temp_root)
    temp_root.mkdir(parents=True)
    report = {"task": "07", "date": "2026-11-05", "scope": "C++ migration acceptance",
              "original_project": str(original), "build": str(build),
              "python_validation": check_python(original),
              "cpp_release_ctest": check_ctest(build),
              "scene_cli": check_scene_runs(build, temp_root),
              "task05b_fixture_sources": check_fixture_sources(original),
              "limitations": ["arbitrary topology not migrated", "C01-C08 remain planned",
                              "GPU/CUDA not executed", "cloud execution not executed",
                              "Linux/MSVC/sanitizer/Docker not executed"]}
    report["passed"] = all(report[key]["passed"] for key in
                            ("python_validation", "cpp_release_ctest", "scene_cli",
                             "task05b_fixture_sources"))
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"passed": report["passed"],
                      "python_tests": report["python_validation"]["passed_tests"],
                      "cpp_tests": report["cpp_release_ctest"]["tests"],
                      "scene_paths": len(report["scene_cli"]["runs"]),
                      "fixture_sources": report["task05b_fixture_sources"]["count"],
                      "output": str(output)}, ensure_ascii=False, indent=2))
    return 0 if report["passed"] else 1

if __name__ == "__main__":
    raise SystemExit(main())
