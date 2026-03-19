#!/usr/bin/env python3

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS = [
    REPO_ROOT / "output" / "jupyter-notebook" / "phase3a-closed-system-dynamics-walkthrough.ipynb",
    REPO_ROOT / "output" / "jupyter-notebook" / "phase3a-circuit-hamiltonian-dynamics.ipynb",
]


def extract_code_cells(notebook_path: Path) -> str:
    notebook = json.loads(notebook_path.read_text())
    code_cells = []

    for index, cell in enumerate(notebook.get("cells", [])):
        if cell.get("cell_type") != "code":
            continue

        source = "".join(cell.get("source", []))
        if not source.strip():
            continue

        code_cells.append(f"# --- {notebook_path.name} cell {index} ---\n{source.rstrip()}\n")

    if not code_cells:
        raise ValueError(f"{notebook_path} does not contain any executable code cells.")

    return "\n".join(code_cells) + "\n"


def run_notebook(notebook_path: Path) -> None:
    julia = os.environ.get("JULIA") or shutil.which("julia") or "julia"
    script_text = extract_code_cells(notebook_path)

    with tempfile.NamedTemporaryFile("w", suffix=".jl", delete=False) as handle:
        handle.write(script_text)
        script_path = Path(handle.name)

    command = [julia, "--project=.", str(script_path)]
    print(f"==> Smoke-running {notebook_path.name}")

    try:
        subprocess.run(command, cwd=REPO_ROOT, check=True)
    finally:
        if script_path.exists():
            script_path.unlink()


def main() -> int:
    missing = [path for path in NOTEBOOKS if not path.is_file()]
    if missing:
        for path in missing:
            print(f"Missing notebook: {path}", file=sys.stderr)
        return 1

    try:
        for notebook_path in NOTEBOOKS:
            run_notebook(notebook_path)
    except (subprocess.CalledProcessError, ValueError) as exc:
        print(exc, file=sys.stderr)
        return 1

    print("Phase 3A notebook smoke verification passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
