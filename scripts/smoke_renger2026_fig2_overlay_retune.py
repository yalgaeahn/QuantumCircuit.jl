#!/usr/bin/env python3

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
REPRO_NOTEBOOK = REPO_ROOT / "output" / "jupyter-notebook" / "renger2026-figure2-cdef-reproduction.ipynb"
RETUNE_NOTEBOOK = REPO_ROOT / "output" / "jupyter-notebook" / "renger2026-figure2-ef-retune.ipynb"


def extract_code(notebook_path: Path) -> str:
    notebook = json.loads(notebook_path.read_text())
    chunks: list[str] = []
    for cell in notebook["cells"]:
        if cell.get("cell_type") == "code":
            source = "".join(cell.get("source", []))
            if source.strip():
                chunks.append(source)
    return "\n\n".join(chunks)


def run_julia(code_path: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["julia", f"--project={REPO_ROOT}", str(code_path)],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )


def run_julia_expr(expr: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["julia", f"--project={REPO_ROOT}", "-e", expr],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )


def ensure_ok(result: subprocess.CompletedProcess[str], label: str) -> None:
    if result.returncode == 0:
        return
    raise SystemExit(
        f"{label} failed with exit code {result.returncode}\n"
        f"STDOUT:\n{result.stdout}\n"
        f"STDERR:\n{result.stderr}"
    )


def build_retune_tiny_runner(tmpdir: Path) -> Path:
    notebook = json.loads(RETUNE_NOTEBOOK.read_text())
    code = "\n\n".join("".join(notebook["cells"][idx]["source"]) for idx in (2, 4, 6))
    replacements = {
        "Pkg.instantiate()": "nothing",
        "const RETUNE_PROFILE = :compact": "const RETUNE_PROFILE = :tiny",
        "const RETUNE_MAX_ROUNDS = 2": "const RETUNE_MAX_ROUNDS = 1",
        "const RETUNE_SHORTLIST_SIZE = 6": "const RETUNE_SHORTLIST_SIZE = 1",
        "const RETUNE_COARSE_Q_POINTS = 9": "const RETUNE_COARSE_Q_POINTS = 3",
        "const RETUNE_COARSE_TC_POINTS = 9": "const RETUNE_COARSE_TC_POINTS = 3",
        "const RETUNE_FINE_Q_POINTS = 21": "const RETUNE_FINE_Q_POINTS = 5",
        "const RETUNE_FINE_TC_POINTS = 21": "const RETUNE_FINE_TC_POINTS = 5",
        "const RUN_REFINEMENT = true": "const RUN_REFINEMENT = false",
    }
    for old, new in replacements.items():
        code = code.replace(old, new)

    code += """

shortlist = joint_shortlist(baseline_snapshot, working_snapshot, working_gate_times; profile = RETUNE_PROFILE, shortlist_size = 1)
@assert length(shortlist) == 1
candidate = shortlist[1]

candidate_snapshot = joint_candidate_snapshot(baseline_snapshot, candidate.params)
move_panel = render_branch_panel(candidate_snapshot, :move, candidate.params["move"]["t_gate"]; charge_cutoff = 2)
cz_panel = render_branch_panel(candidate_snapshot, :cz, candidate.params["cz"]["t_gate"]; charge_cutoff = 2)

fig = Figure(size = (900, 320))
ax_move = Axis(fig[1, 1]; title = "MOVE tiny smoke", xlabel = "omega_QB,01 (GHz)", ylabel = "omega_TC (GHz)")
hm_move = heatmap!(ax_move, move_panel.x, move_panel.y, move_panel.matrix; colorrange = (0.0, 1.0), colormap = :viridis)
xlims!(ax_move, move_panel.display.xlims...)
ylims!(ax_move, move_panel.display.ylims...)
scatter!(ax_move, [move_panel.best_point.q_ghz], [move_panel.best_point.tc_ghz]; color = :white, strokecolor = :black, markersize = 10)

ax_cz = Axis(fig[1, 2]; title = "CZ tiny smoke", xlabel = "omega_QB,01 (GHz)", ylabel = "omega_TC (GHz)")
hm_cz = heatmap!(ax_cz, cz_panel.x, cz_panel.y, cz_panel.matrix; colorrange = (0.0, 1.0), colormap = :viridis)
xlims!(ax_cz, cz_panel.display.xlims...)
ylims!(ax_cz, cz_panel.display.ylims...)
scatter!(ax_cz, [cz_panel.best_point.q_ghz], [cz_panel.best_point.tc_ghz]; color = :white, strokecolor = :black, markersize = 10)

Colorbar(fig[:, 3], hm_cz; label = "final qubit P_e")
tiny_gallery = save_figure(fig, repo_root, "figure2ef_retune_tiny_smoke")

overlay = overlay_dict_from_candidate(baseline_snapshot, candidate)
@assert Set(keys(overlay["devices"])) == Set(["TC1", "TC2"])
@assert Set(keys(overlay["targets"])) == Set(["beta_qc_qb1", "beta_qc_qb2", "beta_cr", "fig2_move_t_gate_ns", "fig2_cz_t_gate_ns"])

println("tiny-gallery=", tiny_gallery)
println("shortlist-loss=", round(candidate.loss; digits = 6))
println("move-best=", candidate.move.best_point)
println("cz-best=", candidate.cz.best_point)
println("overlay-target-keys=", sort(collect(keys(overlay["targets"]))))
"""
    script_path = tmpdir / "retune_tiny_smoke.jl"
    script_path.write_text(code)
    return script_path


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="renger2026-overlay-retune-") as tmp:
        tmpdir = Path(tmp)
        repro_code = extract_code(REPRO_NOTEBOOK)
        retune_code = extract_code(RETUNE_NOTEBOOK)

        repro_path = tmpdir / "repro_extracted.jl"
        retune_path = tmpdir / "retune_extracted.jl"
        repro_path.write_text(repro_code)
        retune_path.write_text(retune_code)

        parse_expr = (
            f'Meta.parseall(read("{repro_path}", String)); '
            f'Meta.parseall(read("{retune_path}", String)); '
            'println("notebook-parse-ok")'
        )
        ensure_ok(run_julia_expr(parse_expr), "Notebook parse check")

        tiny_runner = build_retune_tiny_runner(tmpdir)
        result = run_julia(tiny_runner)
        ensure_ok(result, "Tiny retune smoke")
        print(result.stdout, end="")


if __name__ == "__main__":
    main()
