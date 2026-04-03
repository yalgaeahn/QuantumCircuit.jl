#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
RETUNE_NOTEBOOK = REPO_ROOT / "output" / "jupyter-notebook" / "renger2026-figure2-ef-retune.ipynb"
REPRO_NOTEBOOK = REPO_ROOT / "output" / "jupyter-notebook" / "renger2026-figure2-cdef-reproduction.ipynb"
SHORTLIST_MARKER = "\nshortlist = joint_shortlist(baseline_snapshot, working_snapshot, working_gate_times)\n"


def extract_code(notebook_path: Path) -> str:
    notebook = json.loads(notebook_path.read_text())
    chunks: list[str] = []
    for cell in notebook["cells"]:
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))
        if source.strip():
            chunks.append(source)
    return "\n\n".join(chunks)


def notebook_head(notebook_path: Path) -> str:
    code = extract_code(notebook_path)
    head, marker, _ = code.partition(SHORTLIST_MARKER)
    if not marker:
        raise SystemExit(f"Could not find retune shortlist marker in {notebook_path}")
    return head


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


def build_refine_current_runner(
    tmpdir: Path,
    *,
    instantiate: bool,
    coarse_q_points: int,
    coarse_tc_points: int,
    coarse_rounds: int,
    full_q_points: int,
    full_tc_points: int,
    full_rounds: int,
    write_overlay: bool,
) -> Path:
    code = notebook_head(RETUNE_NOTEBOOK)
    replacements = {
        "Pkg.instantiate()": "Pkg.instantiate()" if instantiate else "nothing",
        "const RETUNE_COARSE_Q_POINTS = 9": f"const RETUNE_COARSE_Q_POINTS = {coarse_q_points}",
        "const RETUNE_COARSE_TC_POINTS = 9": f"const RETUNE_COARSE_TC_POINTS = {coarse_tc_points}",
        "const RETUNE_MAX_ROUNDS = 2": f"const RETUNE_MAX_ROUNDS = {coarse_rounds}",
        "const WRITE_WORKING_OVERLAY = false": f"const WRITE_WORKING_OVERLAY = {'true' if write_overlay else 'false'}",
    }
    for old, new in replacements.items():
        code = code.replace(old, new)

    code += f"""

function evaluate_full_candidate(base_snapshot, params)
    println("full-eval:move")
    move = branch_score_summary(
        base_snapshot,
        :move,
        merge(copy_branch_params(params["move"]), Dict("beta_cr" => params["beta_cr"]));
        charge_cutoff = RETUNE_CHARGE_CUTOFF,
        q_points = {full_q_points},
        tc_points = {full_tc_points},
        max_rounds = {full_rounds},
    )
    println("full-eval:cz")
    cz = branch_score_summary(
        base_snapshot,
        :cz,
        merge(copy_branch_params(params["cz"]), Dict("beta_cr" => params["beta_cr"]));
        charge_cutoff = RETUNE_CHARGE_CUTOFF,
        q_points = {full_q_points},
        tc_points = {full_tc_points},
        max_rounds = {full_rounds},
    )
    summary = joint_loss_summary(move, cz)
    return (loss = summary.loss, params = params, move = move, cz = cz, joint = summary)
end

println("retune:current-overlay")
current_candidate = evaluate_joint_candidate(
    baseline_snapshot,
    candidate_params_from_snapshot(working_snapshot),
)

println("retune:refine-current")
refined_candidate = refine_joint_candidate(baseline_snapshot, current_candidate)
if RUN_EC_FALLBACK || (
    refined_candidate.move.best_point.q_edge != :interior ||
    refined_candidate.move.best_point.tc_edge != :interior ||
    refined_candidate.cz.best_point.q_edge != :interior ||
    refined_candidate.cz.best_point.tc_edge != :interior
)
    println("retune:ec-fallback")
    refined_candidate = fallback_ec_refine(baseline_snapshot, refined_candidate)
end

println("retune:verify-current")
current_full = evaluate_full_candidate(baseline_snapshot, current_candidate.params)
println("retune:verify-refined")
refined_full = evaluate_full_candidate(baseline_snapshot, refined_candidate.params)

improved = refined_full.loss + 1e-10 < current_full.loss
refined_overlay = overlay_dict_from_candidate(baseline_snapshot, refined_candidate)
refined_overlay_toml = sprint(io -> TOML.print(io, refined_overlay))
written_overlay = if WRITE_WORKING_OVERLAY && improved
    write_overlay(RETUNE_OVERLAY_PATH, refined_overlay)
else
    nothing
end

println("RETUNE_RESULT_BEGIN")
show(stdout, MIME("text/plain"), (
    current_fast_loss = current_candidate.loss,
    refined_fast_loss = refined_candidate.loss,
    current_full_loss = current_full.loss,
    refined_full_loss = refined_full.loss,
    improved = improved,
    written_overlay = written_overlay,
    current_params = current_full.params,
    refined_params = refined_full.params,
    current_move_best = compact_namedtuple(current_full.move.best_point),
    current_cz_best = compact_namedtuple(current_full.cz.best_point),
    refined_move_best = compact_namedtuple(refined_full.move.best_point),
    refined_cz_best = compact_namedtuple(refined_full.cz.best_point),
))
println()
println("RETUNE_OVERLAY_TOML_BEGIN")
print(refined_overlay_toml)
println()
println("RETUNE_OVERLAY_TOML_END")
println("RETUNE_RESULT_END")
"""

    script_path = tmpdir / "retune_refine_current.jl"
    script_path.write_text(code)
    return script_path


def parse_command(_: argparse.Namespace) -> int:
    with tempfile.TemporaryDirectory(prefix="renger2026-retune-parse-") as tmp:
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
        result = run_julia_expr(parse_expr)
        ensure_ok(result, "Notebook parse check")
        print(result.stdout, end="")
        return 0


def refine_current_command(args: argparse.Namespace) -> int:
    with tempfile.TemporaryDirectory(prefix="renger2026-retune-refine-") as tmp:
        tmpdir = Path(tmp)
        runner = build_refine_current_runner(
            tmpdir,
            instantiate=args.instantiate,
            coarse_q_points=args.coarse_q_points,
            coarse_tc_points=args.coarse_tc_points,
            coarse_rounds=args.coarse_rounds,
            full_q_points=args.full_q_points,
            full_tc_points=args.full_tc_points,
            full_rounds=args.full_rounds,
            write_overlay=args.write_overlay,
        )
        result = run_julia(runner)
        ensure_ok(result, "Refine current overlay")
        print(result.stdout, end="")
        return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the Renger 2026 Fig. 2 overlay retune notebook headlessly.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    parse_parser = subparsers.add_parser("parse", help="Parse-check the reproduction and retune notebooks.")
    parse_parser.set_defaults(func=parse_command)

    refine_parser = subparsers.add_parser(
        "refine-current",
        help="Refine the current working overlay at cc=10 and optionally write it if it improves.",
    )
    refine_parser.add_argument("--instantiate", action="store_true", help="Run Pkg.instantiate() before executing Julia.")
    refine_parser.add_argument("--write-overlay", action="store_true", help="Write the overlay file if the refined candidate improves.")
    refine_parser.add_argument("--coarse-q-points", type=int, default=5, help="Reduced q-flux grid size used during refinement.")
    refine_parser.add_argument("--coarse-tc-points", type=int, default=5, help="Reduced TC-flux grid size used during refinement.")
    refine_parser.add_argument("--coarse-rounds", type=int, default=1, help="Adaptive scan rounds used during refinement.")
    refine_parser.add_argument("--full-q-points", type=int, default=9, help="Verification q-flux grid size.")
    refine_parser.add_argument("--full-tc-points", type=int, default=9, help="Verification TC-flux grid size.")
    refine_parser.add_argument("--full-rounds", type=int, default=2, help="Verification adaptive scan rounds.")
    refine_parser.set_defaults(func=refine_current_command)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
