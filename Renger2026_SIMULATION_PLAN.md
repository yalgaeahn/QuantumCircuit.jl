# Renger 2026 Simulation Workflow

This file tracks the implemented Hamiltonian-simulation workflow for reproducing the closed-system parts of
Renger et al. (2026) with QuantumCircuit.jl.

Primary parameter sources:

- [Renger2026_PLAN.md](/Users/yalgaeahn/Research/20_Projects/QuantumCircuit.jl/Renger2026_PLAN.md)
- [paper_local_priors.toml](/Users/yalgaeahn/Research/20_Projects/QuantumCircuit.jl/output/renger2026/paper_local_priors.toml)

## Current package support

- Package-level subspace and gate-analysis APIs:
  - `subspace_spec(...)`
  - `projected_unitary(...)`
  - `strip_local_z_phases(U)`
  - `conditional_phase(U)`
  - `best_move(...)`
  - `best_cz(...)`
- Effective-versus-full comparison helpers:
  - `compare_spectra(...)`
  - `compare_model_spectra(...)`
  - `compare_minimum_gap(...)`
- Paper-aligned system builders:
  - `load_renger2026_snapshot(...)`
  - `renger2026_stage1_qr_system(...)`
  - `renger2026_stage1_qcr_system(...)`
  - `renger2026_reduced_system(...)`
  - `renger2026_model_pair(...)`

## Staged workflow

### Stage 1 - Core MOVE sanity checks

- Build a simple `Q-R` model and verify single-excitation exchange.
- Build a `Q-C-R` chain and verify mediated transfer with explicit leakage accounting.
- Use `best_move(...)` on projected unitaries for calibration-style checks, and `population_trace(...)` for state-level sanity checks.

### Stage 2 - Reduced paper device

- Build the reduced `QB1-TC1-CR-TC2-QB2` chain from the frozen local priors.
- Use the effective model as a fast surrogate and the exact circuit model as the source-of-truth validation model.
- Treat the `beta_*` values in `paper_local_priors.toml` only as seeded coupling guesses.

### Stage 3 - Effective-versus-full comparison

- Compare low-lying spectra with `compare_model_spectra(...)`.
- Compare avoided-crossing locations with `compare_minimum_gap(...)`.
- Accept only comparisons that are backed by explicit truncation sweeps.

### Stage 4 - Pulse-level MOVE

- Use `SubsystemDrive(...)` and `FluxControl(...)` directly; no pulse compiler is introduced in this pass.
- Extract projected unitaries on the chosen logical subspace.
- Calibrate transfer times with `best_move(...)`.

### Stage 5 - Pulse-level CZ

- Use an exact circuit `Q-C-Q` or `QB1-TC1-CR-TC2-QB2` model with flux control.
- Extract logical two-qubit unitaries with `projected_unitary(...)`.
- Remove local phase freedom with `strip_local_z_phases(...)` when needed.
- Measure the entangling phase with `conditional_phase(...)`.
- Calibrate operating points with `best_cz(...)`.

## Parameter interpretation

- `QB1` and `QB2` device fields now store bare local back-out priors chosen so the parked reduced effective model matches the paper dressed qubit frequencies while the isolated local qubit anharmonicity stays near `-0.187 GHz`.
- `devices.QB*.f01_ghz` in the frozen snapshot are isolated local bare frequencies derived from the current exact model.
- `targets.qb*_dressed_f01_ghz` store the paper dressed references separately.
- `targets.qb*_alpha_ghz` store the qubit anharmonicity targets separately.
- Table II qubit `EJ`, `EJ/EC`, and derived `EC` remain paper metadata in `targets`, not hard device constraints.
- `TC1` and `TC2` priors remain weakly identified initialization points, but their parked isolated local targets are explicitly anchored at `f01 = 6.5 GHz` and `alpha = -0.11 GHz`.
- `TC1` and `TC2` keep `flux = 0.0` and `asymmetry = 0.10` in the baseline local-prior solve; only `EC/EJmax` are refit there.
- `CR` uses the bare resonator Hamiltonian frequency from `paper_local_priors.toml` (`4.3 GHz`), while the paper dressed reference (`4.22 GHz`) is stored separately for reporting.
- Effective couplings are seeded from explicit `g_*` values frozen in the snapshot.
- Circuit couplings use the raw `beta_*` values as provisional `G` seeds.
- These seeded couplings are not device truth and must be refit against coupled spectra and gate data.

## Practical compute limits

- Full stored paper truncations are too expensive for routine calibration loops.
- Default reduced workflow:
  - effective qubits: `ncut = 5`
  - effective couplers: `ncut = 4`
  - resonator: `dim = 3`
  - circuit charge cutoff: `2`
- Any accepted spectrum, MOVE, or CZ result should be rechecked with explicit truncation convergence sweeps before being treated as validated.

## Current regression targets

- `Q-R` projected MOVE reaches near-unit transfer in the reduced effective model.
- `Q-C-R` projected MOVE shows mediated transfer with nonzero leakage.
- Reduced `QB1-TC1-CR-TC2-QB2` effective dynamics show stronger resonator response near the seeded paper operating point than under a large resonator detuning shift.
- Exact circuit `Q-C-Q` flux control reaches conditional phase near `pi` with low leakage at a reduced operating point used in tests.
