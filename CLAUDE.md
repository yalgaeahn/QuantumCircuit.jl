# QuantumCircuit.jl

Circuit-QED Hamiltonian simulation software for reproducing published results.
Primary goal: reproduce Renger et al. (2026) Fig. 2 (MOVE/CZ gate dynamics on a QB1-TC1-CR-TC2-QB2 device).

## How to run tests

```bash
julia --project=. -e "using Pkg; Pkg.test()"
```

All 320 tests should pass.

## Key files

- `src/Analysis/Renger2026.jl` — model builders, snapshot loading, `renger2026_model_pair()`
- `output/renger2026/paper_local_priors.toml` — baseline device parameters (QB1/QB2 identified, TC1/TC2 weakly identified)
- `output/renger2026/fig2_ef_retune_working.toml` — retune overlay with corrected TC params (must be applied for Fig. 2)
- `output/jupyter-notebook/renger2026_fig2_branch_helpers.jl` — 1518-line helper for flux sweeps, branch specs, gate dynamics

## Architecture

Dependency flow: `Architecture -> Model -> Simulation -> Analysis`

- **Architecture**: Subsystem types (`Transmon`, `TunableTransmon`, `TunableCoupler`, `Resonator`), coupling types (`CapacitiveCoupling`, `CircuitCapacitiveCoupling`), `CompositeSystem`
- **Model**: Two Hamiltonian specs:
  - `EffectiveHamiltonianSpec` — Duffing approximation (fast, works for 5-body)
  - `CircuitHamiltonianSpec(charge_cutoff=N)` — exact charge-basis (accurate, expensive)
- **Simulation**: `spectrum()`, `simulate_sweep()`, `evolve()`
- **Analysis**: spectral analysis, sweep processing, `Renger2026.jl` helpers

---

# Diagnosis: Renger 2026 Fig. 2 Reproduction Failure

## Summary

The Fig. 2 reproduction failed due to **charge_cutoff being too small** in the circuit Hamiltonian. This caused all subsystem energy levels to be wrong by hundreds of MHz to several GHz, making every downstream result (avoided crossings, gate dynamics) physically incorrect.

## Root Cause 1 (Primary): the old notebook-based Fig. 2 workflow used charge_cutoff=3 and gave completely wrong energy levels

`CircuitHamiltonianSpec(charge_cutoff=N)` truncates the charge basis to `2N+1` states per transmon-like subsystem. For deep transmons (high EJ/EC), many charge states are needed to resolve the cosine potential.

The removed notebook-based Fig. 2 workflow used `charge_cutoff=3`, which was catastrophically inadequate:

### Convergence data (isolated subsystem, flux=0)

| Subsystem | EJmax (GHz) | EC (GHz) | EJ/EC | cc=3 | cc=5 | cc=7 | cc=10 | cc=13 |
|---|---|---|---|---|---|---|---|---|
| QB1 | 17.21 | 0.1705 | ~101 | 5.32 GHz | 4.679 GHz | 4.668 GHz | **4.668 GHz** | 4.668 GHz |
| QB2 | 15.88 | 0.17 | ~93 | — | — | 4.470 GHz | **4.470 GHz** | 4.470 GHz |
| TC1 retune | 35.75 | 0.105 | ~340 | 8.76 GHz | 5.77 GHz | 5.389 GHz | **5.373 GHz** | 5.373 GHz |
| TC2 retune | 32.0 | 0.105 | ~305 | 7.94 GHz | 5.38 GHz | 5.087 GHz | **5.077 GHz** | 5.077 GHz |
| TC1 baseline | 51.02 | 0.107 | ~477 | 12.09 GHz | 7.38 GHz | 6.56 GHz | **6.500 GHz** | 6.500 GHz |

**Bold** = converged (matches cc=13 to < 0.1 MHz).

Key observations:
- **cc=3**: QB1 at 5.32 GHz instead of 4.67 GHz (+650 MHz error). TC1 retune at 8.76 GHz instead of 5.37 GHz (+3385 MHz). Completely wrong physics.
- **cc=7**: QB converged, but TC1 retune still +16 MHz off, TC2 retune +10 MHz off. Not converged.
- **cc=10**: All subsystems converged. Matches cc=13 to < 0.1 MHz.

### 3-body Hilbert space sizes (QB + TC + CR, resonator dim=3)

| charge_cutoff | Local states per transmon | 3-body Hilbert space |
|---|---|---|
| 3 | 7 | 147 states |
| 5 | 11 | 363 states |
| 7 | 15 | 675 states |
| **10** | **21** | **1323 states** |

cc=10 with 1323 states is fast to diagonalize. There is no reason to use a smaller cutoff.

## Root Cause 2: TC1/TC2 baseline EJmax values are wrong (~51 GHz → ~32-36 GHz)

The `paper_local_priors.toml` has TC1/TC2 with EJmax=51.02 GHz, but this value is **wrong**. It was fitted only from the parked coupler frequency (f01=6.5 GHz) and anharmonicity (alpha=-0.11 GHz), which is underdetermined — many (EJmax, EC, asymmetry) triples produce the same parked f01. The fitting chose an EJmax ~50% too high.

### Baseline vs Retune parameter comparison

| Parameter | TC1 baseline | TC1 retune | TC2 baseline | TC2 retune |
|---|---|---|---|---|
| EJmax (GHz) | 51.02 | **35.75** | 51.02 | **32.0** |
| EC (GHz) | 0.107 | **0.105** | 0.107 | **0.105** |
| asymmetry | 0.0 | **0.14** | 0.0 | **0.04** |
| EJ/EC | ~477 | ~340 | ~477 | ~305 |

The retune overlay (`fig2_ef_retune_working.toml`) reduces TC EJmax by 30-37%, which fundamentally changes the coupler's flux-tuning curve and avoided crossing positions. These retune values came from the removed notebook-based retune workflow.

**Critical caveat**: The retune fitting was performed with charge_cutoff=3, which gave wrong energy levels (see Root Cause 1). The retune parameters may therefore need re-fitting with charge_cutoff=10 for self-consistent results.

**The retune overlay must always be applied for Fig. 2 reproduction.** The baseline TC parameters are placeholders, not physical truth.

## Root Cause 3: Effective vs Circuit Hamiltonian inconsistency

The effective model (Duffing approximation with `CapacitiveCoupling(g=...)`) and the circuit model (`CircuitCapacitiveCoupling(G=beta)`) use different physical representations:
- Effective: g-parameters (Jaynes-Cummings style exchange coupling)
- Circuit: beta (capacitive ratios) with charge-charge interaction

With the baseline snapshot at charge_cutoff=2 (the default in `renger2026_model_pair`), the 5-body effective vs circuit spectra diverged by 2-6 GHz at levels 3+. This was entirely due to the inadequate charge_cutoff — not a fundamental model mismatch.

## Root Cause 4: 5-body circuit model is impractical

The full QB1-TC1-CR-TC2-QB2 system with charge_cutoff=10 would have 21^4 x 3 = 583,443 states — too large for direct diagonalization. The `renger2026_model_pair()` defaults to charge_cutoff=2 for structural testing only.

For accurate physics: use 3-body QCR branch models (QB+TC+CR, 1323 states at cc=10) with `CircuitHamiltonianSpec`, or use the 5-body model with `EffectiveHamiltonianSpec` (Duffing approximation, 1200 states).

## Root Cause 5: QB1/QB2 back-out uses approximate effective model

The qubit parameters (EJmax, EC) were inferred from paper-reported dressed frequencies using a "parked reduced effective model" (Duffing approximation). The residuals are small (<1 MHz), so this is not a major issue. However, the inference was done with charge_cutoff=10 for the circuit model validation, while the removed notebook workflow used charge_cutoff=3 — making that validation meaningless.

---

# Remaining work

1. **Rebuild the Fig. 2 reproduction stack end-to-end** with charge_cutoff=10 and verify:
   - Fig. 2(c)(d) avoided crossings appear at qualitatively correct flux positions
   - Fig. 2(e) MOVE gate shows population transfer QB1 -> CR -> QB1
   - Fig. 2(f) CZ gate shows entangling phase accumulation
2. **Consider adding charge_cutoff convergence table to ARCHITECTURE.md**
3. **TC1/TC2 retune may need further iteration** — current retune overlay was found with cc=3, which was wrong. The retune should be re-done with cc=10 for self-consistent results.

## Rules for future work

- **Never use charge_cutoff < 5** for any real device simulation
- **Use charge_cutoff=10** for Fig. 2 branch models (fully converged, 1323 states)
- **Always apply the retune overlay** (`fig2_ef_retune_working.toml`) for Fig. 2 work
- The `renger2026_model_pair()` default charge_cutoff=2 is for structural testing only — do not use it for physics
- For 5-body models, use `EffectiveHamiltonianSpec` (Duffing approximation)
