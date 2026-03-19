# Renger 2026 Plan

## Local Parameter Inference Scope

This file tracks the paper-aligned local subsystem parameter inference task for the reduced device
`QB1-TC1-CR-TC2-QB2`.

- Use notebook-local code only.
- Do not change exported constructors or public API.
- Infer `EC`, `EJmax`, `flux`, and `asymmetry` for `QB1`, `TC1`, `QB2`, and `TC2`.
- Treat qubit fits as the primary outputs.
- Treat coupler fits as paper-aligned local priors for later coupled-model fitting.

## Paper Targets

Primary targets from Table II:

- `QB1`: `f01 = 4.67 GHz`, `EJ/h = 14.8 GHz`, `EJ/EC = 74.2`
- `QB2`: `f01 = 4.47 GHz`, `EJ/h = 13.8 GHz`, `EJ/EC = 68.8`
- `CR`: `f = 4.22 GHz`

Derived qubit charging energies:

- `QB1`: `EC = 14.8 / 74.2 = 0.19946091644204852 GHz`
- `QB2`: `EC = 13.8 / 68.8 = 0.20058139534883723 GHz`

Secondary constraints retained for later joint fitting:

- `beta_qc(QB1) = 0.017`
- `beta_qc(QB2) = 0.022`
- `beta_cr = 0.0225`
- `beta_qr = 0.0020`
- Coupler anharmonicity prior: `alpha_c ~= -0.11 GHz`

If Appendix C conflicts with Table II, Table II wins for this local-parameter task.

## Fit Rules

The notebook `output/jupyter-notebook/renger2026-parameter-inference.ipynb` uses:

- Exact single-device spectra from `CircuitHamiltonianSpec`
- A deterministic coarse-to-fine search over `flux`, `asymmetry`, `EC`, and `EJmax`
- Strong sweet-spot priors at `flux ~= 0`
- Strong matching on qubit `f01` and derived `EC`
- Broad matching on qubit parked effective `EJ`
- Strong matching on coupler anharmonicity and only broad parked-frequency plausibility

Interpretation rules:

- Qubit fits are considered identified enough to freeze as local priors.
- Coupler fits are only weakly identified because local single-device targets do not uniquely determine their SQUID parameters.
- `TC1` and `TC2` currently share the same accepted local prior because this task does not use branch-specific coupled information.

## Accepted Local Prior Set

The frozen snapshot is stored in `output/renger2026/paper_local_priors.toml` as `paper_local_priors_v1`.

Accepted values:

- `QB1`: `EJmax = 14.92`, `EC = 0.19946091644204852`, `flux = 0.0`, `asymmetry = 0.10`
- `QB2`: `EJmax = 13.65`, `EC = 0.20058139534883723`, `flux = 0.0`, `asymmetry = 0.10`
- `TC1`: `EJmax = 34.25`, `EC = 0.105`, `flux = 0.0`, `asymmetry = 0.10`
- `TC2`: `EJmax = 34.25`, `EC = 0.105`, `flux = 0.0`, `asymmetry = 0.10`

Acceptance notes:

- The qubit fits reproduce the Table II `f01` targets within `20 MHz`.
- The qubit `EC` values match the derived paper values within `5 MHz`.
- All accepted parking points sit at the flux sweet spot in the local model.
- The coupler parameters reproduce `alpha_c ~= -0.11 GHz` but remain only weakly identified until coupled avoided-crossing data is used.
