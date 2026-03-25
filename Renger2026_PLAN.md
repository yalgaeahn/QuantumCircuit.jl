# Renger 2026 Plan

## Local Parameter Inference Scope

This file tracks the paper-aligned local subsystem parameter inference task for the reduced device
`QB1-TC1-CR-TC2-QB2`.

- Use notebook-local code only.
- Do not change exported constructors or public API.
- Back out `QB1` and `QB2` bare local priors from the paper dressed qubit frequencies and the isolated local qubit anharmonicity target.
- Infer parked local priors for `TC1`/`TC2`.
- Treat coupler fits as paper-aligned local priors for later coupled-model fitting.

## Paper Targets

Primary targets from Table II and local bare/dressed interpretation:

- `QB1` dressed paper reference: `f = 4.67 GHz`, `alpha = -0.187 GHz`, `EJ/h = 14.8 GHz`, `EJ/EC = 74.2`
- `QB2` dressed paper reference: `f = 4.47 GHz`, `alpha = -0.187 GHz`, `EJ/h = 13.8 GHz`, `EJ/EC = 68.8`
- `CR` dressed paper reference: `f = 4.22 GHz`
- `CR` bare Hamiltonian input used in this repo: `f = 4.3 GHz`

Paper Table II qubit charging energies retained as metadata:

- `QB1`: `EC = 14.8 / 74.2 = 0.19946091644204852 GHz`
- `QB2`: `EC = 13.8 / 68.8 = 0.20058139534883723 GHz`

Secondary constraints retained for later joint fitting:

- `beta_qc(QB1) = 0.017`
- `beta_qc(QB2) = 0.022`
- `beta_cr = 0.0225`
- `beta_qr = 0.0020`
- Coupler parked isolated `f01 = 6.5 GHz`
- Coupler anharmonicity prior: `alpha_c ~= -0.11 GHz`
- Coupler parked flux is fixed at `0.0`
- Coupler asymmetry remains prior-held at `0.10`
- Explicit effective `g_*` seeds are frozen into the snapshot as provisional coupled-fit initialization values.

If Appendix C conflicts with Table II, Table II wins for this local-parameter task.

## Fit Rules

The notebook `output/jupyter-notebook/renger2026-parameter-inference.ipynb` uses:

- Exact single-device spectra from `CircuitHamiltonianSpec`
- A compact parked reduced effective-model back-out for qubits
- A parked `EC/EJmax` exact refit for couplers with `flux = 0.0` and `asymmetry = 0.10` held fixed
- Strong matching on dressed qubit references `f = 4.67/4.47 GHz` in the parked reduced effective model
- Strong matching on isolated local qubit anharmonicity `alpha = -0.187 GHz`
- Strong matching on coupler parked isolated `f01 = 6.5 GHz`
- Strong matching on coupler anharmonicity `alpha_c ~= -0.11 GHz`

Interpretation rules:

- `QB1` and `QB2` are treated as dressed-anchor bare back-out priors.
- Table II qubit `EJ`, `EJ/EC`, and derived `EC` remain paper metadata, not hard device constraints.
- Coupler fits remain weakly identified because local single-device targets still do not uniquely determine their full SQUID parameters, even though parked `f01` and `alpha` are now explicit anchors.
- `TC1` and `TC2` currently share the same accepted local prior because this task does not use branch-specific coupled information.

## Accepted Local Prior Set

The frozen snapshot is stored in `output/renger2026/paper_local_priors.toml` as `paper_local_priors_v1`.

Accepted values:

- `QB1`: `EJmax = 17.21`, `EC = 0.1705`, `flux = 0.0`, `asymmetry = 0.10`
- `QB2`: `EJmax = 15.88`, `EC = 0.17`, `flux = 0.0`, `asymmetry = 0.10`
- `TC1`: `EJmax = 51.02`, `EC = 0.107`, `flux = 0.0`, `asymmetry = 0.10`
- `TC2`: `EJmax = 51.02`, `EC = 0.107`, `flux = 0.0`, `asymmetry = 0.10`

Acceptance notes:

- The qubit device `f01_ghz` values are the isolated local bare frequencies derived from the current exact model.
- The qubit dressed paper references are stored separately in snapshot targets and are matched by the parked reduced effective model.
- The qubit device `anharmonicity_ghz` values are matched directly to the local isolated target `-0.187 GHz`.
- The qubit device `EJ/EC` values no longer equal the Table II entries; those remain snapshot metadata only.
- All accepted parking points sit at the flux sweet spot in the local model.
- The coupler parameters reproduce parked isolated `f01 ~= 6.5 GHz`.
- The coupler parameters reproduce `alpha_c ~= -0.11 GHz`.
- Coupler asymmetry remains prior-held at `0.10`; only `EC/EJmax` are newly anchored by the parked local solve.
- Weak isolated coupler `f01` values must not be reused to regenerate downstream effective couplings.
