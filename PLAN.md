# QuantumCircuit.jl Plan

## Goal

QuantumCircuit.jl is a superconducting quantum circuit Hamiltonian simulator
built on top of QuantumToolbox.jl.

The package should let users describe device-level circuit architectures and
compute quantum properties such as:

- energy spectrum
- transition frequencies
- anharmonicity
- dispersive shifts
- driven time evolution
- open-system dynamics

Primary device families:
- `Transmon`
- `TunableTransmon`
- `Resonator`
- `TunableCoupler`

Primary interaction family:
- `CapacitiveCoupling`

The focus is architecture-level modeling, not building a general-purpose solver
library from scratch.

## Document Roles

- `PLAN.md`: project scope, phase boundaries, milestone criteria
- `ARCHITECTURE.md`: layer ownership, domain model, API boundaries
- `ROADMAP.md`: implementation order and task sequence
- `AGENT.md`: entrypoint and rules for coding agents
- `PHASE3A.md`: detailed design for the first dynamics milestone

## Project Status

Status as of 2026-03-13:
- Phase 1: complete
- Phase 2: complete
- Phase 2.5: complete
- Phase 3A: in progress

## Phase Plan

### Phase 1: Static Core

Goal:
- establish a stable fixed-parameter workflow

In scope:
- `Transmon`
- `Resonator`
- `CapacitiveCoupling`
- generic `CompositeSystem`
- model-layer Hamiltonian construction
- Duffing-approximate transmon modeling from `EJ` and `EC`
- QuantumToolbox-compatible Hamiltonian outputs
- spectrum, transition-frequency, and anharmonicity workflows

Out of scope:
- `TunableTransmon`
- `TunableCoupler`
- Lindblad dynamics
- driven evolution
- hardware or pulse-control integration

### Phase 2: Tunable Devices

Goal:
- extend the same architecture to flux-tunable systems

In scope:
- `TunableTransmon`
- `TunableCoupler` as a SQUID-based independent subsystem
- mixed fixed/tunable `CompositeSystem`
- flux-dependent Hamiltonian regeneration
- parameter sweeps over tunable parameters

Out of scope:
- pulse scheduling
- calibration workflows
- automatic circuit quantization

### Phase 2.5: Stabilization and Analysis

Goal:
- stabilize the Phase-1/2 API surface and complete result-driven analysis for static and tunable sweeps

In scope:
- public API cleanup for `CompositeSystem`, `SweepSpec`, `spectrum`, and `simulate_sweep`
- analysis helpers derived from `SpectrumResult` and `SweepResult`
- tuning-curve and summary-style helpers for fixed, tunable, and mixed systems
- short user-facing docs for static, tunable, and mixed-system workflows
- validation and error-message cleanup
- explicit documentation of deferred model fidelity features such as `ng` support

Out of scope:
- new solver capabilities
- branch tracking that depends on eigenstate continuity
- avoided-crossing automation that requires state-overlap logic
- dynamics or open-system analysis
- pulse scheduling or hardware integration

### Phase 3A: Closed-System Dynamics

Goal:
- add time-domain closed-system simulation on top of the stabilized static and tunable models

In scope:
- Schrödinger time evolution through stable public APIs
- reusable typed result objects for time-domain simulations
- named initial-state helpers for composite systems
- named observable helpers and expectation-value traces
- narrowly scoped driven closed-system evolution built on model-layer Hamiltonian assembly

Out of scope:
- Lindblad dynamics
- collapse-operator construction
- steady-state workflows
- hardware execution backends
- waveform compiler infrastructure

### Phase 3B: Open-System Dynamics

Goal:
- extend Phase 3A dynamics to dissipative and steady-state workflows without changing the system model

In scope:
- Lindblad dynamics
- collapse-operator construction where physically supported
- steady-state workflows where useful
- typed result objects that stay aligned with the Phase-3A dynamics conventions

Out of scope:
- stochastic trajectory tooling beyond what is needed for the first Lindblad milestone
- hardware execution backends
- waveform compiler infrastructure

### Phase 4: Device Analysis and Calibration

Goal:
- extend the result-driven analysis layer toward device-characterization workflows

In scope:
- avoided-crossing helpers once state-overlap tracking exists
- dispersive-shift extraction where the underlying model supports it
- effective-coupling and calibration-oriented summaries
- table-friendly export adapters for notebook workflows

Out of scope:
- real hardware calibration loops
- pulse compiler infrastructure

### Phase 5: Model Fidelity and Extensibility

Goal:
- improve model fidelity and widen subsystem coverage without breaking the composition API

In scope:
- higher-fidelity transmon modeling beyond the initial Duffing approximation
- additional interaction families where physically justified
- performance work for larger Hilbert spaces
- future subsystem families such as `Fluxonium`

Out of scope:
- replacing the layer architecture
- introducing solver-specific parallel object models

## Milestone Criteria

Phase 1 is complete when:
- a user can define a transmon-resonator `CompositeSystem`
- the package can build a QuantumToolbox-compatible Hamiltonian
- low-lying spectrum, transition frequencies, and anharmonicity are available from stable APIs
- the main workflow is covered by tests and short docs

Phase 2 is complete when:
- a user can define `TunableTransmon` and `TunableCoupler` systems with explicit bias parameters
- tunable and non-tunable subsystems coexist inside one `CompositeSystem`
- flux sweeps regenerate Hamiltonians through stable APIs
- representative tunable workflows are covered by tests

Phase 2.5 is complete when:
- the public static and tunable APIs are stable enough to document without caveats
- `SpectrumResult` and `SweepResult` support standard summary-style analysis helpers
- the package has short docs for fixed, tunable, and mixed-system examples
- deferred or unsupported physics features are documented explicitly

Phase 3A is complete when:
- closed-system evolution runs through stable public APIs
- at least one documented dynamics example exists for a static or tunable system
- time-domain results follow the same typed-result conventions as spectrum and sweep workflows
- named initial-state and observable helpers cover the documented workflows

Phase 3B is complete when:
- closed- and open-system evolution run through stable public APIs
- at least one documented dynamics example exists for a tunable or coupled system
- dynamics results follow the same typed-result conventions as spectrum and sweep workflows

Phase 4 is complete when:
- advanced spectrum/sweep analysis is available through stable APIs
- at least one calibration-oriented analysis workflow is documented
- advanced analysis remains layered on typed results without rebuilding Hamiltonians

Phase 5 is complete when:
- higher-fidelity subsystem models fit the existing `CompositeSystem -> Model -> Simulation -> Analysis` flow
- at least one new subsystem or interaction family extends the existing architecture cleanly
- fidelity improvements preserve backward-compatible user workflows where feasible
