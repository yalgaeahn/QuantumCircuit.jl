# QuantumCircuit.jl Roadmap

## Purpose

This roadmap translates the project plan into implementation order.
It is the source of truth for sequencing, dependencies, and milestone-sized task
groups.

## Current Status

Status as of 2026-03-19:
- Phase 1: complete
- Phase 2: complete
- Phase 2.5: complete
- Phase 3A: complete

## Phase 1: Static Core

Goal:
- deliver a stable fixed-parameter modeling workflow for transmon-resonator systems

Primary scope:
- `Transmon`
- `Resonator`
- `CapacitiveCoupling`
- generic `CompositeSystem`
- model-layer `hamiltonian(sys)`
- spectrum / transition / anharmonicity workflows

Task sequence:

1. Package bootstrap
- create the module skeleton under `src/`
- create the test layout under `test/`
- add formatting and CI basics if needed

2. Architecture layer
- implement `Transmon`
- implement `Resonator`
- implement `CapacitiveCoupling`
- implement `CompositeSystem`
- validate names, dimensions, and coupling references

3. Generic composition contract
- ensure `CompositeSystem` can hold future tunable subsystems without API redesign
- avoid fixed-only assumptions in internal storage
- validate couplings by symbolic subsystem name

4. Model layer for static systems
- define truncation rules for transmon and resonator subsystems
- implement local operator builders
- implement capacitive interaction terms
- implement a Duffing-approximate transmon Hamiltonian from `EJ` and `EC`
- implement `build_model(sys)`
- implement `hamiltonian(sys)` returning a QuantumToolbox-compatible object

5. Simulation layer for static systems
- implement `spectrum(sys; levels=...)`
- define a typed `SpectrumResult`
- route all Hamiltonian access through the Model layer

6. Analysis layer for static systems
- implement transition frequency extraction
- implement anharmonicity calculation
- derive analysis outputs from typed results

7. Tests and docs
- add regression tests for a transmon-only example
- add regression tests for a transmon-resonator example
- add one short quickstart example

Dependencies:
- model work depends on stable architecture types
- simulation depends on Model-layer Hamiltonian outputs
- analysis depends on typed simulation results

Exit criteria:
- a user can define a transmon-resonator `CompositeSystem`
- `hamiltonian(sys)` returns a QuantumToolbox-compatible Hamiltonian
- low-lying spectrum is computed through a stable API
- transition frequencies and anharmonicity are covered by tests

## Phase 2: Tunable Devices

Goal:
- extend the static core to flux-tunable systems without breaking the Phase-1 API surface

Primary scope:
- `TunableTransmon`
- `TunableCoupler` as a SQUID-based subsystem
- flux-dependent Hamiltonian regeneration
- parameter sweeps for tunable systems
- mixed fixed/tunable `CompositeSystem`

Task sequence:

1. Tunable architecture types
- implement `TunableTransmon`
- implement `TunableCoupler`
- define explicit flux and asymmetry-related parameters
- validate tunable parameters at construction time

2. Mixed-system composition
- confirm fixed and tunable subsystems coexist in one `CompositeSystem`
- ensure couplings reference tunable subsystems exactly like fixed ones
- avoid any tunable-only container

3. Flux-dependent model layer
- implement effective parameter updates from bias values
- extend `build_model(sys)` to bias-dependent systems
- extend `hamiltonian(sys)` to regenerate QuantumToolbox-compatible Hamiltonians under flux changes

4. Sweep support
- define sweep input specs
- implement parameter sweep execution for flux and coupling parameters
- define a typed `SweepResult`

5. Tunable-system analysis
- support basic tuning-curve summaries derived from typed sweep results

6. Tests and docs
- add regression tests for `TunableTransmon`
- add regression tests for `TunableCoupler`
- add mixed-system examples inside one `CompositeSystem`

Dependencies:
- requires Phase-1 `CompositeSystem` to be generic already
- tunable simulation depends on bias-dependent Hamiltonian generation
- tunable analysis depends on typed sweep or spectrum results

Exit criteria:
- tunable and non-tunable subsystems coexist inside one `CompositeSystem`
- flux-dependent Hamiltonians regenerate through stable APIs
- parameter sweeps produce structured results
- representative tunable workflows are covered by tests

## Phase 2.5: Stabilization and Analysis

Goal:
- stabilize the existing public API and complete the analysis layer for spectrum and sweep workflows

Primary scope:
- public API cleanup for `CompositeSystem`, `SweepSpec`, `spectrum`, and `simulate_sweep`
- typed analysis helpers for `SpectrumResult` and `SweepResult`
- short fixed, tunable, and mixed-system quickstarts
- explicit documentation of deferred features and approximation limits

Task sequence:

1. API surface review
- confirm the exported static and tunable APIs that should remain stable
- document deferred semantics such as the current status of `Transmon.ng`
- clean up error messages for invalid sweep targets and unsupported parameter combinations

2. Analysis layer expansion
- split `Analysis` into spectrum, sweep, and summary-oriented entry points
- add tuning-curve helpers derived purely from `SweepResult`
- add summary helpers that are plotting- and table-friendly
- keep advanced branch-tracking work out of this phase

3. Documentation
- add a fixed transmon quickstart
- add a tunable sweep quickstart
- add a mixed fixed/tunable system quickstart
- state the current Duffing-based modeling assumptions clearly

4. Regression coverage
- add tests for new analysis helpers
- add tests for documented sweep workflows
- keep Phase-1 and Phase-2 workflows green while refactoring analysis

Dependencies:
- requires stable Phase-1 and Phase-2 result types
- analysis work must not introduce new solver-layer dependencies

Exit criteria:
- documented examples cover fixed, tunable, and mixed-system workflows
- `SpectrumResult` and `SweepResult` have stable summary-style analysis helpers
- users can extract sweep curves without touching raw eigenspectra

## Phase 3A: Closed-System Dynamics

Goal:
- add time-domain closed-system simulation on top of the stabilized device models

Primary scope:
- named product-state preparation for composite systems
- reusable observable specifications and expectation traces
- static and time-dependent closed-system evolution via `sesolve`
- typed dynamics result objects
- subsystem-local driven evolution
- supported uncoupled circuit-mode closed-system workflows
- time-dependent flux control for supported circuit-mode and single-device nonadiabatic effective workflows
- one baseline dynamics walkthrough and one dedicated circuit-mode walkthrough

Task sequence:

1. Dynamics-ready model utilities
- expose model-layer helpers needed to build named product states
- expose model-layer helpers needed to build subsystem-local observables
- add a narrow drive representation for subsystem-local closed-system drives
- keep all Hamiltonian assembly in the Model layer

2. Simulation layer
- implement `evolve(system, ψ0, tlist; kwargs...)`
- implement overloads for model objects and direct Hamiltonian handoff where appropriate
- delegate static or driven Hamiltonian construction to the Model layer
- wrap `QuantumToolbox.sesolve` in a typed `DynamicsResult`

3. Result and analysis integration
- normalize solver expectation outputs into label-addressable traces
- add basic helpers such as `observable_trace(result, label)` and `final_state(result)`
- keep analysis free of direct solver calls

4. Tests and docs
- add regression tests for free evolution of a static system
- add regression tests for driven, circuit-mode, and flux-control closed-system examples
- add baseline and circuit-mode notebook walkthroughs for the documented APIs
- keep unsupported open-system and coupled circuit-mode features explicit

Dependencies:
- requires stable Phase-1/2 models and completed Phase-2.5 analysis boundaries
- relies on `QuantumToolbox.sesolve` as the only solver backend in this milestone

Exit criteria:
- `evolve(...)` runs closed-system dynamics through stable public APIs
- named initial-state and observable helpers cover the documented examples
- dynamics results are typed, documented, and test-covered
- supported uncoupled circuit-mode and single-device flux-control workflows are regression-verified
- the baseline and circuit-mode Phase 3A notebooks pass smoke verification

See also:
- `PHASE3A.md` for the implementation-level API and type design

## Phase 3B: Open-System Dynamics

Goal:
- add dissipative and steady-state workflows on top of the Phase-3A dynamics foundation

Primary scope:
- Lindblad dynamics
- steady-state workflows where useful
- collapse-operator construction
- typed open-system result objects aligned with `DynamicsResult`

Task sequence:

1. Dynamics-ready model support
- define collapse-operator construction where physically supported
- reuse the Phase-3A drive and observable conventions where possible

2. Simulation layer
- implement Lindblad evolution
- add steady-state support if it fits QuantumToolbox cleanly
- extend the typed dynamics result story without creating a parallel API

3. Result and analysis integration
- align dynamics results with spectrum and sweep result conventions
- add basic summaries for populations and observables
- keep analysis free of direct solver calls

4. Tests and docs
- add at least one documented open-system example
- add at least one documented driven-evolution example
- add regression tests for dynamics APIs

Dependencies:
- requires stable Phase-3A interfaces
- open-system workflows must reuse existing system and model entry points

Exit criteria:
- time evolution runs through stable public APIs
- Lindblad workflows are documented and tested
- dynamics results follow the same typed-result conventions as spectrum and sweep workflows

## Phase 4: Device Analysis and Calibration

Goal:
- extend analysis from basic summaries to calibration-oriented device characterization

Primary scope:
- avoided-crossing helpers built on explicit state-tracking support
- dispersive-shift extraction where the model fidelity supports it
- effective interaction summaries for coupled systems
- notebook-friendly summary adapters

Task sequence:

1. State-tracking prerequisites
- enrich result types with the metadata needed for state-overlap-based tracking
- keep tracking logic in analysis rather than simulation

2. Advanced analysis helpers
- implement avoided-crossing summaries
- implement dispersive-shift helpers
- implement higher-level coupled-system summaries

3. Documentation and coverage
- add calibration-style notebook examples
- add regression tests for tracked sweep analysis

Dependencies:
- requires stable sweep analysis from Phase 2.5
- may require richer eigenstate outputs from Simulation

Exit criteria:
- advanced analysis stays result-driven and test-covered
- calibration-oriented summaries are available without solver reimplementation

## Phase 5: Model Fidelity and Extensibility

Goal:
- raise model fidelity and add subsystem families without disturbing the public architecture

Primary scope:
- exact circuit-mode capacitive coupling through `CircuitCapacitiveCoupling(...; G = ...)`
- supported coupled exact circuit-mode local-drive and `FluxControl` workflows
- higher-fidelity transmon variants
- additional coupling models
- performance work for larger composite Hilbert spaces
- future subsystem families such as `Fluxonium`

Task sequence:

1. Fidelity upgrades
- add optional higher-fidelity subsystem models behind stable APIs
- preserve the current Phase-1/2 workflows as the default path where appropriate

2. Exact circuit coupling extension
- add `CircuitCapacitiveCoupling(...; G = ...)` as the exact circuit-mode capacitive interaction declaration
- keep `CapacitiveCoupling(...; g = ...)` as the effective-mode exchange declaration
- support transmon-like↔transmon-like and resonator↔transmon-like pair families under `CircuitHamiltonianSpec`
- keep resonator↔resonator exact circuit coupling explicitly unsupported until a later milestone

3. Coupled circuit workflows
- support `:G` sweeps through the existing sweep path
- support coupled exact circuit-mode local-drive dynamics and coupled `FluxControl` on the supported pair families
- keep paper-derived effective time-dependent couplings such as `g(t)` and `ḡ(t)` out of scope for this milestone

4. Scaling and validation
- improve performance for larger systems where measurements justify it
- add reference-style tests for exact circuit couplings, supported pair families, and coupled exact circuit flux behavior
- add comparison tests between approximate and higher-fidelity models when possible

Dependencies:
- requires the Phase-2.5 analysis boundary to remain clean
- should build on the same `CompositeSystem -> Model -> Simulation -> Analysis` flow

Exit criteria:
- exact circuit-mode capacitive coupling and supported coupled circuit workflows fit the existing architecture cleanly
- new subsystem or interaction families fit the existing architecture cleanly
- fidelity upgrades remain compatible with the package's public workflow

## Deferred Work

Not part of the current roadmap:
- hardware execution backends
- pulse scheduling or waveform compiler infrastructure
- automatic circuit quantization from schematic descriptions
- advanced calibration workflows
- direct hardware calibration loops

## Agent Checklist

Before starting a task:
- verify the target phase
- verify lower-layer dependencies are already implemented
- prefer extending existing types over introducing parallel abstractions
- add or update tests in the same change as the feature
