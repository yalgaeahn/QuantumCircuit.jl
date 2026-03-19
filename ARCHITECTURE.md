# QuantumCircuit.jl Architecture

## Purpose

This document defines the implementation architecture for QuantumCircuit.jl.
It is the source of truth for layer ownership, domain-model boundaries, and
extension rules.

## Core Invariants

- dependency flow is one-way: `Architecture -> Model -> Simulation -> Analysis`
- `CompositeSystem` is the single system-composition entry point
- fixed and tunable subsystems must coexist inside the same `CompositeSystem`
- `TunableCoupler` is a SQUID-based independent subsystem, not a coupling record
- couplings are interaction declarations between named subsystems
- `hamiltonian(sys)` belongs to the Model layer
- model-layer Hamiltonians must be compatible with QuantumToolbox.jl
- simulation code consumes models or Hamiltonians; it does not rebuild model logic
- analysis code consumes typed results; it does not rebuild Hamiltonians

## Domain Model

### Subsystems

Subsystems represent physical degrees of freedom.

Examples:
- `Transmon`
- `TunableTransmon`
- `Resonator`
- `TunableCoupler`

### Couplings

Couplings represent interactions between named subsystems.

Examples:
- `CapacitiveCoupling`

### System Container

`CompositeSystem` stores:
- subsystem declarations
- coupling declarations
- topology metadata

`CompositeSystem` does not store:
- solver state
- cached eigensolutions by default
- separate fixed/tunable variants

## Unit Convention

Phase 1 uses a shared frequency-unit convention with `ħ = 1`.

- `EJ`, `EC`, resonator `ω`, and coupling `g` must all use the same unit
- Hamiltonian eigenvalues and derived transition frequencies use that same unit
- the example notebooks use `GHz`
- `ncut` and `dim` are dimensionless truncation sizes
- tunable `flux` is a dimensionless reduced flux interpreted as `Φ/Φ0`
- `asymmetry` is a dimensionless SQUID asymmetry parameter
- although the resonator parameter is written as `ω`, the current implementation treats it numerically in the same unit convention as the other Hamiltonian parameters

## Layer Ownership

### 1. Architecture Layer

Owns:
- subsystem parameter types
- coupling types
- `CompositeSystem`
- constructor-time validation
- topology and naming rules

Public examples:
- `Transmon`
- `TunableTransmon`
- `Resonator`
- `TunableCoupler`
- `CapacitiveCoupling`
- `CompositeSystem`

Must not own:
- Hamiltonian assembly
- eigensolver calls
- time evolution
- physics-result analysis

### 2. Model Layer

Owns:
- truncation logic
- local operator construction
- interaction-term construction
- Hamiltonian assembly
- normalized internal model representation
- collapse-operator construction when supported

Public examples:
- `build_model(sys)`
- `hamiltonian(sys)`
- `hamiltonian(model)`

Contract:
- `hamiltonian(...)` returns a QuantumToolbox-compatible Hamiltonian object
- the returned object is the canonical handoff from Model to Simulation
- the default effective transmon model uses a Duffing approximation derived from `EJ` and `EC`
- circuit-mode transmon-like models use an explicit charge-basis Hamiltonian with a separate charge cutoff
- circuit-mode capacitive couplings are intentionally unsupported until a dedicated circuit-coupling model exists

Must not own:
- eigensolver execution
- parameter sweep orchestration
- high-level analysis helpers

### 3. Simulation Layer

Owns:
- eigenspectrum execution
- closed-system time evolution
- Lindblad evolution
- steady-state execution
- parameter sweep execution
- typed result objects

Public examples:
- `spectrum(sys; levels=6)`
- `evolve(sys, ψ0, tlist; kwargs...)`
- `steady_state(sys; kwargs...)`
- `simulate_sweep(system, spec)`

Rules:
- simulation may accept a `CompositeSystem` or model object
- if a `CompositeSystem` is passed, simulation must delegate Hamiltonian construction to the Model layer
- do not duplicate Hamiltonian formulas in simulation code

### 4. Analysis Layer

Owns:
- transition-frequency extraction
- anharmonicity computation
- dispersive-shift extraction
- sweep summarization
- table-friendly output adapters

Public examples:
- `transition_frequencies(result)`
- `anharmonicity(result)`
- `dispersive_shift(result)`

Must not own:
- Hamiltonian rebuilding
- solver orchestration

## Data Flow

Expected flow:

1. User creates subsystems and couplings.
2. User assembles them into one `CompositeSystem`.
3. Model layer builds a normalized model and QuantumToolbox-compatible Hamiltonian.
4. Simulation layer runs spectrum, sweep, or dynamics workflows.
5. Analysis layer derives physics quantities from typed results.

Reference flow:

```text
Subsystems + Couplings
    -> CompositeSystem
    -> Model / Hamiltonian
    -> SimulationResult
    -> AnalysisResult
```

## CompositeSystem Contract

Requirements:
- accepts fixed and tunable subsystems together
- accepts couplings that refer to subsystem names
- remains phase-stable as new subsystem families are added
- does not encode solver/backend details

Example:

```julia
sys = CompositeSystem(
    TunableTransmon(:q1, EJmax=20.0, EC=0.25, flux=0.0, asymmetry=0.0, ncut=15),
    Resonator(:r1, ω=6.8, dim=10),
    TunableCoupler(:c1, EJmax=15.0, EC=0.30, flux=0.1, asymmetry=0.0, ncut=11),
    CapacitiveCoupling(:q1, :c1, g=0.08),
    CapacitiveCoupling(:c1, :r1, g=0.05),
)
```

## Suggested Module Layout

```text
src/
  QuantumCircuit.jl
  Architecture/
    Architecture.jl
    Subsystems.jl
    Couplings.jl
    CompositeSystem.jl
    Validation.jl
  Model/
    Model.jl
    Truncation.jl
    Operators.jl
    Hamiltonians.jl
    CollapseOperators.jl
  Simulation/
    Simulation.jl
    Spectrum.jl
    Dynamics.jl
    SteadyState.jl
    Sweeps.jl
    Results.jl
  Analysis/
    Analysis.jl
    Transitions.jl
    Anharmonicity.jl
    Dispersive.jl
    Summaries.jl
```

The exact file split can change. The layer boundaries above should not.

## Type Strategy

Recommended public type categories:

- subsystem types: `Transmon`, `TunableTransmon`, `Resonator`, `TunableCoupler`
- coupling types: `CapacitiveCoupling`
- composition type: `CompositeSystem`
- sweep types: `SweepSpec`, `SubsystemSweepTarget`, `CouplingSweepTarget`
- result types: `SpectrumResult`, `SweepResult`, `DynamicsResult`

Guidelines:
- use immutable structs unless mutation is clearly justified
- use symbolic names like `:q1`, `:r1`, `:c1` for topology references
- keep parameter names explicit and physics-facing
- do not store derived operators inside architecture-layer objects

## Extension Rules

To add a new subsystem family:

1. add an Architecture-layer type
2. add Model-layer operators and Hamiltonian terms
3. reuse Simulation-layer workflows unless a new solver interface is required
4. add Analysis-layer helpers only if new observables are needed

This means future subsystems such as `Fluxonium` should extend the existing
composition and model flow, not introduce a parallel architecture.

## Phase Mapping

Phase 1:
- `Transmon`
- `Resonator`
- `CapacitiveCoupling`
- generic `CompositeSystem`
- static model construction
- `hamiltonian(sys)` and spectrum-oriented results

Phase 2:
- `TunableTransmon`
- `TunableCoupler`
- flux-dependent model construction
- mixed fixed/tunable systems
- parameter sweeps for tunable systems

Phase 2.5:
- stable spectrum and sweep analysis helpers
- short docs and notebook walkthroughs
- explicit documentation of deferred model features

Phase 3A:
- named product-state preparation
- named observables for time-domain workflows
- closed-system dynamics
- narrow subsystem-local drives

Phase 3B:
- collapse operators
- Lindblad dynamics
- steady-state and open-system results

## Agent Rules

- do not introduce a second system container for tunable devices
- do not treat `TunableCoupler` as a coupling declaration
- do not place Hamiltonian-building code outside the Model layer
- do not let analysis functions call solvers directly
- add tests at the layer where behavior is introduced
