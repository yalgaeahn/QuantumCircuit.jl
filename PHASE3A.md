# Phase 3A Design

## Purpose

This document defines the first dynamics milestone for QuantumCircuit.jl:
closed-system time evolution on top of the completed Phase-1/2/2.5 workflow.

Phase 3A should be implementation-ready. It is not a research wishlist.

## Goal

Deliver a stable public API for:
- named product-state preparation
- named observable specification
- closed-system time evolution
- narrow driven evolution for subsystem-local drives
- supported uncoupled circuit-mode closed-system workflows
- time-dependent flux control for supported circuit and single-device effective workflows
- typed time-domain results that fit the existing `Model -> Simulation -> Analysis` flow

## Non-Goals

Not part of Phase 3A:
- Lindblad dynamics
- collapse operators
- steady-state solvers
- stochastic trajectories
- circuit-mode capacitive couplings
- coupled circuit-mode dynamics examples
- nonadiabatic effective flux control for multi-subsystem or coupled systems
- pulse scheduling
- calibration loops
- branch-tracking analysis

These belong to Phase 3B or later.

## Confirmed Solver Constraints

The current design should wrap the existing `QuantumToolbox.jl` solver surface
that is already available in the project environment.

Confirmed local facts:
- `sesolve(H, ψ0, tlist; alg, e_ops, params, progress_bar, inplace, kwargs...)` exists
- `sesolve` accepts a static operator Hamiltonian or a `Tuple` of operator-function pairs
- `sesolve` returns a `TimeEvolutionSol` with `times`, `times_states`, `states`, and `expect`
- expectation values are produced by passing `e_ops`

Design implication:
- QuantumCircuit.jl should normalize `sesolve` inputs and outputs
- it should not expose raw solver-specific layout as the primary public result API

## Layer Ownership

### Architecture

Owns:
- no new system container
- no drive logic
- no solver logic

### Model

Owns:
- named subsystem-local operator construction
- named product-state construction from subsystem ordering and dimensions
- conversion of narrow drive and flux-control specifications into a `QuantumToolbox`-compatible time-dependent Hamiltonian representation

### Simulation

Owns:
- `evolve(...)`
- solver selection and `sesolve` delegation
- `DynamicsResult`
- mapping observable specifications to solver `e_ops`

### Analysis

Owns:
- trace extraction
- final-state accessors
- population summaries derived from `DynamicsResult`

## Public API Proposal

### Initial States

Provide a named product-state helper that respects subsystem ordering:

```julia
basis_state(system::CompositeSystem; assignments...)
basis_state(model::StaticSystemModel; assignments...)
```

Usage:

```julia
ψ0 = basis_state(sys; q1 = 1, r1 = 0, c1 = 0)
```

Rules:
- unspecified subsystems default to level `0`
- each requested level must be within the local truncation range
- names must match `subsystem_names(system)`

Rationale:
- this keeps users away from manual tensor-product bookkeeping
- it matches the named-topology design used everywhere else in the package

### Observables

Use a small named observable specification rather than exposing raw embedded
operators as the only high-level interface.

```julia
struct ObservableSpec
    label::Symbol
    target::Symbol
    operator::Symbol
end
```

Supported operator symbols for Phase 3A:
- effective-mode oscillator operators: `:a`, `:adag`, `:n`, `:x`, `:y`
- circuit-mode transmon-like operators: `:charge`, `:cosphi`, `:sinphi`

Operator validity depends on the active Hamiltonian specification. Circuit-mode
examples should use the circuit-native symbols where appropriate.

Examples:

```julia
observables = [
    ObservableSpec(:nq, :q1, :n),
    ObservableSpec(:phi_q, :q1, :cosphi),
]
```

Low-level escape hatch:
- `evolve(...; e_ops = raw_quantumtoolbox_ops)` may still be supported internally
- but the documented API should prefer `ObservableSpec`

### Drives

Keep Phase-3A drives narrow and subsystem-local.

```julia
struct SubsystemDrive{F}
    label::Symbol
    target::Symbol
    operator::Symbol
    coefficient::F
end
```

Rules:
- `target` names a subsystem
- `operator` uses the same symbol set as `ObservableSpec`
- `coefficient` should accept `(params, t)` and return a scalar amplitude
- drives are additive terms on top of the static Hamiltonian

Example:

```julia
drive = SubsystemDrive(
    :q1_x_drive,
    :q1,
    :x,
    (p, t) -> p.Ω * cos(p.ωd * t),
)
```

Non-goals for Phase 3A:
- coupling drives
- pulse objects
- waveform compilation

### Flux Controls

Keep time-dependent flux control separate from operator drives.

```julia
FluxControl(label::Symbol, target::Symbol, modulation; derivative)
```

Rules:
- `target` must name a `TunableTransmon` or `TunableCoupler`
- `modulation(params, t)` returns additive reduced flux `ΔΦ/Φ0`
- `derivative(params, t)` is required for the nonadiabatic effective single-device workflow
- circuit-mode flux-control examples must remain uncoupled
- effective nonadiabatic flux-control examples must remain single-subsystem and uncoupled

Example:

```julia
flux_control = FluxControl(
    :q1_flux_pulse,
    :q1,
    (p, t) -> p.δ * sin(p.ωf * t);
    derivative = (p, t) -> p.δ * p.ωf * cos(p.ωf * t),
)
```

### Evolution API

Provide one main closed-system evolution entry point:

```julia
evolve(system::CompositeSystem, ψ0, tlist; observables = nothing, drives = nothing, flux_controls = nothing, params = NamedTuple(), alg = nothing, progress_bar = Val(false), inplace = Val(true), kwargs...)
evolve(model::StaticSystemModel, ψ0, tlist; kwargs...)
```

Behavior:
- `system` input builds the model first
- `model` input reuses the existing static model directly
- Hamiltonian assembly for `drives` and `flux_controls` stays in the Model layer
- simulation delegates to `QuantumToolbox.sesolve`

Important wrapper rule:
- if states are meant to be returned, the wrapper should force `saveat = tlist`
- this avoids solver defaults that may only save the final state when expectation operators are present

### Result Types

Normalize raw solver output into typed result objects.

```julia
struct ObservableTrace{T}
    label::Symbol
    times::Vector{Float64}
    values::Vector{T}
end

struct DynamicsResult{M,State,Solver}
    model::M
    times::Vector{Float64}
    states::Vector{State}
    observables::Vector{ObservableTrace}
    solver_result::Solver
end
```

Notes:
- `solver_result` preserves access to raw `TimeEvolutionSol` when needed
- `observables` should hide the raw storage layout returned by `sesolve`
- `states` should be present in Phase 3A by default

## Internal Model Utilities

Phase 3A needs model-layer helpers that are currently implicit.

Recommended additions:
- `annihilation_operator(model, name::Symbol)`
- `creation_operator(model, name::Symbol)`
- `number_operator(model, name::Symbol)`
- `_quadrature_operator(model, name::Symbol, axis::Symbol)` for internal drive resolution
- `basis_state(model; assignments...)`
- `_time_dependent_hamiltonian(model, drives, flux_controls)` returning either a static operator or a `Tuple` compatible with `sesolve`

Design rule:
- operator embedding stays a Model concern
- Simulation should never rebuild local operators itself

## Supported Circuit-Mode Workflows

Phase 3A supports these circuit-mode workflows:
- uncoupled `CircuitHamiltonianSpec(charge_cutoff = ...)`
- spectrum and flux sweeps on transmon-like subsystems
- circuit-native operators and observables through `charge`, `cosphi`, and `sinphi`
- basis-state preparation and closed-system dynamics on supported uncoupled systems
- tunable-device flux control in the supported uncoupled setting

Deferred circuit-mode features remain out of scope:
- capacitive couplings in circuit mode
- coupled circuit-mode dynamics examples
- circuit-mode reductions or Schrieffer-Wolff-style derived models

## Expected File Layout

```text
src/
  Model/
    Model.jl
    Operators.jl
    States.jl
    Drives.jl
  Simulation/
    Simulation.jl
    Dynamics.jl
    Results.jl
  Analysis/
    Analysis.jl
    Dynamics.jl
```

Exact filenames can vary, but the ownership split should stay the same.

## Analysis Surface for Phase 3A

Keep the first dynamics analysis layer minimal.

Recommended helpers:

```julia
observable_trace(result::DynamicsResult, label::Symbol)
final_state(result::DynamicsResult)
```

Optional if states are available:

```julia
population_trace(result::DynamicsResult, subsystem::Symbol, level::Int)
```

Do not add in Phase 3A:
- FFT or frequency-domain toolkits
- calibration heuristics
- open-system summaries

## Example User Flow

```julia
q1 = TunableTransmon(:q1; EJmax = 20.0, EC = 0.25, flux = 0.10, asymmetry = 0.20, ng = 0.05, ncut = 6)
sys = CompositeSystem(q1)
circuit_spec = CircuitHamiltonianSpec(charge_cutoff = 3)

ψ0 = basis_state(sys; hamiltonian_spec = circuit_spec, q1 = 0)
tlist = collect(range(0.0, 20.0; length = 201))

flux_control = FluxControl(
    :q1_flux_pulse,
    :q1,
    (p, t) -> p.δ * sin(p.ωf * t);
    derivative = (p, t) -> p.δ * p.ωf * cos(p.ωf * t),
)

result = evolve(
    sys,
    ψ0,
    tlist;
    hamiltonian_spec = circuit_spec,
    flux_controls = [flux_control],
    observables = [ObservableSpec(:charge_q, :q1, :charge)],
    params = (; δ = 0.06, ωf = 5.5),
)

trace = observable_trace(result, :charge_q)
```

## Test Plan

Mandatory regression coverage:

1. Static free evolution
- evolve a stationary state under a diagonal static Hamiltonian
- verify the state norm is preserved
- verify a conserved number observable stays constant

2. Named basis-state preparation
- verify subsystem names map to the correct tensor-product ordering
- verify invalid subsystem names and out-of-range levels throw helpful errors

3. Observable mapping
- verify `ObservableSpec(:nq, :q1, :n)` matches manual expectation values
- verify labels remain stable and lookup works

4. Driven closed-system example
- run a small driven two-level-style example
- verify observable traces have the same length as `tlist`
- verify the wrapper path with `drives` reaches `sesolve` successfully

5. Circuit-mode closed-system example
- run an uncoupled circuit-mode example with `CircuitHamiltonianSpec(charge_cutoff = ...)`
- verify circuit-native observables such as `:charge` or `:cosphi` can be traced
- verify the documented notebook flow stays within supported uncoupled circuit-mode limits

6. Flux-control example
- run an uncoupled tunable-device example with `FluxControl`
- verify the circuit-mode flux-control path reaches `sesolve`
- verify the single-device nonadiabatic effective flux path is documented with its explicit limits

7. Result shape
- verify `DynamicsResult.times`, `states`, and `observables` stay aligned
- verify `solver_result` is preserved for debugging

## Documentation Deliverables

Phase 3A should ship with:
- one README section or short note pointing to dynamics support
- one general notebook walkthrough for closed-system dynamics
- one dedicated notebook walkthrough for supported circuit-Hamiltonian dynamics
- one regression test mirroring the documented example

## Risks

Primary risks:
- overdesigning drive abstractions before a real use case exists
- leaking `QuantumToolbox` result layout into the public API
- putting operator construction into Simulation instead of Model
- letting open-system concerns distort the first closed-system milestone

Mitigation:
- keep drives subsystem-local in Phase 3A
- always normalize solver output through `DynamicsResult`
- build model utilities first, then wire `evolve`

## Recommended Implementation Order

1. Add named basis-state helpers.
2. Add named observable/operator helpers in the Model layer.
3. Implement static `evolve(...; observables=...)` without drives.
4. Define `DynamicsResult` and minimal analysis helpers.
5. Add narrow `SubsystemDrive` support for time-dependent closed-system evolution.
6. Add notebook and regression tests.
