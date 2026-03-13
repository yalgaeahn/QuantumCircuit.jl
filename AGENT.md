# QuantumCircuit.jl Agent Guide

## Read Order

When working in this repository, read documents in this order:

1. `AGENT.md`
2. `PLAN.md`
3. `ARCHITECTURE.md`
4. `ROADMAP.md`

## Source of Truth

Use each document for its intended purpose:

- `PLAN.md`: product scope, phase boundaries, milestone definitions
- `ARCHITECTURE.md`: domain model, layer ownership, API boundaries
- `ROADMAP.md`: implementation order and task decomposition

If documents appear to conflict:
- `ARCHITECTURE.md` wins for code structure and API ownership
- `ROADMAP.md` wins for sequencing
- `PLAN.md` wins for scope and phase intent
- update the docs to remove the conflict rather than silently choosing one

## Project Invariants

- `CompositeSystem` is the single system-composition API
- fixed and tunable subsystems must coexist in one `CompositeSystem`
- `TunableCoupler` is a SQUID-based independent subsystem
- `CapacitiveCoupling` is an interaction declaration, not a subsystem
- `hamiltonian(sys)` belongs to the Model layer
- `hamiltonian(sys)` must return a QuantumToolbox-compatible object
- the initial Phase-1 transmon path uses a Duffing approximation from `EJ` and `EC`
- Simulation-layer code may consume systems or models, but must delegate Hamiltonian construction to the Model layer
- Analysis-layer code must consume typed results and must not call solvers directly

## Working Rules

- implement lower layers before higher layers
- do not introduce a tunable-only system container
- do not duplicate Hamiltonian formulas across layers
- prefer explicit typed structs over ad hoc dictionaries and tuples
- add tests in the same change where behavior is introduced
- if a design decision changes the domain model or layer boundary, update the docs in the same change

## Phase Discipline

- Phase 1: static core only
- Phase 2: tunable devices on top of the existing core
- Phase 3: dynamics and open-system workflows

Do not pull Phase-2 or Phase-3 abstractions into Phase-1 code unless they are
required to keep a stable generic interface.

## Preferred Workflow

1. Identify the target phase in `PLAN.md`.
2. Confirm layer ownership in `ARCHITECTURE.md`.
3. Choose the next dependency-safe task from `ROADMAP.md`.
4. Implement the smallest coherent slice.
5. Add or update tests.
6. If the implementation changes design intent, update the relevant docs.
