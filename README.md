# QuantumCircuit.jl

QuantumCircuit.jl is a superconducting quantum-circuit Hamiltonian simulator
built on top of `QuantumToolbox.jl`.

Current scope:
- fixed and flux-tunable subsystem declarations
- composite-system Hamiltonian construction
- explicit effective and circuit Hamiltonian specifications
- exact circuit-mode capacitive coupling through `CircuitCapacitiveCoupling(...; G = ...)` for transmon-like↔transmon-like and resonator↔transmon-like pairs
- low-lying eigenspectra
- parameter sweeps
- result-driven spectrum and sweep analysis
- closed-system time evolution
- subsystem-local driven closed-system evolution
- supported time-dependent flux control for exact circuit-mode workflows and the single-device nonadiabatic effective workflow

Not yet in scope:
- Lindblad or steady-state workflows
- exact circuit-mode resonator↔resonator coupling
- paper-derived effective time-dependent couplings such as `g(t)` and `ḡ(t)`
- automatic circuit quantization

## Project Status

Status as of 2026-03-19:
- Phase 1: complete
- Phase 2: complete
- Phase 2.5: complete
- Phase 3A: complete

## Quickstart

The examples below use a shared frequency unit such as `GHz`.

## Notebook Walkthroughs

- `output/jupyter-notebook/phase1-static-core-walkthrough.ipynb`: Phase 1 fixed-system walkthrough.
- `output/jupyter-notebook/phase2-tunable-sweep-walkthrough.ipynb`: Phase 2 tunable sweep walkthrough with the current Analysis helpers.
- `output/jupyter-notebook/phase2-analysis-recipes.ipynb`: focused recipes for reading `SweepResult`, `SweepSeries`, and `SweepSummary`.
- `output/jupyter-notebook/phase3a-closed-system-dynamics-walkthrough.ipynb`: Phase 3A closed-system dynamics walkthrough with `basis_state`, `ObservableSpec`, `SubsystemDrive`, and `evolve`.
- `output/jupyter-notebook/phase3a-circuit-hamiltonian-dynamics.ipynb`: supported uncoupled circuit-Hamiltonian walkthrough with `CircuitHamiltonianSpec`, circuit-native operators, basis states, dynamics, and `FluxControl`.
- `output/jupyter-notebook/phase5-circuit-capacitive-coupling-walkthrough.ipynb`: coupled exact circuit-Hamiltonian walkthrough with `CircuitCapacitiveCoupling`, `:G` sweeps, local-drive transfer, and coupled `FluxControl`.
- `output/jupyter-notebook/phase3a-drive-dynamics-recipes.ipynb`: focused drive-dynamics recipes for comparing undriven and driven closed-system evolution.
- `output/jupyter-notebook/phase3a-detuning-and-coupled-response.ipynb`: single-mode drive detuning notebook with a short coupled-response extension.
- `output/jupyter-notebook/phase3a-fixed-qubit-coupler-resonator-dynamics.ipynb`: mixed `Transmon - TunableCoupler - Resonator` dynamics walkthrough with static flux selection and transferred resonator response.
- `output/jupyter-notebook/phase3a-tunable-coupler-two-qubit-dynamics.ipynb`: mixed `TunableTransmon - TunableCoupler - Transmon` dynamics and population-transfer walkthrough.

The Phase 3A dynamics notebooks use `UnicodePlots.jl` for lightweight inline plots inside the notebooks.
Use `phase3a-closed-system-dynamics-walkthrough.ipynb` for the baseline dynamics API, `phase3a-circuit-hamiltonian-dynamics.ipynb` for uncoupled charge-basis workflows and circuit-native operators, `phase5-circuit-capacitive-coupling-walkthrough.ipynb` for exact coupled circuit-mode capacitive workflows, `phase3a-detuning-and-coupled-response.ipynb` for drive-detuning intuition, `phase3a-fixed-qubit-coupler-resonator-dynamics.ipynb` for fixed-qubit transfer into a resonator, and `phase3a-tunable-coupler-two-qubit-dynamics.ipynb` for mixed two-qubit transfer.
Smoke-check the two canonical Phase 3A notebooks with `python3 scripts/smoke_phase3a_notebooks.py`.

### 1. Static transmon spectrum

```julia
using QuantumCircuit

q1 = Transmon(:q1; EJ = 20.0, EC = 0.25, ncut = 6)
sys = CompositeSystem(q1)

spec = spectrum(sys; levels = 4)
energies = spec.energies
freqs = transition_frequencies(spec)
alpha = anharmonicity(spec)

(; energies, ω01 = freqs[1], α = alpha)
```

### 2. Tunable transmon flux sweep

```julia
using QuantumCircuit

tq = TunableTransmon(:tq; EJmax = 20.0, EC = 0.25, flux = 0.0, asymmetry = 0.0, ncut = 6)
sys = CompositeSystem(tq)

sweep = SweepSpec(:tq, :flux, [0.0, 0.15, 0.30]; levels = 4)
result = simulate_sweep(sys, sweep)

ω01_curve = transition_curve(result)
alpha_curve = anharmonicity_curve(result)
gap = minimum_gap(result)
summary = sweep_summary(result)

(;
    flux_values = summary.values,
    ω01 = ω01_curve.data,
    α = alpha_curve.data,
    minimum_gap = gap,
)
```

### 3. Mixed fixed/tunable system

```julia
using QuantumCircuit

q_fixed = Transmon(:qf; EJ = 19.0, EC = 0.24, ncut = 5)
coupler = TunableCoupler(:c1; EJmax = 15.0, EC = 0.30, flux = 0.0, asymmetry = 0.1, ncut = 4)
resonator = Resonator(:r1; ω = 6.7, dim = 3)

sys = CompositeSystem(
    q_fixed,
    coupler,
    resonator,
    CapacitiveCoupling(:qf, :c1; g = 0.07),
    CapacitiveCoupling(:c1, :r1; g = 0.05),
)

spec = spectrum(sys; levels = 5)
flux_sweep = SweepSpec(:c1, :flux, [0.0, 0.10, 0.20]; levels = 5)
result = simulate_sweep(sys, flux_sweep)
summary = sweep_summary(result)

(;
    low_lying_energies = spec.energies,
    coupler_flux_values = summary.values,
    transition_01 = summary.transition_01,
)
```

### 4. Coupling sweep in a static system

```julia
using QuantumCircuit

q1 = Transmon(:q1; EJ = 20.0, EC = 0.25, ncut = 5)
r1 = Resonator(:r1; ω = 6.8, dim = 4)

sys = CompositeSystem(
    q1,
    r1,
    CapacitiveCoupling(:q1, :r1; g = 0.02),
)

g_sweep = SweepSpec(:q1, :r1, :g, [0.02, 0.08, 0.14]; levels = 5)
result = simulate_sweep(sys, g_sweep)
ω01_curve = transition_curve(result)
summary = sweep_summary(result)

(;
    coupling_values = summary.values,
    ω01 = ω01_curve.data,
    weakest_coupling_spectrum = result.spectra[1].energies,
)
```

### 5. Driven closed-system dynamics

```julia
using QuantumCircuit

resonator = Resonator(:r1; ω = 1.0, dim = 2)
sys = CompositeSystem(resonator)
ψ0 = basis_state(sys; r1 = 0)
tlist = collect(range(0.0, 6.0; length = 61))

drive = SubsystemDrive(
    :r1_x_drive,
    :r1,
    :x,
    (p, t) -> p.Ω * cos(p.ωd * t),
)

result = evolve(
    sys,
    ψ0,
    tlist;
    drives = [drive],
    observables = [ObservableSpec(:nr, :r1, :n)],
    params = (; Ω = 0.35, ωd = 1.0),
)

trace = observable_trace(result, :nr)
p0 = population_trace(result, :r1, 0)
p1 = population_trace(result, :r1, 1)

(;
    times = trace.times,
    driven_population = real.(trace.values),
    ground_population = p0.values,
    excited_population = p1.values,
    final_state = final_state(result),
)
```

## Reading Analysis Results

- `transition_curve(result)` returns a `SweepSeries` with `values`, `data`, and a descriptive `label`.
- `anharmonicity_curve(result)` returns the same shape, specialized for anharmonicity.
- `minimum_gap(result)` returns the smallest transition gap and the sweep value where it occurs.
- `sweep_summary(result)` returns a compact `SweepSummary` with the default `0 -> 1` transition and anharmonicity vectors.
- `observable_trace(result, label)` returns a named time trace from `DynamicsResult`.
- `population_trace(result, subsystem, level)` returns a level population trace derived from the saved dynamics states.
- `final_state(result)` returns the last saved quantum state from a dynamics run.

## Modeling Notes

- The default transmon and tunable-device path uses a Duffing-style effective Hamiltonian through `EffectiveHamiltonianSpec()`.
- `CircuitHamiltonianSpec(charge_cutoff = N)` builds a charge-basis Hamiltonian with Hilbert dimension `2N + 1` for each transmon-like subsystem and supports exact circuit-mode couplings through `CircuitCapacitiveCoupling(...; G = ...)`.
- `EJ`, `EC`, resonator `ω`, effective coupling `g`, and exact circuit coupling `G` should share the same frequency unit.
- `flux` is interpreted as reduced flux `Phi/Phi0`.
- subsystem `ncut` remains the effective-model Hilbert dimension; circuit-mode charge cutoff is configured separately on `CircuitHamiltonianSpec`.
- `ng` is used in circuit mode for `Transmon`, `TunableTransmon`, and `TunableCoupler`.
- `FluxControl(...)` is available for supported exact circuit-mode workflows and for the single-device nonadiabatic effective workflow.
- `CapacitiveCoupling(...; g = ...)` remains the effective-mode exchange coupling; use `CircuitCapacitiveCoupling(...; G = ...)` for exact circuit mode.
- supported exact circuit-mode pair families are transmon-like↔transmon-like and resonator↔transmon-like; resonator↔resonator remains out of scope.
- `output/jupyter-notebook/phase3a-circuit-hamiltonian-dynamics.ipynb` is the primary uncoupled circuit-mode walkthrough, and `output/jupyter-notebook/phase5-circuit-capacitive-coupling-walkthrough.ipynb` is the primary coupled exact-circuit walkthrough.
