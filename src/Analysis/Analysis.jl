module Analysis

using QuantumToolbox: Ket, QuantumObject, dag, eigenstates, kron
using ..Simulation:
    DynamicsResult,
    EigensystemResult,
    EigensystemSweepResult,
    ObservableTrace,
    SpectrumResult,
    SweepResult,
    eigensystem
using ..Architecture: CompositeSystem, name, subsystems
using ..Model: AbstractHamiltonianSpec, EffectiveHamiltonianSpec, StaticSystemModel, build_model

export SweepSeries,
    BareProductBasis,
    BareSubsystemBasis,
    LabelMap,
    LabeledEigensystemSweepResult,
    SweepSummary,
    SpectrumComparisonResult,
    SubspaceSpec,
    UnitaryTraceResult,
    anharmonicity,
    anharmonicity_curve,
    bare_label,
    bare_product_basis,
    best_cz,
    best_move,
    compare_minimum_gap,
    compare_model_spectra,
    compare_spectra,
    conditional_phase,
    dressed_index,
    dressed_state_components,
    energy_by_bare_label,
    energy_curve,
    final_state,
    label_sweep,
    load_renger2026_snapshot,
    minimum_gap,
    observable_trace,
    population_trace,
    projected_unitary,
    renger2026_model_pair,
    renger2026_reduced_system,
    renger2026_stage1_qcr_system,
    renger2026_stage1_qr_system,
    strip_local_z_phases,
    subspace_spec,
    sweep_summary,
    transition_curve,
    transition_frequencies

include("Spectrum.jl")
include("Sweeps.jl")
include("Labeling.jl")
include("Summaries.jl")
include("Dynamics.jl")
include("Unitaries.jl")
include("Renger2026.jl")

end
