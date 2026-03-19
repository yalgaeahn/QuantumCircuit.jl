module Analysis

using ..Simulation: DynamicsResult, ObservableTrace, SpectrumResult, SweepResult
using ..Architecture: CompositeSystem
using ..Model: AbstractHamiltonianSpec, EffectiveHamiltonianSpec, StaticSystemModel

export SweepSeries,
    SweepSummary,
    SpectrumComparisonResult,
    SubspaceSpec,
    UnitaryTraceResult,
    anharmonicity,
    anharmonicity_curve,
    best_cz,
    best_move,
    compare_minimum_gap,
    compare_model_spectra,
    compare_spectra,
    conditional_phase,
    final_state,
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
include("Summaries.jl")
include("Dynamics.jl")
include("Unitaries.jl")
include("Renger2026.jl")

end
