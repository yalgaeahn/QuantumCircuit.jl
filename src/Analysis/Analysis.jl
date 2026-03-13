module Analysis

using ..Simulation: DynamicsResult, ObservableTrace, SpectrumResult, SweepResult

export SweepSeries,
    SweepSummary,
    anharmonicity,
    anharmonicity_curve,
    final_state,
    minimum_gap,
    observable_trace,
    population_trace,
    sweep_summary,
    transition_curve,
    transition_frequencies

include("Spectrum.jl")
include("Sweeps.jl")
include("Summaries.jl")
include("Dynamics.jl")

end
