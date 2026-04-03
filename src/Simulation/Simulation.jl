module Simulation

using QuantumToolbox: Ket, QuantumObject, eigenenergies, eigenstates, sesolve
using ..Architecture: CompositeSystem, with_coupling_parameter, with_subsystem_parameter
using ..Model:
    AbstractHamiltonianSpec,
    EffectiveHamiltonianSpec,
    FluxControl,
    ObservableSpec,
    StaticSystemModel,
    SubsystemDrive,
    _embedded_operator,
    _time_dependent_hamiltonian,
    basis_state,
    build_model,
    hamiltonian

export CouplingSweepTarget,
    DynamicsResult,
    EigensystemResult,
    EigensystemSweepResult,
    ObservableTrace,
    SpectrumResult,
    SubsystemSweepTarget,
    SweepResult,
    SweepSpec,
    eigensystem,
    evolve,
    simulate_eigensystem_sweep,
    simulate_sweep,
    spectrum

struct SpectrumResult{M}
    model::M
    energies::Vector{Float64}
end

struct EigensystemResult{M, S}
    model::M
    energies::Vector{Float64}
    states::Vector{S}
end

abstract type AbstractSweepTarget end

struct SubsystemSweepTarget <: AbstractSweepTarget
    name::Symbol
end

struct CouplingSweepTarget <: AbstractSweepTarget
    source::Symbol
    target::Symbol
end

struct SweepSpec{Target<:AbstractSweepTarget, T}
    target::Target
    parameter::Symbol
    values::Vector{T}
    levels::Int
end

struct SweepResult{Spec, T, Result}
    base_system::CompositeSystem
    spec::Spec
    values::Vector{T}
    spectra::Vector{Result}
end

struct EigensystemSweepResult{Spec, T, Result}
    base_system::CompositeSystem
    spec::Spec
    values::Vector{T}
    spectra::Vector{Result}
end

struct ObservableTrace{T}
    label::Symbol
    times::Vector{Float64}
    values::Vector{T}
end

struct DynamicsResult{M, State, Trace, Solver}
    model::M
    times::Vector{Float64}
    states::Vector{State}
    observables::Vector{Trace}
    solver_result::Solver
end

function SweepSpec(target::SubsystemSweepTarget, parameter::Symbol, values::AbstractVector; levels::Integer = 6)
    isempty(values) && throw(ArgumentError("SweepSpec requires at least one sweep value."))
    levels > 0 || throw(ArgumentError("levels must be a positive integer."))
    collected_values = collect(values)
    return SweepSpec(target, parameter, collected_values, Int(levels))
end

function SweepSpec(target::CouplingSweepTarget, parameter::Symbol, values::AbstractVector; levels::Integer = 6)
    isempty(values) && throw(ArgumentError("SweepSpec requires at least one sweep value."))
    levels > 0 || throw(ArgumentError("levels must be a positive integer."))
    collected_values = collect(values)
    return SweepSpec(target, parameter, collected_values, Int(levels))
end

SweepSpec(name::Symbol, parameter::Symbol, values::AbstractVector; levels::Integer = 6) =
    SweepSpec(SubsystemSweepTarget(name), parameter, values; levels = levels)

function SweepSpec(source::Symbol, target::Symbol, parameter::Symbol, values::AbstractVector; levels::Integer = 6)
    return SweepSpec(CouplingSweepTarget(source, target), parameter, values; levels = levels)
end

function spectrum(model::StaticSystemModel; levels::Integer = 6)
    levels > 0 || throw(ArgumentError("levels must be a positive integer."))
    energies = sort(Float64.(eigenenergies(hamiltonian(model))))
    levels <= length(energies) || throw(ArgumentError("levels exceeds the available Hilbert-space dimension."))
    return SpectrumResult(model, energies[1:levels])
end

function spectrum(
    system::CompositeSystem;
    levels::Integer = 6,
    hamiltonian_spec::AbstractHamiltonianSpec = EffectiveHamiltonianSpec(),
)
    return spectrum(build_model(system; hamiltonian_spec = hamiltonian_spec); levels = levels)
end

function eigensystem(model::StaticSystemModel; levels::Integer = 6)
    levels > 0 || throw(ArgumentError("levels must be a positive integer."))
    eigensolve = eigenstates(hamiltonian(model))
    energies = Float64.(real.(eigensolve.values))
    order = sortperm(energies)
    levels <= length(order) || throw(ArgumentError("levels exceeds the available Hilbert-space dimension."))
    dims = basis_state(model).dims
    ordered_states = [QuantumObject(eigensolve.vectors[:, index]; type = Ket(), dims = dims) for index in order[1:levels]]
    return EigensystemResult(model, energies[order][1:levels], ordered_states)
end

function eigensystem(
    system::CompositeSystem;
    levels::Integer = 6,
    hamiltonian_spec::AbstractHamiltonianSpec = EffectiveHamiltonianSpec(),
)
    return eigensystem(build_model(system; hamiltonian_spec = hamiltonian_spec); levels = levels)
end

function simulate_sweep(
    system::CompositeSystem,
    spec::SweepSpec;
    hamiltonian_spec::AbstractHamiltonianSpec = EffectiveHamiltonianSpec(),
)
    spectra = [
        spectrum(
            _apply_sweep(system, spec.target, spec.parameter, value);
            levels = spec.levels,
            hamiltonian_spec = hamiltonian_spec,
        ) for value in spec.values
    ]
    return SweepResult(system, spec, copy(spec.values), spectra)
end

function simulate_eigensystem_sweep(
    system::CompositeSystem,
    spec::SweepSpec;
    hamiltonian_spec::AbstractHamiltonianSpec = EffectiveHamiltonianSpec(),
)
    spectra = [
        eigensystem(
            _apply_sweep(system, spec.target, spec.parameter, value);
            levels = spec.levels,
            hamiltonian_spec = hamiltonian_spec,
        ) for value in spec.values
    ]
    return EigensystemSweepResult(system, spec, copy(spec.values), spectra)
end

function _apply_sweep(system::CompositeSystem, target::SubsystemSweepTarget, parameter::Symbol, value)
    return with_subsystem_parameter(system, target.name, parameter, value)
end

function _apply_sweep(system::CompositeSystem, target::CouplingSweepTarget, parameter::Symbol, value)
    return with_coupling_parameter(system, target.source, target.target, parameter, value)
end

function evolve(
    system::CompositeSystem,
    ψ0,
    tlist::AbstractVector;
    hamiltonian_spec::AbstractHamiltonianSpec = EffectiveHamiltonianSpec(),
    observables = nothing,
    drives = nothing,
    flux_controls = nothing,
    params = NamedTuple(),
    alg = nothing,
    progress_bar::Union{Val, Bool} = Val(false),
    inplace::Union{Val, Bool} = Val(true),
    kwargs...,
)
    return evolve(
        build_model(system; hamiltonian_spec = hamiltonian_spec),
        ψ0,
        tlist;
        observables = observables,
        drives = drives,
        flux_controls = flux_controls,
        params = params,
        alg = alg,
        progress_bar = progress_bar,
        inplace = inplace,
        kwargs...,
    )
end

function evolve(
    model::StaticSystemModel,
    ψ0,
    tlist::AbstractVector;
    observables = nothing,
    drives = nothing,
    flux_controls = nothing,
    params = NamedTuple(),
    alg = nothing,
    progress_bar::Union{Val, Bool} = Val(false),
    inplace::Union{Val, Bool} = Val(true),
    kwargs...,
)
    isempty(tlist) && throw(ArgumentError("tlist must contain at least one time point."))
    haskey(kwargs, :saveat) &&
        throw(ArgumentError("evolve manages saveat internally to keep states aligned with the requested tlist."))

    observable_specs = _normalize_observables(observables)
    _validate_unique_observable_labels(observable_specs)

    H = _time_dependent_hamiltonian(model, drives, flux_controls)
    e_ops = isempty(observable_specs) ? nothing : [_embedded_operator(model, spec.target, spec.operator) for spec in observable_specs]

    solver_kwargs = (; e_ops, params, progress_bar, inplace, saveat = tlist, kwargs...)
    solver_result =
        isnothing(alg) ? sesolve(H, ψ0, tlist; solver_kwargs...) :
        sesolve(H, ψ0, tlist; alg = alg, solver_kwargs...)

    observable_traces = _build_observable_traces(observable_specs, solver_result)
    states = collect(solver_result.states)
    times = Float64.(solver_result.times_states)

    length(states) == length(times) ||
        throw(ArgumentError("Dynamics result states and saved times must have the same length."))

    return DynamicsResult(model, times, states, observable_traces, solver_result)
end

function _normalize_observables(observables::Nothing)
    return ObservableSpec[]
end

function _normalize_observables(observables::Union{Tuple, AbstractVector})
    normalized = ObservableSpec[]
    for observable in observables
        observable isa ObservableSpec ||
            throw(ArgumentError("observables must contain only ObservableSpec values, got $(typeof(observable))."))
        push!(normalized, observable)
    end
    return normalized
end

function _validate_unique_observable_labels(observables::Vector{ObservableSpec})
    labels = [observable.label for observable in observables]
    length(unique(labels)) == length(labels) ||
        throw(ArgumentError("Observable labels must be unique within one evolve call."))
    return nothing
end

function _build_observable_traces(observables::Vector{ObservableSpec}, solver_result)
    isempty(observables) && return ObservableTrace[]
    solver_result.expect === nothing &&
        throw(ArgumentError("Expected solver expectation values because observables were requested."))

    traces = ObservableTrace[]
    expect_values = solver_result.expect
    times = Float64.(solver_result.times)

    for (index, observable) in enumerate(observables)
        push!(traces, ObservableTrace(observable.label, times, collect(vec(expect_values[index, :]))))
    end

    return traces
end

end
