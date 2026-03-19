using LinearAlgebra: Diagonal
using QuantumToolbox: Ket, QuantumObject, dag, eigenstates
using ..Architecture: CompositeSystem
using ..Model:
    AbstractHamiltonianSpec,
    EffectiveHamiltonianSpec,
    StaticSystemModel,
    build_model,
    hamiltonian,
    basis_state
using ..Simulation: spectrum
using ..Simulation: evolve

struct SubspaceSpec{State}
    basis::Symbol
    targets::Vector{Symbol}
    level_sets::Vector{Vector{Int}}
    model_subsystem_order::Vector{Symbol}
    model_dimensions::Dict{Symbol, Int}
    labels::Vector{String}
    assignments::Vector{Dict{Symbol, Int}}
    reference_states::Vector{State}
    reference_energies::Vector{Float64}
    reference_overlaps::Vector{Float64}
end

struct UnitaryTraceResult{M, S}
    model::M
    spec::S
    times::Vector{Float64}
    unitaries::Vector{Matrix{ComplexF64}}
    leakages::Matrix{Float64}
    phase_convention::Symbol
end

struct SpectrumComparisonResult
    reference::SpectrumResult
    candidate::SpectrumResult
    level_deltas::Vector{Float64}
    transition_deltas::Vector{Float64}
end

function subspace_spec(
    system::CompositeSystem;
    subsystem_levels,
    hamiltonian_spec::AbstractHamiltonianSpec = EffectiveHamiltonianSpec(),
    basis::Symbol = :bare,
)
    return subspace_spec(build_model(system; hamiltonian_spec = hamiltonian_spec); subsystem_levels = subsystem_levels, basis = basis)
end

function subspace_spec(model::StaticSystemModel; subsystem_levels, basis::Symbol = :bare)
    basis == :bare || basis == :dressed_static ||
        throw(ArgumentError("basis must be :bare or :dressed_static, got $basis."))

    targets, level_sets = _normalize_subsystem_levels(model, subsystem_levels)
    assignments = _enumerate_assignments(targets, level_sets)
    labels = [_assignment_label(targets, assignment) for assignment in assignments]
    bare_states = [_assignment_state(model, assignment) for assignment in assignments]

    reference_states, reference_energies, reference_overlaps =
        basis == :bare ? _bare_reference_bundle(bare_states) : _dressed_reference_bundle(model, bare_states)

    return SubspaceSpec(
        basis,
        targets,
        level_sets,
        copy(model.subsystem_order),
        copy(model.dimensions),
        labels,
        assignments,
        reference_states,
        reference_energies,
        reference_overlaps,
    )
end

function projected_unitary(
    system::CompositeSystem,
    spec::SubspaceSpec,
    tlist::AbstractVector;
    hamiltonian_spec::AbstractHamiltonianSpec = EffectiveHamiltonianSpec(),
    kwargs...,
)
    return projected_unitary(build_model(system; hamiltonian_spec = hamiltonian_spec), spec, tlist; kwargs...)
end

function projected_unitary(
    model::StaticSystemModel,
    spec::SubspaceSpec,
    tlist::AbstractVector;
    drives = nothing,
    flux_controls = nothing,
    params = NamedTuple(),
    alg = nothing,
    progress_bar::Union{Val, Bool} = Val(false),
    inplace::Union{Val, Bool} = Val(true),
    phase_convention::Symbol = :lab,
    kwargs...,
)
    _validate_subspace_spec(model, spec)
    phase_convention == :lab || phase_convention == :strip_local_z ||
        throw(ArgumentError("phase_convention must be :lab or :strip_local_z, got $phase_convention."))
    phase_convention == :lab || length(spec.reference_states) == 4 ||
        throw(ArgumentError(":strip_local_z phase convention requires a four-state subspace."))

    basis_results = [
        evolve(
            model,
            ψ0,
            tlist;
            drives = drives,
            flux_controls = flux_controls,
            params = params,
            alg = alg,
            progress_bar = progress_bar,
            inplace = inplace,
            kwargs...,
        ) for ψ0 in spec.reference_states
    ]

    times = copy(first(basis_results).times)
    subspace_dimension = length(spec.reference_states)
    unitaries = Vector{Matrix{ComplexF64}}(undef, length(times))
    leakages = Matrix{Float64}(undef, subspace_dimension, length(times))

    for time_index in eachindex(times)
        projected = Matrix{ComplexF64}(undef, subspace_dimension, subspace_dimension)

        for column in 1:subspace_dimension
            state = basis_results[column].states[time_index]
            projected_weight = 0.0

            for row in 1:subspace_dimension
                amplitude = dag(spec.reference_states[row]) * state
                projected[row, column] = amplitude
                projected_weight += abs2(amplitude)
            end

            leakages[column, time_index] = max(0.0, 1.0 - projected_weight)
        end

        unitaries[time_index] = phase_convention == :lab ? projected : strip_local_z_phases(projected)
    end

    return UnitaryTraceResult(model, spec, times, unitaries, leakages, phase_convention)
end

function strip_local_z_phases(U::AbstractMatrix)
    _validate_two_qubit_unitary(U)

    φ00 = _diagonal_phase(U[1, 1], 1)
    φ01 = _diagonal_phase(U[2, 2], 2)
    φ10 = _diagonal_phase(U[3, 3], 3)
    correction_phases = [φ00, φ01, φ10, φ01 + φ10 - φ00]
    correction = Diagonal(exp.(-1im .* correction_phases))

    return Matrix{ComplexF64}(U * correction)
end

function conditional_phase(U::AbstractMatrix)
    _validate_two_qubit_unitary(U)

    φ00 = _diagonal_phase(U[1, 1], 1)
    φ01 = _diagonal_phase(U[2, 2], 2)
    φ10 = _diagonal_phase(U[3, 3], 3)
    φ11 = _diagonal_phase(U[4, 4], 4)
    return _wrap_phase(φ00 + φ11 - φ01 - φ10)
end

function compare_spectra(reference::SpectrumResult, candidate::SpectrumResult)
    level_count = min(length(reference.energies), length(candidate.energies))
    level_count > 0 || throw(ArgumentError("Spectrum comparison requires at least one level in each result."))

    reference_levels = reference.energies[1:level_count]
    candidate_levels = candidate.energies[1:level_count]
    level_deltas = candidate_levels .- reference_levels
    transition_deltas = diff(candidate_levels) .- diff(reference_levels)

    return SpectrumComparisonResult(reference, candidate, level_deltas, transition_deltas)
end

function compare_model_spectra(
    system::CompositeSystem;
    reference_hamiltonian_spec::AbstractHamiltonianSpec,
    candidate_hamiltonian_spec::AbstractHamiltonianSpec = EffectiveHamiltonianSpec(),
    levels::Integer = 6,
)
    return compare_model_spectra(
        system,
        system;
        reference_hamiltonian_spec = reference_hamiltonian_spec,
        candidate_hamiltonian_spec = candidate_hamiltonian_spec,
        levels = levels,
    )
end

function compare_model_spectra(
    reference_system::CompositeSystem,
    candidate_system::CompositeSystem;
    reference_hamiltonian_spec::AbstractHamiltonianSpec,
    candidate_hamiltonian_spec::AbstractHamiltonianSpec = EffectiveHamiltonianSpec(),
    levels::Integer = 6,
)
    reference_result = spectrum(reference_system; levels = levels, hamiltonian_spec = reference_hamiltonian_spec)
    candidate_result = spectrum(candidate_system; levels = levels, hamiltonian_spec = candidate_hamiltonian_spec)
    return compare_spectra(reference_result, candidate_result)
end

function compare_minimum_gap(reference::SweepResult, candidate::SweepResult; level_pair::Tuple{<:Integer, <:Integer} = (1, 2))
    reference_gap = minimum_gap(reference; level_pair = level_pair)
    candidate_gap = minimum_gap(candidate; level_pair = level_pair)
    return (
        reference = reference_gap,
        candidate = candidate_gap,
        gap_delta = candidate_gap.gap - reference_gap.gap,
        sweep_value_delta = candidate_gap.sweep_value - reference_gap.sweep_value,
    )
end

function best_move(result::UnitaryTraceResult; source_index::Integer, target_index::Integer)
    subspace_dimension = size(first(result.unitaries), 1)
    1 <= source_index <= subspace_dimension || throw(ArgumentError("source_index must satisfy 1 <= source_index <= $subspace_dimension."))
    1 <= target_index <= subspace_dimension || throw(ArgumentError("target_index must satisfy 1 <= target_index <= $subspace_dimension."))

    probabilities = [abs2(U[target_index, source_index]) for U in result.unitaries]
    best_index = argmax(probabilities)

    return (
        index = best_index,
        time = result.times[best_index],
        transfer_probability = probabilities[best_index],
        leakage = result.leakages[source_index, best_index],
        unitary = result.unitaries[best_index],
    )
end

function best_cz(result::UnitaryTraceResult; target_phase::Real = π)
    size(first(result.unitaries)) == (4, 4) ||
        throw(ArgumentError("best_cz requires a four-state subspace."))

    phase_errors = Float64[]
    max_leakages = Float64[]

    for time_index in eachindex(result.unitaries)
        push!(phase_errors, abs(_wrap_phase(conditional_phase(result.unitaries[time_index]) - Float64(target_phase))))
        push!(max_leakages, maximum(view(result.leakages, :, time_index)))
    end

    best_index = first(eachindex(result.unitaries))
    best_metric = phase_errors[best_index] + max_leakages[best_index]

    for time_index in eachindex(result.unitaries)
        metric = phase_errors[time_index] + max_leakages[time_index]
        if metric < best_metric
            best_metric = metric
            best_index = time_index
        end
    end

    return (
        index = best_index,
        time = result.times[best_index],
        conditional_phase = conditional_phase(result.unitaries[best_index]),
        phase_error = phase_errors[best_index],
        leakage = max_leakages[best_index],
        unitary = result.unitaries[best_index],
    )
end

function _normalize_subsystem_levels(model::StaticSystemModel, subsystem_levels)
    if subsystem_levels isa NamedTuple
        targets = collect(keys(subsystem_levels))
        level_sets = [_normalize_level_set(model, target, getfield(subsystem_levels, target)) for target in targets]
        return _validated_subsystem_level_spec(targets, level_sets)
    elseif subsystem_levels isa AbstractVector{<:Pair}
        targets = Symbol[]
        level_sets = Vector{Vector{Int}}()

        for entry in subsystem_levels
            target = Symbol(first(entry))
            push!(targets, target)
            push!(level_sets, _normalize_level_set(model, target, last(entry)))
        end

        return _validated_subsystem_level_spec(targets, level_sets)
    elseif subsystem_levels isa AbstractDict
        for target in keys(subsystem_levels)
            Symbol(target) in model.subsystem_order ||
                throw(ArgumentError("No subsystem named $(Symbol(target)) exists in the model."))
        end
        targets = [target for target in model.subsystem_order if haskey(subsystem_levels, target)]
        isempty(targets) && throw(ArgumentError("subsystem_levels did not match any subsystem in the model."))
        level_sets = [_normalize_level_set(model, target, subsystem_levels[target]) for target in targets]
        return _validated_subsystem_level_spec(targets, level_sets)
    end

    throw(ArgumentError("subsystem_levels must be a NamedTuple, a vector of pairs, or a dictionary."))
end

function _normalize_level_set(model::StaticSystemModel, target::Symbol, levels)
    target in model.subsystem_order || throw(ArgumentError("No subsystem named $target exists in the model."))

    raw_levels =
        levels isa Integer ? [Int(levels)] :
        levels isa AbstractVector ? Int[level for level in levels] :
        levels isa AbstractRange ? collect(Int, levels) :
        throw(ArgumentError("Levels for subsystem $target must be an integer, vector, or range."))

    isempty(raw_levels) && throw(ArgumentError("Levels for subsystem $target must not be empty."))
    length(unique(raw_levels)) == length(raw_levels) ||
        throw(ArgumentError("Levels for subsystem $target must be unique."))

    dimension = model.dimensions[target]
    for level in raw_levels
        0 <= level < dimension ||
            throw(ArgumentError("Levels for subsystem $target must satisfy 0 <= level < $dimension."))
    end

    return sort(raw_levels)
end

function _enumerate_assignments(targets::Vector{Symbol}, level_sets::Vector{Vector{Int}})
    assignments = Dict{Symbol, Int}[]
    current = Dict{Symbol, Int}()

    function _walk(index::Int)
        if index > length(targets)
            push!(assignments, copy(current))
            return nothing
        end

        target = targets[index]
        for level in level_sets[index]
            current[target] = level
            _walk(index + 1)
        end

        delete!(current, target)
        return nothing
    end

    _walk(1)
    return assignments
end

function _assignment_label(targets::Vector{Symbol}, assignment::Dict{Symbol, Int})
    levels = [assignment[target] for target in targets]
    if all(level -> 0 <= level <= 9, levels)
        return "|" * join(string.(levels)) * ">"
    end

    return "|" * join(["$(targets[index])=$(levels[index])" for index in eachindex(targets)], ",") * ">"
end

function _assignment_state(model::StaticSystemModel, assignment::Dict{Symbol, Int})
    kwargs = (; assignment...)
    return basis_state(model; kwargs...)
end

function _bare_reference_bundle(bare_states::Vector)
    reference_energies = fill(NaN, length(bare_states))
    reference_overlaps = ones(Float64, length(bare_states))
    return bare_states, reference_energies, reference_overlaps
end

function _dressed_reference_bundle(model::StaticSystemModel, bare_states::Vector)
    eigensolve = eigenstates(hamiltonian(model))
    eigenvectors = eigensolve.vectors
    eigenvalues = Float64.(real.(eigensolve.values))
    dims = first(bare_states).dims
    eigenstates_kets = [QuantumObject(eigenvectors[:, index]; type = Ket(), dims = dims) for index in axes(eigenvectors, 2)]

    overlap_matrix = Matrix{Float64}(undef, length(bare_states), length(eigenstates_kets))
    for row in eachindex(bare_states), column in eachindex(eigenstates_kets)
        overlap_matrix[row, column] = abs2(dag(bare_states[row]) * eigenstates_kets[column])
    end

    matched_states = similar(bare_states)
    matched_energies = zeros(Float64, length(bare_states))
    matched_overlaps = zeros(Float64, length(bare_states))
    assigned_rows = falses(length(bare_states))
    assigned_columns = falses(length(eigenstates_kets))

    for _ in eachindex(bare_states)
        best_row = 0
        best_column = 0
        best_overlap = -1.0

        for row in eachindex(bare_states)
            assigned_rows[row] && continue
            for column in eachindex(eigenstates_kets)
                assigned_columns[column] && continue
                overlap = overlap_matrix[row, column]
                if overlap > best_overlap
                    best_overlap = overlap
                    best_row = row
                    best_column = column
                end
            end
        end

        best_row > 0 || throw(ArgumentError("Failed to match dressed states to the requested bare subspace."))
        matched_states[best_row] = eigenstates_kets[best_column]
        matched_energies[best_row] = eigenvalues[best_column]
        matched_overlaps[best_row] = best_overlap
        assigned_rows[best_row] = true
        assigned_columns[best_column] = true
    end

    return matched_states, matched_energies, matched_overlaps
end

function _validate_subspace_spec(model::StaticSystemModel, spec::SubspaceSpec)
    spec.model_subsystem_order == model.subsystem_order ||
        throw(ArgumentError("SubspaceSpec was built for a different subsystem order than the provided model."))
    spec.model_dimensions == model.dimensions ||
        throw(ArgumentError("SubspaceSpec dimensions do not match the provided model."))
    return nothing
end

function _validate_two_qubit_unitary(U::AbstractMatrix)
    size(U, 1) == size(U, 2) || throw(ArgumentError("Unitary must be square."))
    size(U, 1) == 4 || throw(ArgumentError("Two-qubit phase utilities require a 4x4 matrix."))
    return nothing
end

function _diagonal_phase(value, index::Integer)
    abs(value) > 1e-12 || throw(ArgumentError("Diagonal element $index is too small to define a stable phase."))
    return angle(value)
end

_wrap_phase(value::Real) = mod(value + π, 2π) - π

function _validated_subsystem_level_spec(targets::Vector{Symbol}, level_sets::Vector{Vector{Int}})
    isempty(targets) && throw(ArgumentError("subsystem_levels must select at least one subsystem."))
    length(unique(targets)) == length(targets) ||
        throw(ArgumentError("subsystem_levels must not contain duplicate subsystem names."))
    return targets, level_sets
end
