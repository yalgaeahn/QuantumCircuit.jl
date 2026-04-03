struct BareSubsystemBasis
    name::Symbol
    energies::Vector{Float64}
    states::Vector
end

struct BareProductBasis
    subsystems::Vector{BareSubsystemBasis}
    labels::Vector{Tuple}
    states::Vector
end

struct LabelMap
    bare_to_dressed::Vector{Union{Int, Missing}}
    confidences::Vector{Float64}
end

struct LabeledEigensystemSweepResult{Spec, T, Result}
    sweep::EigensystemSweepResult{Spec, T, Result}
    basis::BareProductBasis
    label_maps::Vector{LabelMap}
end

function bare_product_basis(
    system::CompositeSystem;
    hamiltonian_spec::AbstractHamiltonianSpec = EffectiveHamiltonianSpec(),
    subsystem_levels = nothing,
)
    subsystem_list = subsystems(system)
    isempty(subsystem_list) && throw(ArgumentError("bare_product_basis requires at least one subsystem."))
    level_map = _normalize_subsystem_levels(subsystem_levels, subsystem_list)
    subsystem_bases = BareSubsystemBasis[]

    for subsystem in subsystem_list
        subsystem_name = name(subsystem)
        local_model = build_model(CompositeSystem(subsystem); hamiltonian_spec = hamiltonian_spec)
        max_levels = local_model.dimensions[subsystem_name]
        level_count = get(level_map, subsystem_name, max_levels)
        level_count <= max_levels || throw(
            ArgumentError(
                "Requested $level_count levels for subsystem $subsystem_name, but only $max_levels are available.",
            ),
        )
        local_esys = eigensystem(local_model; levels = level_count)
        push!(subsystem_bases, BareSubsystemBasis(subsystem_name, local_esys.energies, collect(local_esys.states)))
    end

    state_products = Vector{Any}()
    label_products = Tuple[]
    index_ranges = (0:(length(subsystem_basis.states) - 1) for subsystem_basis in subsystem_bases)

    for label in Iterators.product(index_ranges...)
        bare_label = Tuple(Int(level) for level in label)
        local_states = [subsystem_bases[index].states[bare_label[index] + 1] for index in eachindex(subsystem_bases)]
        push!(label_products, bare_label)
        push!(state_products, length(local_states) == 1 ? only(local_states) : kron(local_states...))
    end

    return BareProductBasis(subsystem_bases, label_products, state_products)
end

function label_sweep(
    sweep::EigensystemSweepResult,
    basis::BareProductBasis;
    overlap_threshold::Real = 0.5,
)
    0.0 <= overlap_threshold <= 1.0 ||
        throw(ArgumentError("overlap_threshold must satisfy 0 <= overlap_threshold <= 1."))
    label_maps = [_label_eigensystem(esys, basis; overlap_threshold = Float64(overlap_threshold)) for esys in sweep.spectra]
    return LabeledEigensystemSweepResult(sweep, basis, label_maps)
end

function dressed_index(labeled::LabeledEigensystemSweepResult, bare_label, sweep_index::Integer)
    map = _label_map_at(labeled, sweep_index)
    position = _lookup_position(labeled.basis, bare_label)
    return map.bare_to_dressed[position]
end

function bare_label(labeled::LabeledEigensystemSweepResult, dressed_state_index::Integer, sweep_index::Integer)
    map = _label_map_at(labeled, sweep_index)
    dressed_state_index > 0 || throw(ArgumentError("dressed_state_index must be positive."))
    position = findfirst(==(dressed_state_index), map.bare_to_dressed)
    return isnothing(position) ? missing : labeled.basis.labels[position]
end

function energy_by_bare_label(
    labeled::LabeledEigensystemSweepResult,
    bare_label,
    sweep_index::Integer;
    subtract_ground::Bool = false,
)
    assigned_dressed_index = dressed_index(labeled, bare_label, sweep_index)
    ismissing(assigned_dressed_index) && return NaN
    energies = labeled.sweep.spectra[sweep_index].energies
    energy = energies[assigned_dressed_index]
    return subtract_ground ? energy - energies[1] : energy
end

function energy_curve(
    labeled::LabeledEigensystemSweepResult,
    bare_label;
    subtract_ground::Bool = false,
)
    normalized_label = _normalize_bare_label(labeled.basis, bare_label)
    data = [
        energy_by_bare_label(labeled, normalized_label, sweep_index; subtract_ground = subtract_ground) for
        sweep_index in eachindex(labeled.sweep.spectra)
    ]
    return SweepSeries(
        copy(labeled.sweep.values),
        data,
        Symbol("bare_" * join(string.(normalized_label), "_")),
    )
end

function dressed_state_components(
    labeled::LabeledEigensystemSweepResult,
    state_label,
    sweep_index::Integer;
    top::Union{Nothing, Integer} = nothing,
)
    top === nothing || top > 0 || throw(ArgumentError("top must be a positive integer when provided."))
    eigensystem_result = _eigensystem_at(labeled, sweep_index)
    dressed_state_index = if state_label isa Integer
        state_label
    else
        assigned_dressed_index = dressed_index(labeled, state_label, sweep_index)
        ismissing(assigned_dressed_index) && throw(
            ArgumentError("No dressed state is assigned to bare label $state_label at sweep index $sweep_index."),
        )
        assigned_dressed_index
    end
    1 <= dressed_state_index <= length(eigensystem_result.states) ||
        throw(ArgumentError("dressed_state_index $dressed_state_index is out of bounds for the available eigensystem."))

    dressed_state = eigensystem_result.states[dressed_state_index]
    components = Pair{Tuple, Float64}[]
    for (position, bare_state) in enumerate(labeled.basis.states)
        push!(components, labeled.basis.labels[position] => abs2(dag(bare_state) * dressed_state))
    end
    sort!(components; by = last, rev = true)
    return isnothing(top) ? components : components[1:min(top, length(components))]
end

function _label_eigensystem(
    eigensystem_result::EigensystemResult,
    basis::BareProductBasis;
    overlap_threshold::Float64,
)
    bare_to_dressed = Union{Int, Missing}[missing for _ in eachindex(basis.labels)]
    confidences = zeros(Float64, length(basis.labels))
    unclaimed = trues(length(basis.labels))

    for (dressed_state_index, dressed_state) in enumerate(eigensystem_result.states)
        best_position = 0
        best_overlap = -1.0
        for position in eachindex(basis.states)
            unclaimed[position] || continue
            overlap = abs2(dag(basis.states[position]) * dressed_state)
            if overlap > best_overlap
                best_overlap = overlap
                best_position = position
            end
        end
        best_position == 0 && continue
        best_overlap >= overlap_threshold || continue
        bare_to_dressed[best_position] = dressed_state_index
        confidences[best_position] = best_overlap
        unclaimed[best_position] = false
    end

    return LabelMap(bare_to_dressed, confidences)
end

function _normalize_subsystem_levels(subsystem_levels, subsystem_list)
    subsystem_names = Set(name(subsystem) for subsystem in subsystem_list)
    subsystem_levels === nothing && return Dict{Symbol, Int}()
    subsystem_levels isa AbstractDict ||
        throw(ArgumentError("subsystem_levels must be an AbstractDict keyed by subsystem Symbol names."))

    normalized = Dict{Symbol, Int}()
    for (subsystem_name, level_count) in pairs(subsystem_levels)
        subsystem_name isa Symbol ||
            throw(ArgumentError("subsystem_levels keys must be Symbols, got $(typeof(subsystem_name))."))
        subsystem_name in subsystem_names ||
            throw(ArgumentError("subsystem_levels references missing subsystem $subsystem_name."))
        level_count isa Integer ||
            throw(ArgumentError("subsystem_levels[$subsystem_name] must be an integer."))
        level_count > 0 ||
            throw(ArgumentError("subsystem_levels[$subsystem_name] must be positive."))
        normalized[subsystem_name] = Int(level_count)
    end

    return normalized
end

function _normalize_bare_label(basis::BareProductBasis, bare_label)
    bare_label isa Tuple || throw(ArgumentError("bare_label must be provided as a tuple of subsystem excitation levels."))
    expected_length = length(basis.subsystems)
    length(bare_label) == expected_length ||
        throw(ArgumentError("bare_label must have length $expected_length, got $(length(bare_label))."))
    normalized = Tuple(Int(level) for level in bare_label)
    all(level isa Integer for level in normalized) ||
        throw(ArgumentError("bare_label must contain only integer excitation levels."))
    return normalized
end

function _lookup_position(basis::BareProductBasis, bare_label)
    normalized_label = _normalize_bare_label(basis, bare_label)
    position = findfirst(==(normalized_label), basis.labels)
    isnothing(position) && throw(ArgumentError("bare label $normalized_label is not present in the bare product basis."))
    return position
end

function _label_map_at(labeled::LabeledEigensystemSweepResult, sweep_index::Integer)
    1 <= sweep_index <= length(labeled.label_maps) ||
        throw(ArgumentError("sweep_index $sweep_index is out of bounds for the labeled sweep."))
    return labeled.label_maps[sweep_index]
end

function _eigensystem_at(labeled::LabeledEigensystemSweepResult, sweep_index::Integer)
    1 <= sweep_index <= length(labeled.sweep.spectra) ||
        throw(ArgumentError("sweep_index $sweep_index is out of bounds for the labeled sweep."))
    return labeled.sweep.spectra[sweep_index]
end
