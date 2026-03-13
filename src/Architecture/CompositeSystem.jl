struct CompositeSystem
    subsystems::Vector{AbstractSubsystem}
    couplings::Vector{AbstractCoupling}

    function CompositeSystem(subsystems::Vector{AbstractSubsystem}, couplings::Vector{AbstractCoupling})
        isempty(subsystems) && throw(ArgumentError("CompositeSystem requires at least one subsystem."))

        names = Symbol[name(subsystem) for subsystem in subsystems]
        length(unique(names)) == length(names) ||
            throw(ArgumentError("Subsystem names in CompositeSystem must be unique."))

        valid_names = Set(names)
        for coupling in couplings
            source, target = coupling_endpoints(coupling)
            source in valid_names ||
                throw(ArgumentError("Coupling source $source does not match any subsystem in CompositeSystem."))
            target in valid_names ||
                throw(ArgumentError("Coupling target $target does not match any subsystem in CompositeSystem."))
        end

        new(subsystems, couplings)
    end
end

function CompositeSystem(
    subsystem_items::AbstractVector{<:AbstractSubsystem},
    coupling_items::AbstractVector{<:AbstractCoupling} = AbstractCoupling[],
)
    return CompositeSystem(AbstractSubsystem[subsystem_items...], AbstractCoupling[coupling_items...])
end

function CompositeSystem(items...)
    subsystem_items = AbstractSubsystem[]
    coupling_items = AbstractCoupling[]

    for item in items
        if item isa AbstractSubsystem
            push!(subsystem_items, item)
        elseif item isa AbstractCoupling
            push!(coupling_items, item)
        else
            throw(
                ArgumentError(
                    "CompositeSystem only accepts subsystems and couplings, got $(typeof(item)).",
                ),
            )
        end
    end

    return CompositeSystem(subsystem_items, coupling_items)
end

subsystems(system::CompositeSystem) = system.subsystems
couplings(system::CompositeSystem) = system.couplings
subsystem_names(system::CompositeSystem) = Symbol[name(subsystem) for subsystem in system.subsystems]

function with_subsystem_parameter(system::CompositeSystem, subsystem_name::Symbol, parameter::Symbol, value)
    updated_subsystems = AbstractSubsystem[]
    found = false

    for subsystem in subsystems(system)
        if name(subsystem) == subsystem_name
            push!(updated_subsystems, _with_parameter(subsystem, parameter, value))
            found = true
        else
            push!(updated_subsystems, subsystem)
        end
    end

    found || throw(ArgumentError("No subsystem named $subsystem_name exists in CompositeSystem."))
    return CompositeSystem(updated_subsystems, couplings(system))
end

function with_coupling_parameter(system::CompositeSystem, source::Symbol, target::Symbol, parameter::Symbol, value)
    updated_couplings = AbstractCoupling[]
    match_count = 0

    for coupling in couplings(system)
        if _coupling_matches(coupling, source, target)
            push!(updated_couplings, _with_parameter(coupling, parameter, value))
            match_count += 1
        else
            push!(updated_couplings, coupling)
        end
    end

    match_count == 1 ||
        throw(
            ArgumentError(
                "Expected exactly one coupling between $source and $target, found $match_count.",
            ),
        )

    return CompositeSystem(subsystems(system), updated_couplings)
end

function _with_parameter(subsystem::AbstractSubsystem, parameter::Symbol, value)
    fields = fieldnames(typeof(subsystem))
    parameter in fields || throw(ArgumentError("$(typeof(subsystem)) has no field named $parameter."))
    parameter != :name || throw(ArgumentError("Subsystem names cannot be changed through sweep helpers."))

    kwargs = (; (field => (field == parameter ? value : getfield(subsystem, field)) for field in fields if field != :name)...)
    return typeof(subsystem)(name(subsystem); kwargs...)
end

function _with_parameter(coupling::AbstractCoupling, parameter::Symbol, value)
    fields = fieldnames(typeof(coupling))
    parameter in fields || throw(ArgumentError("$(typeof(coupling)) has no field named $parameter."))
    (parameter != :source && parameter != :target) ||
        throw(ArgumentError("Coupling endpoints cannot be changed through sweep helpers."))

    source, target = coupling_endpoints(coupling)
    kwargs = (; (field => (field == parameter ? value : getfield(coupling, field)) for field in fields if field != :source && field != :target)...)
    return typeof(coupling)(source, target; kwargs...)
end

function _coupling_matches(coupling::AbstractCoupling, source::Symbol, target::Symbol)
    coupling_source, coupling_target = coupling_endpoints(coupling)
    return (coupling_source == source && coupling_target == target) ||
        (coupling_source == target && coupling_target == source)
end
