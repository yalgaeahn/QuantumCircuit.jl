function basis_state(model::StaticSystemModel; assignments...)
    requested_levels = Dict{Symbol, Int}()

    for (subsystem_name, level) in pairs(assignments)
        subsystem_name in model.subsystem_order ||
            throw(ArgumentError("No subsystem named $subsystem_name exists in the model."))
        level isa Integer ||
            throw(ArgumentError("Basis-state level for subsystem $subsystem_name must be an integer."))

        dimension = model.dimensions[subsystem_name]
        0 <= level < dimension ||
            throw(
                ArgumentError(
                    "Basis-state level for subsystem $subsystem_name must satisfy 0 <= level < $dimension.",
                ),
            )

        requested_levels[subsystem_name] = Int(level)
    end

    local_states = [basis(model.dimensions[subsystem_name], get(requested_levels, subsystem_name, 0)) for subsystem_name in model.subsystem_order]
    return length(local_states) == 1 ? only(local_states) : kron(local_states...)
end

basis_state(system::CompositeSystem; assignments...) = basis_state(build_model(system); assignments...)
