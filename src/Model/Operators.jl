struct ObservableSpec
    label::Symbol
    target::Symbol
    operator::Symbol

    function ObservableSpec(label::Symbol, target::Symbol, operator::Symbol)
        _validate_nonempty_symbol(label, "observable label")
        _validate_nonempty_symbol(target, "observable target")
        _validate_nonempty_symbol(operator, "observable operator")
        return new(label, target, operator)
    end
end

function annihilation_operator(model::StaticSystemModel, target::Symbol)
    return _lookup_operator(model, target, :a)
end

function annihilation_operator(
    system::CompositeSystem,
    target::Symbol;
    hamiltonian_spec::AbstractHamiltonianSpec = EffectiveHamiltonianSpec(),
)
    return annihilation_operator(build_model(system; hamiltonian_spec = hamiltonian_spec), target)
end

creation_operator(model::StaticSystemModel, target::Symbol) = _lookup_operator(model, target, :adag)

function creation_operator(
    system::CompositeSystem,
    target::Symbol;
    hamiltonian_spec::AbstractHamiltonianSpec = EffectiveHamiltonianSpec(),
)
    return creation_operator(build_model(system; hamiltonian_spec = hamiltonian_spec), target)
end

number_operator(model::StaticSystemModel, target::Symbol) = _lookup_operator(model, target, :n)

function number_operator(
    system::CompositeSystem,
    target::Symbol;
    hamiltonian_spec::AbstractHamiltonianSpec = EffectiveHamiltonianSpec(),
)
    return number_operator(build_model(system; hamiltonian_spec = hamiltonian_spec), target)
end

function quadrature_operator(model::StaticSystemModel, target::Symbol, axis::Symbol)
    (axis == :x || axis == :y) || throw(ArgumentError("quadrature axis must be :x or :y, got $axis."))
    return _lookup_operator(model, target, axis)
end

function quadrature_operator(
    system::CompositeSystem,
    target::Symbol,
    axis::Symbol;
    hamiltonian_spec::AbstractHamiltonianSpec = EffectiveHamiltonianSpec(),
)
    return quadrature_operator(build_model(system; hamiltonian_spec = hamiltonian_spec), target, axis)
end

charge_operator(model::StaticSystemModel, target::Symbol) = _lookup_operator(model, target, :charge)

function charge_operator(
    system::CompositeSystem,
    target::Symbol;
    hamiltonian_spec::AbstractHamiltonianSpec = EffectiveHamiltonianSpec(),
)
    return charge_operator(build_model(system; hamiltonian_spec = hamiltonian_spec), target)
end

cosphi_operator(model::StaticSystemModel, target::Symbol) = _lookup_operator(model, target, :cosphi)

function cosphi_operator(
    system::CompositeSystem,
    target::Symbol;
    hamiltonian_spec::AbstractHamiltonianSpec = EffectiveHamiltonianSpec(),
)
    return cosphi_operator(build_model(system; hamiltonian_spec = hamiltonian_spec), target)
end

sinphi_operator(model::StaticSystemModel, target::Symbol) = _lookup_operator(model, target, :sinphi)

function sinphi_operator(
    system::CompositeSystem,
    target::Symbol;
    hamiltonian_spec::AbstractHamiltonianSpec = EffectiveHamiltonianSpec(),
)
    return sinphi_operator(build_model(system; hamiltonian_spec = hamiltonian_spec), target)
end

function _embedded_operator(model::StaticSystemModel, target::Symbol, operator::Symbol)
    return _lookup_operator(model, target, operator)
end

function _lookup_operator(model::StaticSystemModel, target::Symbol, operator::Symbol)
    subsystem_name = _subsystem_name(model, target)
    cache_key = (subsystem_name, operator)

    if haskey(model.operator_cache, cache_key)
        return model.operator_cache[cache_key]
    end

    supported_operators = _supported_operator_symbols(model, subsystem_name)
    throw(
        ArgumentError(
            "Operator $operator is not supported for subsystem $subsystem_name under $(nameof(typeof(model.hamiltonian_spec))). Supported operators: $(supported_operators).",
        ),
    )
end

function _supported_operator_symbols(model::StaticSystemModel, target::Symbol)
    symbols = Symbol[key[2] for key in keys(model.operator_cache) if key[1] == target]
    sort!(symbols; by = string)
    return Tuple(symbols)
end

function _projector_operator(model::StaticSystemModel, target::Symbol, level::Int)
    subsystem_name = _subsystem_name(model, target)
    dimension = model.dimensions[subsystem_name]
    0 <= level < dimension ||
        throw(ArgumentError("Population level for subsystem $subsystem_name must satisfy 0 <= level < $dimension."))

    ket = basis(dimension, level)
    return _embed_operator(ket * dag(ket), _subsystem_dimensions(model), _subsystem_index(model, subsystem_name))
end

function _subsystem_name(model::StaticSystemModel, target::Symbol)
    target in keys(model.dimensions) || throw(ArgumentError("No subsystem named $target exists in the model."))
    return target
end

function _subsystem_index(model::StaticSystemModel, target::Symbol)
    index = findfirst(==(target), model.subsystem_order)
    isnothing(index) && throw(ArgumentError("No subsystem named $target exists in the model."))
    return index
end

_subsystem_dimensions(model::StaticSystemModel) = [model.dimensions[subsystem_name] for subsystem_name in model.subsystem_order]
