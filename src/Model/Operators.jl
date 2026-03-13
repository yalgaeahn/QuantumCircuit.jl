const _SUPPORTED_LOCAL_OPERATOR_SYMBOLS = (:a, :adag, :n, :x, :y)

struct ObservableSpec
    label::Symbol
    target::Symbol
    operator::Symbol

    function ObservableSpec(label::Symbol, target::Symbol, operator::Symbol)
        _validate_nonempty_symbol(label, "observable label")
        _validate_nonempty_symbol(target, "observable target")
        _validate_operator_symbol(operator, "observable operator")
        return new(label, target, operator)
    end
end

function _validate_nonempty_symbol(value::Symbol, label::AbstractString)
    value == Symbol("") && throw(ArgumentError("$label must not be empty."))
    return value
end

function _validate_operator_symbol(operator::Symbol, label::AbstractString)
    operator in _SUPPORTED_LOCAL_OPERATOR_SYMBOLS ||
        throw(ArgumentError("$label must be one of $(_SUPPORTED_LOCAL_OPERATOR_SYMBOLS), got $operator."))
    return operator
end

function annihilation_operator(model::StaticSystemModel, target::Symbol)
    return model.annihilation_operators[_subsystem_name(model, target)]
end

annihilation_operator(system::CompositeSystem, target::Symbol) = annihilation_operator(build_model(system), target)

creation_operator(model::StaticSystemModel, target::Symbol) = dag(annihilation_operator(model, target))
creation_operator(system::CompositeSystem, target::Symbol) = creation_operator(build_model(system), target)

number_operator(model::StaticSystemModel, target::Symbol) = creation_operator(model, target) * annihilation_operator(model, target)
number_operator(system::CompositeSystem, target::Symbol) = number_operator(build_model(system), target)

function quadrature_operator(model::StaticSystemModel, target::Symbol, axis::Symbol)
    a = annihilation_operator(model, target)
    adag = dag(a)
    axis == :x && return a + adag
    axis == :y && return -1im * (a - adag)
    throw(ArgumentError("quadrature axis must be :x or :y, got $axis."))
end

quadrature_operator(system::CompositeSystem, target::Symbol, axis::Symbol) =
    quadrature_operator(build_model(system), target, axis)

function _embedded_operator(model::StaticSystemModel, target::Symbol, operator::Symbol)
    operator == :a && return annihilation_operator(model, target)
    operator == :adag && return creation_operator(model, target)
    operator == :n && return number_operator(model, target)
    (operator == :x || operator == :y) && return quadrature_operator(model, target, operator)
    _validate_operator_symbol(operator, "operator")
    error("Unreachable operator symbol branch: $operator")
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
