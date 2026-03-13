struct SubsystemDrive{F}
    label::Symbol
    target::Symbol
    operator::Symbol
    coefficient::F

    function SubsystemDrive(label::Symbol, target::Symbol, operator::Symbol, coefficient::F) where {F}
        _validate_nonempty_symbol(label, "drive label")
        _validate_nonempty_symbol(target, "drive target")
        _validate_operator_symbol(operator, "drive operator")
        coefficient isa Function || throw(ArgumentError("drive coefficient must be callable."))
        return new{F}(label, target, operator, coefficient)
    end
end

function _time_dependent_hamiltonian(model::StaticSystemModel, drives::Nothing)
    return hamiltonian(model)
end

function _time_dependent_hamiltonian(model::StaticSystemModel, drives::Union{Tuple, AbstractVector})
    isempty(drives) && return hamiltonian(model)

    terms = Any[hamiltonian(model)]
    for drive in drives
        drive isa SubsystemDrive ||
            throw(ArgumentError("drives must contain only SubsystemDrive values, got $(typeof(drive))."))
        operator = _embedded_operator(model, drive.target, drive.operator)
        push!(terms, (operator, (params, t) -> drive.coefficient(params, t)))
    end

    return Tuple(terms)
end
