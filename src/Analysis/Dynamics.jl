using QuantumToolbox: expect
using ..Model: _projector_operator

function observable_trace(result::DynamicsResult, label::Symbol)
    for trace in result.observables
        trace.label == label && return trace
    end

    throw(ArgumentError("No observable trace labeled $label exists in the dynamics result."))
end

function final_state(result::DynamicsResult)
    isempty(result.states) && throw(ArgumentError("DynamicsResult does not contain any saved states."))
    return result.states[end]
end

function population_trace(result::DynamicsResult, subsystem::Symbol, level::Int)
    isempty(result.states) && throw(ArgumentError("DynamicsResult does not contain any saved states."))

    projector = _projector_operator(result.model, subsystem, level)
    values = Float64[real(expect(projector, state)) for state in result.states]
    label = Symbol("population_", subsystem, "_", level)

    return ObservableTrace(label, result.times, values)
end
