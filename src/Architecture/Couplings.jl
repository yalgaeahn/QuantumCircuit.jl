abstract type AbstractCoupling end

function coupling_endpoints(coupling::AbstractCoupling)
    throw(ArgumentError("coupling_endpoints is not implemented for $(typeof(coupling))."))
end

struct CapacitiveCoupling <: AbstractCoupling
    source::Symbol
    target::Symbol
    g::Float64

    function CapacitiveCoupling(source::Symbol, target::Symbol; g::Real)
        validate_identifier(source, "coupling source")
        validate_identifier(target, "coupling target")
        source != target || throw(ArgumentError("Coupling endpoints must be distinct subsystem names."))

        new(source, target, validate_positive(g, "g"))
    end
end

coupling_endpoints(coupling::CapacitiveCoupling) = (coupling.source, coupling.target)

struct CircuitCapacitiveCoupling <: AbstractCoupling
    source::Symbol
    target::Symbol
    G::Float64

    function CircuitCapacitiveCoupling(source::Symbol, target::Symbol; G::Real)
        validate_identifier(source, "coupling source")
        validate_identifier(target, "coupling target")
        source != target || throw(ArgumentError("Coupling endpoints must be distinct subsystem names."))

        new(source, target, validate_positive(G, "G"))
    end
end

coupling_endpoints(coupling::CircuitCapacitiveCoupling) = (coupling.source, coupling.target)
