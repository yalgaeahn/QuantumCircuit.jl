module Model

using QuantumToolbox: basis, dag, destroy, eye, kron
using ..Architecture:
    AbstractCoupling,
    AbstractSubsystem,
    CapacitiveCoupling,
    CompositeSystem,
    Resonator,
    TunableCoupler,
    TunableTransmon,
    Transmon,
    coupling_endpoints,
    couplings,
    name,
    subsystems

export ObservableSpec,
    StaticSystemModel,
    SubsystemDrive,
    annihilation_operator,
    basis_state,
    build_model,
    creation_operator,
    hamiltonian,
    number_operator,
    quadrature_operator

struct StaticSystemModel{QT}
    system::CompositeSystem
    subsystem_order::Vector{Symbol}
    dimensions::Dict{Symbol, Int}
    hamiltonian::QT
    annihilation_operators::Dict{Symbol, Any}
end

hamiltonian(model::StaticSystemModel) = model.hamiltonian

function build_model(system::CompositeSystem)
    subsystem_list = subsystems(system)
    dims = [_local_dimension(subsystem) for subsystem in subsystem_list]
    order = Symbol[name(subsystem) for subsystem in subsystem_list]

    annihilation_ops = Dict{Symbol, Any}()
    H = _zero_operator(dims)

    for (index, subsystem) in enumerate(subsystem_list)
        local_annihilator = destroy(dims[index])
        annihilation_ops[name(subsystem)] = _embed_operator(local_annihilator, dims, index)
        H += _embedded_local_hamiltonian(subsystem, dims, index, local_annihilator)
    end

    for coupling in couplings(system)
        H += _coupling_hamiltonian(coupling, annihilation_ops)
    end

    return StaticSystemModel(
        system,
        order,
        Dict(order[i] => dims[i] for i in eachindex(order)),
        H,
        annihilation_ops,
    )
end

hamiltonian(system::CompositeSystem) = hamiltonian(build_model(system))

_local_dimension(subsystem::Transmon) = subsystem.ncut
_local_dimension(subsystem::Resonator) = subsystem.dim
_local_dimension(subsystem::TunableTransmon) = subsystem.ncut
_local_dimension(subsystem::TunableCoupler) = subsystem.ncut
_local_dimension(subsystem::AbstractSubsystem) = throw(ArgumentError("Unsupported subsystem type $(typeof(subsystem))."))

function _zero_operator(dims::Vector{Int})
    identity_factors = [eye(dim) for dim in dims]
    return 0.0 * kron(identity_factors...)
end

function _embed_operator(local_operator, dims::Vector{Int}, index::Int)
    operators = [position == index ? local_operator : eye(dims[position]) for position in eachindex(dims)]
    return kron(operators...)
end

function _embedded_local_hamiltonian(subsystem::Transmon, dims::Vector{Int}, index::Int, annihilator)
    # Phase 1 uses a Duffing approximation so `ncut` acts as the local Hilbert-space truncation size.
    local_hamiltonian = _duffing_local_hamiltonian(subsystem.EJ, subsystem.EC, dims[index], annihilator; label = string(name(subsystem)))
    return _embed_operator(local_hamiltonian, dims, index)
end

function _embedded_local_hamiltonian(subsystem::Resonator, dims::Vector{Int}, index::Int, annihilator)
    number_operator = dag(annihilator) * annihilator
    local_hamiltonian = subsystem.ω * number_operator
    return _embed_operator(local_hamiltonian, dims, index)
end

function _embedded_local_hamiltonian(subsystem::TunableTransmon, dims::Vector{Int}, index::Int, annihilator)
    effective_EJ = _effective_josephson_energy(subsystem)
    local_hamiltonian =
        _duffing_local_hamiltonian(effective_EJ, subsystem.EC, dims[index], annihilator; label = string(name(subsystem)))
    return _embed_operator(local_hamiltonian, dims, index)
end

function _embedded_local_hamiltonian(subsystem::TunableCoupler, dims::Vector{Int}, index::Int, annihilator)
    effective_EJ = _effective_josephson_energy(subsystem)
    local_hamiltonian =
        _duffing_local_hamiltonian(effective_EJ, subsystem.EC, dims[index], annihilator; label = string(name(subsystem)))
    return _embed_operator(local_hamiltonian, dims, index)
end

function _embedded_local_hamiltonian(subsystem::AbstractSubsystem, dims::Vector{Int}, index::Int, annihilator)
    throw(ArgumentError("Unsupported subsystem type $(typeof(subsystem))."))
end

function _coupling_hamiltonian(coupling::CapacitiveCoupling, annihilation_ops::Dict{Symbol, Any})
    source, target = coupling_endpoints(coupling)
    a = get(annihilation_ops, source, nothing)
    b = get(annihilation_ops, target, nothing)
    a === nothing && throw(ArgumentError("Missing annihilation operator for subsystem $source."))
    b === nothing && throw(ArgumentError("Missing annihilation operator for subsystem $target."))

    return coupling.g * (dag(a) * b + a * dag(b))
end

function _coupling_hamiltonian(coupling::AbstractCoupling, annihilation_ops::Dict{Symbol, Any})
    throw(ArgumentError("Unsupported coupling type $(typeof(coupling))."))
end

function _effective_josephson_energy(subsystem::Union{TunableTransmon, TunableCoupler})
    return subsystem.EJmax * sqrt(cospi(subsystem.flux)^2 + subsystem.asymmetry^2 * sinpi(subsystem.flux)^2)
end

function _duffing_local_hamiltonian(EJ::Float64, EC::Float64, dimension::Int, annihilator; label::AbstractString)
    number_operator = dag(annihilator) * annihilator
    identity_operator = eye(dimension)
    ω = sqrt(8 * EJ * EC) - EC
    ω > 0 || throw(ArgumentError("Duffing approximation for $label requires a positive local frequency; choose parameters with larger effective EJ."))
    α = -EC
    return ω * number_operator + (α / 2) * number_operator * (number_operator - identity_operator)
end

include("Operators.jl")
include("States.jl")
include("Drives.jl")

end
