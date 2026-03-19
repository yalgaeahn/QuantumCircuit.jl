module Model

using QuantumToolbox: Operator, QuantumObject, basis, dag, destroy, eye, kron
using SparseArrays: sparse, spdiagm, spzeros
using ..Architecture:
    AbstractCoupling,
    AbstractSubsystem,
    CapacitiveCoupling,
    CircuitCapacitiveCoupling,
    CompositeSystem,
    Resonator,
    TunableCoupler,
    TunableTransmon,
    Transmon,
    coupling_endpoints,
    couplings,
    name,
    subsystems

export AbstractEffectiveMethod,
    AbstractHamiltonianSpec,
    CircuitHamiltonianSpec,
    DuffingEffectiveMethod,
    EffectiveHamiltonianSpec,
    FluxControl,
    NonadiabaticDuffingEffectiveMethod,
    ObservableSpec,
    StaticSystemModel,
    SubsystemDrive,
    annihilation_operator,
    basis_state,
    build_model,
    charge_operator,
    cosphi_operator,
    creation_operator,
    hamiltonian,
    number_operator,
    quadrature_operator,
    sinphi_operator

abstract type AbstractHamiltonianSpec end
abstract type AbstractEffectiveMethod end

struct DuffingEffectiveMethod <: AbstractEffectiveMethod end
struct NonadiabaticDuffingEffectiveMethod <: AbstractEffectiveMethod end

struct EffectiveHamiltonianSpec{M<:AbstractEffectiveMethod} <: AbstractHamiltonianSpec
    method::M
end

EffectiveHamiltonianSpec() = EffectiveHamiltonianSpec(DuffingEffectiveMethod())

struct CircuitHamiltonianSpec <: AbstractHamiltonianSpec
    charge_cutoff::Int
    charge_cutoffs::Dict{Symbol, Int}

    function CircuitHamiltonianSpec(; charge_cutoff::Integer, charge_cutoffs::AbstractDict{Symbol, <:Integer} = Dict{Symbol, Int}())
        validated_cutoffs = Dict{Symbol, Int}()
        for (subsystem_name, cutoff) in pairs(charge_cutoffs)
            _validate_nonempty_symbol(subsystem_name, "charge cutoff override subsystem")
            validated_cutoffs[subsystem_name] =
                _validate_nonnegative_integer(cutoff, "charge cutoff override for subsystem $subsystem_name")
        end

        return new(
            _validate_nonnegative_integer(charge_cutoff, "charge_cutoff"),
            validated_cutoffs,
        )
    end
end

const _CircuitChargeSubsystem = Union{Transmon, TunableTransmon, TunableCoupler}
const _DuffingLikeEffectiveMethod = Union{DuffingEffectiveMethod, NonadiabaticDuffingEffectiveMethod}
const _TunableSubsystem = Union{TunableTransmon, TunableCoupler}
const _OperatorCache = Dict{Tuple{Symbol, Symbol}, Any}
const _LocalOperatorCache = Dict{Symbol, Any}

struct StaticSystemModel{HS<:AbstractHamiltonianSpec, QT}
    system::CompositeSystem
    hamiltonian_spec::HS
    subsystem_order::Vector{Symbol}
    dimensions::Dict{Symbol, Int}
    hamiltonian::QT
    operator_cache::_OperatorCache
end

hamiltonian(model::StaticSystemModel) = model.hamiltonian

function _validate_nonempty_symbol(value::Symbol, label::AbstractString)
    value == Symbol("") && throw(ArgumentError("$label must not be empty."))
    return value
end

function _validate_nonnegative_integer(value, label::AbstractString)
    value isa Integer || throw(ArgumentError("$label must be an integer, got $(typeof(value))."))
    value >= 0 || throw(ArgumentError("$label must be nonnegative."))
    return Int(value)
end

function build_model(system::CompositeSystem; hamiltonian_spec::AbstractHamiltonianSpec = EffectiveHamiltonianSpec())
    subsystem_list = subsystems(system)
    _validate_hamiltonian_spec(system, subsystem_list, hamiltonian_spec)

    dims = [_local_dimension(subsystem, hamiltonian_spec) for subsystem in subsystem_list]
    order = Symbol[name(subsystem) for subsystem in subsystem_list]
    operator_cache = _OperatorCache()
    H = _zero_operator(dims)

    for (index, subsystem) in enumerate(subsystem_list)
        local_operators = _local_operator_bundle(subsystem, hamiltonian_spec, dims[index])
        local_hamiltonian = _local_hamiltonian(subsystem, hamiltonian_spec, dims[index], local_operators)
        H += _embed_operator(local_hamiltonian, dims, index)
        _cache_embedded_operators!(operator_cache, local_operators, subsystem, dims, index)
    end

    for coupling in couplings(system)
        H += _coupling_hamiltonian(coupling, system, operator_cache, hamiltonian_spec)
    end

    return StaticSystemModel(
        system,
        hamiltonian_spec,
        order,
        Dict(order[i] => dims[i] for i in eachindex(order)),
        H,
        operator_cache,
    )
end

function hamiltonian(system::CompositeSystem; hamiltonian_spec::AbstractHamiltonianSpec = EffectiveHamiltonianSpec())
    return hamiltonian(build_model(system; hamiltonian_spec = hamiltonian_spec))
end

function _validate_hamiltonian_spec(
    system::CompositeSystem,
    ::Vector{<:AbstractSubsystem},
    ::EffectiveHamiltonianSpec,
)
    for coupling in couplings(system)
        coupling isa CapacitiveCoupling && continue
        coupling isa CircuitCapacitiveCoupling &&
            throw(ArgumentError("CircuitCapacitiveCoupling(...; G = ...) requires CircuitHamiltonianSpec(...)."))
        throw(ArgumentError("Unsupported coupling type $(typeof(coupling))."))
    end

    return nothing
end

function _validate_hamiltonian_spec(
    system::CompositeSystem,
    subsystem_list::Vector{<:AbstractSubsystem},
    hamiltonian_spec::CircuitHamiltonianSpec,
)
    subsystem_by_name = Dict(name(subsystem) => subsystem for subsystem in subsystem_list)
    for subsystem_name in keys(hamiltonian_spec.charge_cutoffs)
        subsystem = get(subsystem_by_name, subsystem_name, nothing)
        subsystem === nothing &&
            throw(ArgumentError("Charge cutoff override references missing subsystem $subsystem_name."))
        subsystem isa _CircuitChargeSubsystem ||
            throw(
                ArgumentError(
                    "Charge cutoff override for subsystem $subsystem_name is only valid for transmon-like subsystems in circuit mode.",
                ),
            )
    end

    for coupling in couplings(system)
        _validate_circuit_coupling(coupling, subsystem_by_name)
    end

    return nothing
end

function _validate_circuit_coupling(
    coupling::CapacitiveCoupling,
    subsystem_by_name::Dict{Symbol, <:AbstractSubsystem},
)
    throw(
        ArgumentError(
            "CapacitiveCoupling(...; g = ...) is only supported with EffectiveHamiltonianSpec(). Use CircuitCapacitiveCoupling(...; G = ...) with CircuitHamiltonianSpec().",
        ),
    )
end

function _validate_circuit_coupling(
    coupling::CircuitCapacitiveCoupling,
    subsystem_by_name::Dict{Symbol, <:AbstractSubsystem},
)
    source, target = coupling_endpoints(coupling)
    source_subsystem = subsystem_by_name[source]
    target_subsystem = subsystem_by_name[target]

    source_subsystem isa Resonator && target_subsystem isa Resonator &&
        throw(ArgumentError("CircuitCapacitiveCoupling between resonators is not implemented yet."))

    _supports_exact_circuit_coupling(source_subsystem, target_subsystem) ||
        throw(
            ArgumentError(
                "CircuitCapacitiveCoupling under CircuitHamiltonianSpec supports only transmon-like↔transmon-like and resonator↔transmon-like pairs.",
            ),
        )

    return nothing
end

function _validate_circuit_coupling(
    coupling::AbstractCoupling,
    subsystem_by_name::Dict{Symbol, <:AbstractSubsystem},
)
    throw(ArgumentError("Unsupported coupling type $(typeof(coupling)) under CircuitHamiltonianSpec."))
end

function _supports_exact_circuit_coupling(source::AbstractSubsystem, target::AbstractSubsystem)
    return (source isa _CircuitChargeSubsystem && target isa _CircuitChargeSubsystem) ||
        (source isa Resonator && target isa _CircuitChargeSubsystem) ||
        (source isa _CircuitChargeSubsystem && target isa Resonator)
end

_local_dimension(subsystem::Transmon, ::EffectiveHamiltonianSpec) = subsystem.ncut
_local_dimension(subsystem::Resonator, ::AbstractHamiltonianSpec) = subsystem.dim
_local_dimension(subsystem::TunableTransmon, ::EffectiveHamiltonianSpec) = subsystem.ncut
_local_dimension(subsystem::TunableCoupler, ::EffectiveHamiltonianSpec) = subsystem.ncut

_local_dimension(subsystem::Transmon, hamiltonian_spec::CircuitHamiltonianSpec) = 2 * _charge_cutoff(hamiltonian_spec, subsystem) + 1
_local_dimension(subsystem::TunableTransmon, hamiltonian_spec::CircuitHamiltonianSpec) = 2 * _charge_cutoff(hamiltonian_spec, subsystem) + 1
_local_dimension(subsystem::TunableCoupler, hamiltonian_spec::CircuitHamiltonianSpec) = 2 * _charge_cutoff(hamiltonian_spec, subsystem) + 1

_local_dimension(subsystem::AbstractSubsystem, ::AbstractHamiltonianSpec) =
    throw(ArgumentError("Unsupported subsystem type $(typeof(subsystem))."))

function _charge_cutoff(hamiltonian_spec::CircuitHamiltonianSpec, subsystem::_CircuitChargeSubsystem)
    return get(hamiltonian_spec.charge_cutoffs, name(subsystem), hamiltonian_spec.charge_cutoff)
end

function _zero_operator(dims::Vector{Int})
    identity_factors = [eye(dim) for dim in dims]
    return 0.0 * kron(identity_factors...)
end

function _embed_operator(local_operator, dims::Vector{Int}, index::Int)
    operators = [position == index ? local_operator : eye(dims[position]) for position in eachindex(dims)]
    return kron(operators...)
end

function _cache_embedded_operators!(
    operator_cache::_OperatorCache,
    local_operators::_LocalOperatorCache,
    subsystem::AbstractSubsystem,
    dims::Vector{Int},
    index::Int,
)
    subsystem_name = name(subsystem)
    for (operator_symbol, local_operator) in pairs(local_operators)
        operator_cache[(subsystem_name, operator_symbol)] = _embed_operator(local_operator, dims, index)
    end
    return nothing
end

function _local_operator_bundle(::Union{Transmon, TunableTransmon, TunableCoupler, Resonator}, ::EffectiveHamiltonianSpec, dimension::Int)
    return _oscillator_operator_bundle(dimension)
end

function _local_operator_bundle(::Resonator, ::CircuitHamiltonianSpec, dimension::Int)
    return _oscillator_operator_bundle(dimension)
end

function _local_operator_bundle(::Union{Transmon, TunableTransmon, TunableCoupler}, ::CircuitHamiltonianSpec, dimension::Int)
    return _charge_basis_operator_bundle(dimension)
end

function _local_operator_bundle(subsystem::AbstractSubsystem, ::AbstractHamiltonianSpec, ::Int)
    throw(ArgumentError("Unsupported subsystem type $(typeof(subsystem))."))
end

function _oscillator_operator_bundle(dimension::Int)
    annihilator = destroy(dimension)
    creator = dag(annihilator)
    number_operator = creator * annihilator
    return _LocalOperatorCache(
        :a => annihilator,
        :adag => creator,
        :aa => annihilator * annihilator,
        :adagadag => creator * creator,
        :n => number_operator,
        :x => annihilator + creator,
        :y => -1im * (annihilator - creator),
    )
end

function _charge_basis_operator_bundle(dimension::Int)
    cutoff = (dimension - 1) ÷ 2
    charge_values = ComplexF64.(collect(-cutoff:cutoff))
    cosphi_matrix = spzeros(ComplexF64, dimension, dimension)
    sinphi_matrix = spzeros(ComplexF64, dimension, dimension)

    for index in 1:(dimension - 1)
        cosphi_matrix[index, index + 1] = 0.5
        cosphi_matrix[index + 1, index] = 0.5
        sinphi_matrix[index, index + 1] = -0.5im
        sinphi_matrix[index + 1, index] = 0.5im
    end

    return _LocalOperatorCache(
        :charge => _custom_operator(spdiagm(0 => charge_values), dimension),
        :cosphi => _custom_operator(cosphi_matrix, dimension),
        :sinphi => _custom_operator(sinphi_matrix, dimension),
    )
end

function _custom_operator(matrix, dimension::Int)
    return QuantumObject(sparse(matrix); type = Operator(), dims = eye(dimension).dimensions)
end

function _local_hamiltonian(
    subsystem::Transmon,
    hamiltonian_spec::EffectiveHamiltonianSpec{<:_DuffingLikeEffectiveMethod},
    dimension::Int,
    local_operators::_LocalOperatorCache,
)
    return _duffing_local_hamiltonian(subsystem.EJ, subsystem.EC, dimension, local_operators[:n]; label = string(name(subsystem)))
end

function _local_hamiltonian(
    subsystem::TunableTransmon,
    hamiltonian_spec::EffectiveHamiltonianSpec{<:_DuffingLikeEffectiveMethod},
    dimension::Int,
    local_operators::_LocalOperatorCache,
)
    effective_EJ = _effective_josephson_energy(subsystem)
    return _duffing_local_hamiltonian(
        effective_EJ,
        subsystem.EC,
        dimension,
        local_operators[:n];
        label = string(name(subsystem)),
    )
end

function _local_hamiltonian(
    subsystem::TunableCoupler,
    hamiltonian_spec::EffectiveHamiltonianSpec{<:_DuffingLikeEffectiveMethod},
    dimension::Int,
    local_operators::_LocalOperatorCache,
)
    effective_EJ = _effective_josephson_energy(subsystem)
    return _duffing_local_hamiltonian(
        effective_EJ,
        subsystem.EC,
        dimension,
        local_operators[:n];
        label = string(name(subsystem)),
    )
end

function _local_hamiltonian(
    subsystem::Resonator,
    ::EffectiveHamiltonianSpec{<:_DuffingLikeEffectiveMethod},
    ::Int,
    local_operators::_LocalOperatorCache,
)
    return subsystem.ω * local_operators[:n]
end

function _local_hamiltonian(
    subsystem::Resonator,
    ::CircuitHamiltonianSpec,
    ::Int,
    local_operators::_LocalOperatorCache,
)
    return subsystem.ω * local_operators[:n]
end

function _local_hamiltonian(
    subsystem::Transmon,
    ::CircuitHamiltonianSpec,
    ::Int,
    local_operators::_LocalOperatorCache,
)
    return _circuit_transmon_local_hamiltonian(subsystem.EJ, subsystem.EC, subsystem.ng, local_operators)
end

function _local_hamiltonian(
    subsystem::TunableTransmon,
    ::CircuitHamiltonianSpec,
    ::Int,
    local_operators::_LocalOperatorCache,
)
    return _circuit_tunable_local_hamiltonian(
        subsystem.EJmax,
        subsystem.asymmetry,
        subsystem.flux,
        subsystem.EC,
        subsystem.ng,
        local_operators,
    )
end

function _local_hamiltonian(
    subsystem::TunableCoupler,
    ::CircuitHamiltonianSpec,
    ::Int,
    local_operators::_LocalOperatorCache,
)
    return _circuit_tunable_local_hamiltonian(
        subsystem.EJmax,
        subsystem.asymmetry,
        subsystem.flux,
        subsystem.EC,
        subsystem.ng,
        local_operators,
    )
end

function _local_hamiltonian(
    subsystem::Union{Transmon, TunableTransmon, TunableCoupler},
    hamiltonian_spec::EffectiveHamiltonianSpec,
    ::Int,
    ::_LocalOperatorCache,
)
    throw(ArgumentError("Effective Hamiltonian method $(nameof(typeof(hamiltonian_spec.method))) is not implemented."))
end

function _local_hamiltonian(subsystem::AbstractSubsystem, ::AbstractHamiltonianSpec, ::Int, ::_LocalOperatorCache)
    throw(ArgumentError("Unsupported subsystem type $(typeof(subsystem))."))
end

function _coupling_hamiltonian(
    coupling::CapacitiveCoupling,
    system::CompositeSystem,
    operator_cache::_OperatorCache,
    ::EffectiveHamiltonianSpec,
)
    source, target = coupling_endpoints(coupling)
    a = get(operator_cache, (source, :a), nothing)
    b = get(operator_cache, (target, :a), nothing)
    a === nothing && throw(ArgumentError("Missing annihilation operator for subsystem $source."))
    b === nothing && throw(ArgumentError("Missing annihilation operator for subsystem $target."))
    return coupling.g * (dag(a) * b + a * dag(b))
end

function _coupling_hamiltonian(
    coupling::CapacitiveCoupling,
    system::CompositeSystem,
    operator_cache::_OperatorCache,
    ::CircuitHamiltonianSpec,
)
    throw(
        ArgumentError(
            "CapacitiveCoupling(...; g = ...) is only supported with EffectiveHamiltonianSpec(). Use CircuitCapacitiveCoupling(...; G = ...) with CircuitHamiltonianSpec().",
        ),
    )
end

function _coupling_hamiltonian(
    coupling::CircuitCapacitiveCoupling,
    system::CompositeSystem,
    operator_cache::_OperatorCache,
    ::CircuitHamiltonianSpec,
)
    source, target = coupling_endpoints(coupling)
    source_subsystem = _subsystem(system, source)
    target_subsystem = _subsystem(system, target)

    if source_subsystem isa _CircuitChargeSubsystem && target_subsystem isa _CircuitChargeSubsystem
        source_charge = _cached_operator(operator_cache, source, :charge)
        target_charge = _cached_operator(operator_cache, target, :charge)
        return coupling.G * source_charge * target_charge
    elseif source_subsystem isa Resonator && target_subsystem isa _CircuitChargeSubsystem
        source_x = _cached_operator(operator_cache, source, :x)
        target_charge = _cached_operator(operator_cache, target, :charge)
        return coupling.G * source_x * target_charge
    elseif source_subsystem isa _CircuitChargeSubsystem && target_subsystem isa Resonator
        source_charge = _cached_operator(operator_cache, source, :charge)
        target_x = _cached_operator(operator_cache, target, :x)
        return coupling.G * source_charge * target_x
    elseif source_subsystem isa Resonator && target_subsystem isa Resonator
        throw(ArgumentError("CircuitCapacitiveCoupling between resonators is not implemented yet."))
    end

    throw(
        ArgumentError(
            "CircuitCapacitiveCoupling under CircuitHamiltonianSpec supports only transmon-like↔transmon-like and resonator↔transmon-like pairs.",
        ),
    )
end

function _coupling_hamiltonian(
    coupling::CircuitCapacitiveCoupling,
    system::CompositeSystem,
    operator_cache::_OperatorCache,
    ::EffectiveHamiltonianSpec,
)
    throw(ArgumentError("CircuitCapacitiveCoupling(...; G = ...) requires CircuitHamiltonianSpec(...)."))
end

function _coupling_hamiltonian(
    coupling::AbstractCoupling,
    system::CompositeSystem,
    operator_cache::_OperatorCache,
    ::AbstractHamiltonianSpec,
)
    throw(ArgumentError("Unsupported coupling type $(typeof(coupling))."))
end

function _cached_operator(operator_cache::_OperatorCache, target::Symbol, operator::Symbol)
    embedded_operator = get(operator_cache, (target, operator), nothing)
    embedded_operator === nothing &&
        throw(ArgumentError("Missing embedded operator $operator for subsystem $target."))
    return embedded_operator
end

function _effective_josephson_energy(subsystem::Union{TunableTransmon, TunableCoupler})
    return _effective_josephson_energy(subsystem.EJmax, subsystem.flux, subsystem.asymmetry)
end

function _effective_josephson_energy(EJmax::Float64, flux::Float64, asymmetry::Float64)
    return EJmax * sqrt(cospi(flux)^2 + asymmetry^2 * sinpi(flux)^2)
end

function _effective_josephson_energy(EJmax::Real, flux::Real, asymmetry::Real)
    return Float64(EJmax) * sqrt(cospi(Float64(flux))^2 + Float64(asymmetry)^2 * sinpi(Float64(flux))^2)
end

function _effective_circuit_coefficients(EJmax::Float64, flux::Float64, asymmetry::Float64)
    return EJmax * cospi(flux), EJmax * asymmetry * sinpi(flux)
end

function _effective_circuit_coefficients(subsystem::_TunableSubsystem, flux::Real)
    return _effective_circuit_coefficients(subsystem.EJmax, Float64(flux), subsystem.asymmetry)
end

function _duffing_local_hamiltonian(
    EJ::Float64,
    EC::Float64,
    dimension::Int,
    number_operator;
    label::AbstractString,
)
    identity_operator = eye(dimension)
    ω, α = _duffing_parameters(EJ, EC; label = label)
    return ω * number_operator + (α / 2) * number_operator * (number_operator - identity_operator)
end

function _duffing_parameters(EJ::Float64, EC::Float64; label::AbstractString)
    ω = sqrt(8 * EJ * EC) - EC
    ω > 0 || throw(ArgumentError("Duffing approximation for $label requires a positive local frequency; choose parameters with larger effective EJ."))
    α = -EC
    return ω, α
end

function _circuit_transmon_local_hamiltonian(EJ::Float64, EC::Float64, ng::Float64, local_operators::_LocalOperatorCache)
    charge_operator = local_operators[:charge]
    identity_operator = eye(size(charge_operator.data, 1))
    offset_charge = charge_operator - ng * identity_operator
    return 4 * EC * offset_charge * offset_charge - EJ * local_operators[:cosphi]
end

function _circuit_tunable_local_hamiltonian(
    EJmax::Float64,
    asymmetry::Float64,
    flux::Float64,
    EC::Float64,
    ng::Float64,
    local_operators::_LocalOperatorCache,
)
    charge_operator = local_operators[:charge]
    identity_operator = eye(size(charge_operator.data, 1))
    offset_charge = charge_operator - ng * identity_operator
    josephson_cos, josephson_sin = _effective_circuit_coefficients(EJmax, flux, asymmetry)
    return 4 * EC * offset_charge * offset_charge -
        josephson_cos * local_operators[:cosphi] -
        josephson_sin * local_operators[:sinphi]
end

function _subsystem(model::StaticSystemModel, target::Symbol)
    for subsystem in subsystems(model.system)
        name(subsystem) == target && return subsystem
    end

    throw(ArgumentError("No subsystem named $target exists in the model."))
end

function _subsystem(system::CompositeSystem, target::Symbol)
    for subsystem in subsystems(system)
        name(subsystem) == target && return subsystem
    end

    throw(ArgumentError("No subsystem named $target exists in the system."))
end

include("Operators.jl")
include("States.jl")
include("Drives.jl")

end
