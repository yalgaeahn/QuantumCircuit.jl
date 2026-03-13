abstract type AbstractSubsystem end

name(subsystem::AbstractSubsystem) = subsystem.name

struct Transmon <: AbstractSubsystem
    name::Symbol
    EJ::Float64
    EC::Float64
    ng::Float64
    ncut::Int

    function Transmon(name::Symbol; EJ::Real, EC::Real, ng::Real = 0.0, ncut::Integer = 15)
        new(
            validate_identifier(name, "subsystem"),
            validate_positive(EJ, "EJ"),
            validate_positive(EC, "EC"),
            validate_finite(ng, "ng"),
            validate_positive_integer(ncut, "ncut"),
        )
    end
end

struct Resonator <: AbstractSubsystem
    name::Symbol
    ω::Float64
    dim::Int

    function Resonator(name::Symbol; ω::Real, dim::Integer)
        new(
            validate_identifier(name, "subsystem"),
            validate_positive(ω, "ω"),
            validate_positive_integer(dim, "dim"),
        )
    end
end

struct TunableTransmon <: AbstractSubsystem
    name::Symbol
    EJmax::Float64
    EC::Float64
    flux::Float64
    asymmetry::Float64
    ncut::Int

    function TunableTransmon(
        name::Symbol;
        EJmax::Real,
        EC::Real,
        flux::Real = 0.0,
        asymmetry::Real = 0.0,
        ncut::Integer = 15,
    )
        new(
            validate_identifier(name, "subsystem"),
            validate_positive(EJmax, "EJmax"),
            validate_positive(EC, "EC"),
            validate_finite(flux, "flux"),
            validate_abs_leq_one(asymmetry, "asymmetry"),
            validate_positive_integer(ncut, "ncut"),
        )
    end
end

struct TunableCoupler <: AbstractSubsystem
    name::Symbol
    EJmax::Float64
    EC::Float64
    flux::Float64
    asymmetry::Float64
    ncut::Int

    function TunableCoupler(
        name::Symbol;
        EJmax::Real,
        EC::Real,
        flux::Real = 0.0,
        asymmetry::Real = 0.0,
        ncut::Integer = 11,
    )
        new(
            validate_identifier(name, "subsystem"),
            validate_positive(EJmax, "EJmax"),
            validate_positive(EC, "EC"),
            validate_finite(flux, "flux"),
            validate_abs_leq_one(asymmetry, "asymmetry"),
            validate_positive_integer(ncut, "ncut"),
        )
    end
end
