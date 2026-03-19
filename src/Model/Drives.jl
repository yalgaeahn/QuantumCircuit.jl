struct SubsystemDrive{F}
    label::Symbol
    target::Symbol
    operator::Symbol
    coefficient::F

    function SubsystemDrive(label::Symbol, target::Symbol, operator::Symbol, coefficient::F) where {F}
        _validate_nonempty_symbol(label, "drive label")
        _validate_nonempty_symbol(target, "drive target")
        _validate_nonempty_symbol(operator, "drive operator")
        coefficient isa Function || throw(ArgumentError("drive coefficient must be callable."))
        return new{F}(label, target, operator, coefficient)
    end
end

struct FluxControl{F, D}
    label::Symbol
    target::Symbol
    modulation::F
    derivative::D

    function FluxControl(label::Symbol, target::Symbol, modulation::F; derivative = nothing) where {F}
        _validate_nonempty_symbol(label, "flux-control label")
        _validate_nonempty_symbol(target, "flux-control target")
        modulation isa Function || throw(ArgumentError("flux-control modulation must be callable."))
        isnothing(derivative) || derivative isa Function || throw(ArgumentError("flux-control derivative must be callable when provided."))
        return new{F, typeof(derivative)}(label, target, modulation, derivative)
    end
end

function _time_dependent_hamiltonian(model::StaticSystemModel, drives::Nothing)
    return _time_dependent_hamiltonian(model, drives, nothing)
end

function _time_dependent_hamiltonian(model::StaticSystemModel, drives::Union{Tuple, AbstractVector})
    return _time_dependent_hamiltonian(model, drives, nothing)
end

function _time_dependent_hamiltonian(model::StaticSystemModel, drives::Nothing, flux_controls::Nothing)
    return hamiltonian(model)
end

function _time_dependent_hamiltonian(model::StaticSystemModel, drives, flux_controls)
    normalized_drives = _normalize_drives(drives)
    normalized_flux_controls = _normalize_flux_controls(flux_controls)

    isempty(normalized_drives) && isempty(normalized_flux_controls) && return hamiltonian(model)

    _validate_flux_controls(model, normalized_flux_controls)

    terms = Any[hamiltonian(model)]
    _append_drive_terms!(terms, model, normalized_drives)
    _append_flux_control_terms!(terms, model, normalized_flux_controls, model.hamiltonian_spec)

    return Tuple(terms)
end

function _normalize_drives(drives::Nothing)
    return SubsystemDrive[]
end

function _normalize_drives(drives::Union{Tuple, AbstractVector})
    normalized = SubsystemDrive[]
    for drive in drives
        drive isa SubsystemDrive ||
            throw(ArgumentError("drives must contain only SubsystemDrive values, got $(typeof(drive))."))
        push!(normalized, drive)
    end
    return normalized
end

function _normalize_flux_controls(flux_controls::Nothing)
    return FluxControl[]
end

function _normalize_flux_controls(flux_controls::Union{Tuple, AbstractVector})
    normalized = FluxControl[]
    for control in flux_controls
        control isa FluxControl ||
            throw(ArgumentError("flux_controls must contain only FluxControl values, got $(typeof(control))."))
        push!(normalized, control)
    end
    return normalized
end

function _append_drive_terms!(terms::Vector{Any}, model::StaticSystemModel, drives::Vector{SubsystemDrive})
    for drive in drives
        operator = _embedded_operator(model, drive.target, drive.operator)
        push!(terms, (operator, (params, t) -> drive.coefficient(params, t)))
    end

    return nothing
end

function _validate_flux_controls(model::StaticSystemModel, flux_controls::Vector{FluxControl})
    isempty(flux_controls) && return nothing

    labels = [control.label for control in flux_controls]
    length(unique(labels)) == length(labels) ||
        throw(ArgumentError("Flux-control labels must be unique within one evolve call."))

    targets = [control.target for control in flux_controls]
    length(unique(targets)) == length(targets) ||
        throw(ArgumentError("At most one FluxControl may target a given subsystem within one evolve call."))

    for control in flux_controls
        subsystem = _subsystem(model, control.target)
        subsystem isa _TunableSubsystem ||
            throw(ArgumentError("FluxControl target $(control.target) must be a TunableTransmon or TunableCoupler."))
    end

    _validate_flux_controls(model, flux_controls, model.hamiltonian_spec)
    return nothing
end

function _validate_flux_controls(
    model::StaticSystemModel,
    flux_controls::Vector{FluxControl},
    ::CircuitHamiltonianSpec,
)
    return nothing
end

function _validate_flux_controls(
    ::StaticSystemModel,
    flux_controls::Vector{FluxControl},
    ::EffectiveHamiltonianSpec{DuffingEffectiveMethod},
)
    throw(
        ArgumentError(
            "FluxControl requires EffectiveHamiltonianSpec(NonadiabaticDuffingEffectiveMethod()) for effective-mode dynamics.",
        ),
    )
end

function _validate_flux_controls(
    model::StaticSystemModel,
    flux_controls::Vector{FluxControl},
    ::EffectiveHamiltonianSpec{NonadiabaticDuffingEffectiveMethod},
)
    isempty(couplings(model.system)) ||
        throw(ArgumentError("NonadiabaticDuffingEffectiveMethod with flux_controls currently supports only uncoupled systems."))
    length(subsystems(model.system)) == 1 ||
        throw(
            ArgumentError(
                "NonadiabaticDuffingEffectiveMethod with flux_controls currently supports only single-subsystem systems.",
            ),
        )

    for control in flux_controls
        isnothing(control.derivative) &&
            throw(ArgumentError("FluxControl $(control.label) requires a derivative when NonadiabaticDuffingEffectiveMethod() is used."))
    end

    return nothing
end

function _validate_flux_controls(
    ::StaticSystemModel,
    flux_controls::Vector{FluxControl},
    hamiltonian_spec::EffectiveHamiltonianSpec,
)
    throw(
        ArgumentError(
            "FluxControl is not implemented for effective Hamiltonian method $(nameof(typeof(hamiltonian_spec.method))).",
        ),
    )
end

function _append_flux_control_terms!(
    terms::Vector{Any},
    model::StaticSystemModel,
    flux_controls::Vector{FluxControl},
    ::CircuitHamiltonianSpec,
)
    for control in flux_controls
        subsystem = _subsystem(model, control.target)
        _append_circuit_flux_terms!(terms, model, control, subsystem)
    end

    return nothing
end

function _append_flux_control_terms!(
    terms::Vector{Any},
    model::StaticSystemModel,
    flux_controls::Vector{FluxControl},
    ::EffectiveHamiltonianSpec{NonadiabaticDuffingEffectiveMethod},
)
    for control in flux_controls
        subsystem = _subsystem(model, control.target)
        _append_nonadiabatic_flux_terms!(terms, model, control, subsystem)
    end

    return nothing
end

function _append_flux_control_terms!(
    terms::Vector{Any},
    model::StaticSystemModel,
    flux_controls::Vector{FluxControl},
    ::AbstractHamiltonianSpec,
)
    isempty(flux_controls) || throw(ArgumentError("FluxControl is not implemented for $(nameof(typeof(model.hamiltonian_spec)))."))
    return nothing
end

function _append_circuit_flux_terms!(
    terms::Vector{Any},
    model::StaticSystemModel,
    control::FluxControl,
    subsystem::_TunableSubsystem,
)
    base_flux = subsystem.flux
    base_cos, base_sin = _effective_circuit_coefficients(subsystem, base_flux)
    cosphi_operator = _embedded_operator(model, control.target, :cosphi)
    sinphi_operator = _embedded_operator(model, control.target, :sinphi)

    push!(
        terms,
        (
            cosphi_operator,
            (params, t) -> begin
                flux = _modulated_flux(subsystem, control, params, t)
                shifted_cos, _ = _effective_circuit_coefficients(subsystem, flux)
                return -(shifted_cos - base_cos)
            end,
        ),
    )
    push!(
        terms,
        (
            sinphi_operator,
            (params, t) -> begin
                flux = _modulated_flux(subsystem, control, params, t)
                _, shifted_sin = _effective_circuit_coefficients(subsystem, flux)
                return -(shifted_sin - base_sin)
            end,
        ),
    )

    return nothing
end

function _append_nonadiabatic_flux_terms!(
    terms::Vector{Any},
    model::StaticSystemModel,
    control::FluxControl,
    subsystem::_TunableSubsystem,
)
    label = string(name(subsystem))
    base_frequency, _ = _duffing_parameters(_effective_josephson_energy(subsystem), subsystem.EC; label = label)
    number_operator = _embedded_operator(model, control.target, :n)
    quadrature_y = _embedded_operator(model, control.target, :y)
    pair_operator = 1im * (_embedded_operator(model, control.target, :adagadag) - _embedded_operator(model, control.target, :aa))

    push!(terms, (number_operator, (params, t) -> _nonadiabatic_frequency_shift(subsystem, control, base_frequency, params, t)))
    push!(terms, (quadrature_y, (params, t) -> _nonadiabatic_y_coefficient(subsystem, control, params, t)))
    push!(terms, (pair_operator, (params, t) -> _nonadiabatic_pair_coefficient(subsystem, control, params, t)))

    return nothing
end

function _modulated_flux(subsystem::_TunableSubsystem, control::FluxControl, params, t)
    return subsystem.flux + control.modulation(params, t)
end

function _modulated_flux_derivative(control::FluxControl, params, t)
    isnothing(control.derivative) && return 0.0
    return control.derivative(params, t)
end

function _nonadiabatic_frequency_shift(
    subsystem::_TunableSubsystem,
    control::FluxControl,
    base_frequency::Float64,
    params,
    t,
)
    flux = _modulated_flux(subsystem, control, params, t)
    shifted_frequency, _ = _duffing_parameters(
        _effective_josephson_energy(subsystem.EJmax, flux, subsystem.asymmetry),
        subsystem.EC;
        label = string(name(subsystem)),
    )
    return shifted_frequency - base_frequency
end

function _nonadiabatic_y_coefficient(subsystem::_TunableSubsystem, control::FluxControl, params, t)
    flux = _modulated_flux(subsystem, control, params, t)
    flux_rate = _modulated_flux_derivative(control, params, t)
    xi = _effective_xi(_effective_josephson_energy(subsystem.EJmax, flux, subsystem.asymmetry), subsystem.EC)
    return -sqrt(xi / 2) * _effective_phi_eff_derivative(flux, subsystem.asymmetry, flux_rate)
end

function _nonadiabatic_pair_coefficient(subsystem::_TunableSubsystem, control::FluxControl, params, t)
    flux = _modulated_flux(subsystem, control, params, t)
    flux_rate = _modulated_flux_derivative(control, params, t)
    return _effective_xi_log_derivative(flux, subsystem.asymmetry, flux_rate) / 4
end

function _effective_xi(EJ::Real, EC::Real)
    return sqrt(EJ / (8 * EC))
end

function _effective_phi_eff_derivative(flux::Real, asymmetry::Real, flux_rate::Real)
    denominator = cospi(flux)^2 + asymmetry^2 * sinpi(flux)^2
    return pi * asymmetry * flux_rate / denominator
end

function _effective_xi_log_derivative(flux::Real, asymmetry::Real, flux_rate::Real)
    denominator = cospi(flux)^2 + asymmetry^2 * sinpi(flux)^2
    return pi * flux_rate * (asymmetry^2 - 1) * sinpi(2 * flux) / (4 * denominator)
end
