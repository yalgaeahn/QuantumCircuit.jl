using Test
using QuantumCircuit
using QuantumToolbox: expect, eye
using TOML

effective_EJ(EJmax, flux, asymmetry) = EJmax * sqrt(cospi(flux)^2 + asymmetry^2 * sinpi(flux)^2)
transition_01(system::CompositeSystem; kwargs...) = transition_frequencies(spectrum(system; levels = 4, kwargs...))[1]
is_strictly_decreasing(values) = all(diff(values) .< 0)
wrapped_phase_distance(phase, target) = abs(mod(phase - target + π, 2π) - π)

function exact_tunable_circuit_hamiltonian(model, target, subsystem)
    charge = charge_operator(model, target)
    cosphi = cosphi_operator(model, target)
    sinphi = sinphi_operator(model, target)
    identity_operator = eye(size(charge.data, 1))
    offset_charge = charge - subsystem.ng * identity_operator
    return 4 * subsystem.EC * offset_charge * offset_charge -
        subsystem.EJmax * cospi(subsystem.flux) * cosphi -
        subsystem.EJmax * subsystem.asymmetry * sinpi(subsystem.flux) * sinphi
end

function exact_charge_charge_circuit_coupling(model, source, target, G)
    return G * charge_operator(model, source) * charge_operator(model, target)
end

function exact_resonator_charge_circuit_coupling(model, resonator, target, G)
    return G * quadrature_operator(model, resonator, :x) * charge_operator(model, target)
end

@testset "Architecture layer" begin
    q1 = Transmon(:q1; EJ = 20.0, EC = 0.25)
    r1 = Resonator(:r1; ω = 6.8, dim = 10)
    g1 = CapacitiveCoupling(:q1, :r1; g = 0.09)
    tq1 = TunableTransmon(:tq1; EJmax = 22.0, EC = 0.30, flux = 0.1, asymmetry = 0.05)
    c1 = TunableCoupler(:c1; EJmax = 15.0, EC = 0.35, flux = -0.2, asymmetry = 0.10)

    @test name(q1) == :q1
    @test q1.ng == 0.0
    @test q1.ncut == 15
    @test name(r1) == :r1
    @test r1.dim == 10
    @test name(tq1) == :tq1
    @test tq1.EJmax == 22.0
    @test tq1.flux == 0.1
    @test tq1.ng == 0.0
    @test name(c1) == :c1
    @test c1.ng == 0.0
    @test c1.ncut == 11
    @test coupling_endpoints(g1) == (:q1, :r1)

    sys = CompositeSystem(q1, r1, g1)
    @test subsystem_names(sys) == [:q1, :r1]
    @test length(subsystems(sys)) == 2
    @test length(couplings(sys)) == 1

    sys_from_vectors = CompositeSystem([q1, r1], [g1])
    @test subsystem_names(sys_from_vectors) == [:q1, :r1]
end

@testset "CompositeSystem validation" begin
    q1 = Transmon(:q1; EJ = 20.0, EC = 0.25)
    q1_dup = Transmon(:q1; EJ = 18.0, EC = 0.30)
    r1 = Resonator(:r1; ω = 6.8, dim = 10)
    tq1 = TunableTransmon(:tq1; EJmax = 22.0, EC = 0.30, flux = 0.1)
    c1 = TunableCoupler(:c1; EJmax = 15.0, EC = 0.35, flux = -0.2)
    g1 = CapacitiveCoupling(:q1, :c1; g = 0.09)
    g2 = CapacitiveCoupling(:c1, :r1; g = 0.06)

    @test_throws ArgumentError CompositeSystem(q1, q1_dup)
    @test_throws ArgumentError CompositeSystem(CapacitiveCoupling(:q1, :r1; g = 0.09))
    @test_throws ArgumentError CompositeSystem(q1, CapacitiveCoupling(:q1, :missing; g = 0.09))
    @test_throws ArgumentError CompositeSystem(q1, r1, 1)

    mixed_sys = CompositeSystem(q1, tq1, r1, c1, g1, g2)
    @test subsystem_names(mixed_sys) == [:q1, :tq1, :r1, :c1]
    @test length(couplings(mixed_sys)) == 2
end

@testset "Constructor validation" begin
    @test_throws ArgumentError Transmon(Symbol(""); EJ = 20.0, EC = 0.25)
    @test_throws ArgumentError Resonator(:r1; ω = -1.0, dim = 10)
    @test_throws ArgumentError CapacitiveCoupling(:q1, :q1; g = 0.09)
    @test_throws ArgumentError CircuitCapacitiveCoupling(:q1, :q1; G = 0.09)
    @test_throws ArgumentError TunableTransmon(:tq1; EJmax = 20.0, EC = 0.25, asymmetry = 1.2)
    @test_throws ArgumentError TunableCoupler(:c1; EJmax = -15.0, EC = 0.30, flux = 0.0)

    @test TunableTransmon(:tq1; EJmax = 20.0, EC = 0.25, ng = 0.2).ng == 0.2
    @test TunableCoupler(:c1; EJmax = 15.0, EC = 0.30).ng == 0.0
    @test CircuitCapacitiveCoupling(:q1, :q2; G = 0.05).G == 0.05
    @test EffectiveHamiltonianSpec().method isa DuffingEffectiveMethod
    @test CircuitHamiltonianSpec(charge_cutoff = 2).charge_cutoff == 2
    @test_throws UndefKeywordError CircuitHamiltonianSpec()
    @test_throws ArgumentError CircuitHamiltonianSpec(charge_cutoff = -1)
end

@testset "Model and spectrum" begin
    q1 = Transmon(:q1; EJ = 20.0, EC = 0.25, ncut = 6)
    transmon_sys = CompositeSystem(q1)

    model = build_model(transmon_sys)
    H = hamiltonian(model)
    @test size(H.data) == (6, 6)
    @test model.dimensions[:q1] == 6

    spec = spectrum(transmon_sys; levels = 4)
    @test length(spec.energies) == 4
    @test issorted(spec.energies)
    @test transition_frequencies(spec)[1] ≈ (sqrt(8 * q1.EJ * q1.EC) - q1.EC) atol = 1e-10
    @test anharmonicity(spec) ≈ -q1.EC atol = 1e-10
end

@testset "Coupled static Hamiltonian" begin
    q1 = Transmon(:q1; EJ = 20.0, EC = 0.25, ncut = 5)
    r1 = Resonator(:r1; ω = 6.8, dim = 4)
    g1 = CapacitiveCoupling(:q1, :r1; g = 0.09)

    sys = CompositeSystem(q1, r1, g1)
    H = hamiltonian(sys)
    @test size(H.data) == (20, 20)

    spec = spectrum(sys; levels = 5)
    @test length(spec.energies) == 5
end

@testset "Tunable model and mixed systems" begin
    tq_flux0 = TunableTransmon(:tq0; EJmax = 20.0, EC = 0.25, flux = 0.0, asymmetry = 0.0, ncut = 6)
    tq_flux1 = TunableTransmon(:tq1; EJmax = 20.0, EC = 0.25, flux = 0.25, asymmetry = 0.0, ncut = 6)

    spec0 = spectrum(CompositeSystem(tq_flux0); levels = 4)
    spec1 = spectrum(CompositeSystem(tq_flux1); levels = 4)

    ω0 = sqrt(8 * effective_EJ(tq_flux0.EJmax, tq_flux0.flux, tq_flux0.asymmetry) * tq_flux0.EC) - tq_flux0.EC
    ω1 = sqrt(8 * effective_EJ(tq_flux1.EJmax, tq_flux1.flux, tq_flux1.asymmetry) * tq_flux1.EC) - tq_flux1.EC

    @test transition_frequencies(spec0)[1] ≈ ω0 atol = 1e-10
    @test transition_frequencies(spec1)[1] ≈ ω1 atol = 1e-10
    @test transition_frequencies(spec1)[1] < transition_frequencies(spec0)[1]
    @test anharmonicity(spec1) ≈ -tq_flux1.EC atol = 1e-10

    q_fixed = Transmon(:qf; EJ = 19.0, EC = 0.24, ncut = 5)
    coupler = TunableCoupler(:c1; EJmax = 15.0, EC = 0.30, flux = 0.1, asymmetry = 0.1, ncut = 4)
    resonator = Resonator(:r1; ω = 6.7, dim = 3)
    g_qc = CapacitiveCoupling(:qf, :c1; g = 0.07)
    g_cr = CapacitiveCoupling(:c1, :r1; g = 0.05)

    mixed_sys = CompositeSystem(q_fixed, coupler, resonator, g_qc, g_cr)
    mixed_H = hamiltonian(mixed_sys)
    @test size(mixed_H.data) == (60, 60)

    mixed_spec = spectrum(mixed_sys; levels = 5)
    @test length(mixed_spec.energies) == 5
end

@testset "Tunable approximation guard" begin
    near_half_flux = TunableTransmon(:bad; EJmax = 20.0, EC = 0.25, flux = 0.5, asymmetry = 0.0, ncut = 6)
    @test_throws ArgumentError hamiltonian(CompositeSystem(near_half_flux))
end

@testset "Hamiltonian spec API" begin
    q1 = Transmon(:q1; EJ = 20.0, EC = 0.25, ng = 0.0, ncut = 6)
    q2 = Transmon(:q2; EJ = 19.0, EC = 0.24, ng = 0.1, ncut = 4)
    r1 = Resonator(:r1; ω = 6.8, dim = 3)
    circuit_spec = CircuitHamiltonianSpec(charge_cutoff = 2, charge_cutoffs = Dict(:q2 => 4))

    effective_model = build_model(CompositeSystem(q1, q2, r1))
    circuit_model = build_model(CompositeSystem(q1, q2, r1); hamiltonian_spec = circuit_spec)

    @test effective_model.hamiltonian_spec isa EffectiveHamiltonianSpec
    @test circuit_model.hamiltonian_spec isa CircuitHamiltonianSpec
    @test effective_model.dimensions[:q1] == 6
    @test circuit_model.dimensions[:q1] == 5
    @test circuit_model.dimensions[:q2] == 9
    @test circuit_model.dimensions[:r1] == 3
    @test length(basis_state(CompositeSystem(q1); hamiltonian_spec = CircuitHamiltonianSpec(charge_cutoff = 2), q1 = 2).data) == 5
    @test hamiltonian(CompositeSystem(q1)).data ≈ hamiltonian(CompositeSystem(q1); hamiltonian_spec = EffectiveHamiltonianSpec()).data

    @test_throws ArgumentError build_model(CompositeSystem(q1, q2); hamiltonian_spec = CircuitHamiltonianSpec(charge_cutoff = 2, charge_cutoffs = Dict(:missing => 3)))
    @test_throws ArgumentError build_model(CompositeSystem(q1, r1); hamiltonian_spec = CircuitHamiltonianSpec(charge_cutoff = 2, charge_cutoffs = Dict(:r1 => 3)))
    @test_throws ArgumentError build_model(CompositeSystem(q1, r1, CapacitiveCoupling(:q1, :r1; g = 0.02)); hamiltonian_spec = CircuitHamiltonianSpec(charge_cutoff = 2))
end

@testset "Circuit capacitive coupling" begin
    circuit_spec = CircuitHamiltonianSpec(charge_cutoff = 2)

    q1 = Transmon(:q1; EJ = 20.0, EC = 0.25, ng = 0.05, ncut = 5)
    q2 = TunableCoupler(:q2; EJmax = 15.0, EC = 0.30, flux = 0.08, asymmetry = 0.10, ng = 0.0, ncut = 5)
    qq_coupling = CircuitCapacitiveCoupling(:q1, :q2; G = 0.015)
    qq_uncoupled_model = build_model(CompositeSystem(q1, q2); hamiltonian_spec = circuit_spec)
    qq_model = build_model(CompositeSystem(q1, q2, qq_coupling); hamiltonian_spec = circuit_spec)
    qq_expected = hamiltonian(qq_uncoupled_model) + exact_charge_charge_circuit_coupling(qq_model, :q1, :q2, qq_coupling.G)
    @test hamiltonian(qq_model).data ≈ qq_expected.data atol = 1e-10

    q = Transmon(:q; EJ = 20.0, EC = 0.25, ng = 0.05, ncut = 5)
    r = Resonator(:r; ω = 6.2, dim = 3)
    qr_coupling = CircuitCapacitiveCoupling(:q, :r; G = 0.03)
    qr_uncoupled_model = build_model(CompositeSystem(q, r); hamiltonian_spec = circuit_spec)
    qr_model = build_model(CompositeSystem(q, r, qr_coupling); hamiltonian_spec = circuit_spec)
    qr_expected = hamiltonian(qr_uncoupled_model) + exact_resonator_charge_circuit_coupling(qr_model, :r, :q, qr_coupling.G)
    @test hamiltonian(qr_model).data ≈ qr_expected.data atol = 1e-10

    migrated_error =
        try
            build_model(CompositeSystem(q, r, CapacitiveCoupling(:q, :r; g = 0.02)); hamiltonian_spec = circuit_spec)
            nothing
        catch err
            err
        end
    @test migrated_error isa ArgumentError
    @test occursin("CircuitCapacitiveCoupling", sprint(showerror, migrated_error))

    effective_error =
        try
            build_model(CompositeSystem(q, r, qr_coupling))
            nothing
        catch err
            err
        end
    @test effective_error isa ArgumentError
    @test occursin("CircuitHamiltonianSpec", sprint(showerror, effective_error))

    r2 = Resonator(:r2; ω = 6.8, dim = 2)
    rr_error =
        try
            build_model(CompositeSystem(r, r2, CircuitCapacitiveCoupling(:r, :r2; G = 0.02)); hamiltonian_spec = circuit_spec)
            nothing
        catch err
            err
        end
    @test rr_error isa ArgumentError
    @test occursin("between resonators", sprint(showerror, rr_error))

    g_sweep = simulate_sweep(
        CompositeSystem(q, r, qr_coupling),
        SweepSpec(:q, :r, :G, [0.01, 0.05]; levels = 5);
        hamiltonian_spec = circuit_spec,
    )
    @test g_sweep.values == [0.01, 0.05]
    @test length(g_sweep.spectra) == 2
    @test g_sweep.spectra[1].energies != g_sweep.spectra[2].energies

    shifted_sys = with_coupling_parameter(CompositeSystem(q, r, qr_coupling), :q, :r, :G, 0.05)
    @test couplings(shifted_sys)[1].G == 0.05

    tq = TunableTransmon(:tq; EJmax = 20.0, EC = 0.25, flux = 0.12, asymmetry = 0.20, ng = 0.05, ncut = 6)
    mixed_spec = spectrum(
        CompositeSystem(tq, r, CircuitCapacitiveCoupling(:tq, :r; G = 0.025));
        levels = 5,
        hamiltonian_spec = CircuitHamiltonianSpec(charge_cutoff = 3),
    )
    @test length(mixed_spec.energies) == 5
end

@testset "Circuit Hamiltonian and operators" begin
    reference_spec = CircuitHamiltonianSpec(charge_cutoff = 3)
    operator_spec = CircuitHamiltonianSpec(charge_cutoff = 2)
    q_ng0 = Transmon(:q; EJ = 20.0, EC = 0.25, ng = 0.0, ncut = 6)
    q_ng1 = Transmon(:q; EJ = 20.0, EC = 0.25, ng = 0.25, ncut = 6)

    effective_spec0 = spectrum(CompositeSystem(q_ng0); levels = 4)
    effective_spec1 = spectrum(CompositeSystem(q_ng1); levels = 4)
    circuit_f01_ng0 = transition_01(CompositeSystem(q_ng0); hamiltonian_spec = reference_spec)
    circuit_f01_ng1 = transition_01(CompositeSystem(q_ng1); hamiltonian_spec = reference_spec)

    @test effective_spec0.energies ≈ effective_spec1.energies atol = 1e-10
    @test abs(circuit_f01_ng1 - circuit_f01_ng0) > 1e-2

    symmetric_fluxes = [0.0, 0.12, 0.24, 0.36]
    symmetric_f01 = [
        transition_01(
            CompositeSystem(TunableTransmon(:tq; EJmax = 20.0, EC = 0.25, flux = flux, asymmetry = 0.0, ng = 0.0, ncut = 6));
            hamiltonian_spec = reference_spec,
        ) for flux in symmetric_fluxes
    ]
    @test is_strictly_decreasing(symmetric_f01)

    convergence_tq = TunableTransmon(:tq_conv; EJmax = 20.0, EC = 0.25, flux = 0.12, asymmetry = 0.20, ng = 0.05, ncut = 6)
    cutoff_f01 = [
        transition_01(CompositeSystem(convergence_tq); hamiltonian_spec = CircuitHamiltonianSpec(charge_cutoff = cutoff))
        for cutoff in 3:5
    ]
    cutoff_corrections = abs.(diff(cutoff_f01))
    @test cutoff_corrections[2] < cutoff_corrections[1]
    @test abs(cutoff_f01[3] - cutoff_f01[2]) < 0.1

    resonator = Resonator(:r1; ω = 1.25, dim = 4)
    @test spectrum(CompositeSystem(resonator); levels = 4, hamiltonian_spec = operator_spec).energies ≈
        spectrum(CompositeSystem(resonator); levels = 4).energies atol = 1e-10

    tq_asym = TunableTransmon(:tq_asym; EJmax = 21.0, EC = 0.23, flux = 0.17, asymmetry = 0.35, ng = 0.15, ncut = 6)
    tq_asym_model = build_model(CompositeSystem(tq_asym); hamiltonian_spec = reference_spec)
    expected_tq_asym = exact_tunable_circuit_hamiltonian(tq_asym_model, :tq_asym, tq_asym)
    @test hamiltonian(tq_asym_model).data ≈ expected_tq_asym.data atol = 1e-10

    circuit_model = build_model(CompositeSystem(q_ng0); hamiltonian_spec = operator_spec)
    effective_model = build_model(CompositeSystem(q_ng0))
    resonator_model = build_model(CompositeSystem(resonator); hamiltonian_spec = operator_spec)

    @test size(charge_operator(circuit_model, :q).data) == (5, 5)
    @test size(cosphi_operator(circuit_model, :q).data) == (5, 5)
    @test size(sinphi_operator(circuit_model, :q).data) == (5, 5)
    @test_throws ArgumentError annihilation_operator(circuit_model, :q)
    @test_throws ArgumentError number_operator(circuit_model, :q)
    @test_throws ArgumentError charge_operator(effective_model, :q)
    @test size(number_operator(resonator_model, :r1).data) == (4, 4)
end

@testset "FluxControl validation and circuit dynamics" begin
    circuit_spec = CircuitHamiltonianSpec(charge_cutoff = 3)
    tq = TunableTransmon(:q; EJmax = 20.0, EC = 0.25, flux = 0.15, asymmetry = 0.25, ng = 0.05, ncut = 6)
    sys = CompositeSystem(tq)
    ψ0 = basis_state(sys; hamiltonian_spec = circuit_spec, q = 0)
    short_tlist = collect(range(0.0, 2.0; length = 21))
    zero_control = FluxControl(:zero_flux, :q, (p, t) -> 0.0)
    charge_observable = ObservableSpec(:charge, :q, :charge)

    baseline = evolve(
        sys,
        ψ0,
        short_tlist;
        hamiltonian_spec = circuit_spec,
        observables = [charge_observable],
    )
    zero_result = evolve(
        sys,
        ψ0,
        short_tlist;
        hamiltonian_spec = circuit_spec,
        observables = [charge_observable],
        flux_controls = [zero_control],
    )
    @test observable_trace(zero_result, :charge).values ≈ observable_trace(baseline, :charge).values atol = 1e-4

    shifted_flux = 0.07
    constant_control = FluxControl(:constant_flux, :q, (p, t) -> shifted_flux)
    shifted_sys = with_subsystem_parameter(sys, :q, :flux, tq.flux + shifted_flux)
    shifted_result = evolve(
        shifted_sys,
        basis_state(shifted_sys; hamiltonian_spec = circuit_spec, q = 0),
        short_tlist;
        hamiltonian_spec = circuit_spec,
        observables = [charge_observable],
    )
    constant_result = evolve(
        sys,
        ψ0,
        short_tlist;
        hamiltonian_spec = circuit_spec,
        observables = [charge_observable],
        flux_controls = [constant_control],
    )
    @test observable_trace(constant_result, :charge).values ≈ observable_trace(shifted_result, :charge).values atol = 1e-4

    flux_pulse = FluxControl(
        :flux_pulse,
        :q,
        (p, t) -> p.δ * sin(p.ω * t);
        derivative = (p, t) -> p.δ * p.ω * cos(p.ω * t),
    )
    long_tlist = collect(range(0.0, 16.0; length = 161))
    pulsed_result = evolve(
        sys,
        ψ0,
        long_tlist;
        hamiltonian_spec = circuit_spec,
        flux_controls = [flux_pulse],
        observables = [charge_observable],
        params = (; δ = 0.08, ω = 5.6),
    )
    pulsed_population = population_trace(pulsed_result, :q, 1)
    @test maximum(pulsed_population.values) > 1e-4

    qf = Transmon(:qf; EJ = 20.0, EC = 0.25, ncut = 6)
    @test_throws ArgumentError evolve(
        CompositeSystem(qf),
        basis_state(CompositeSystem(qf); hamiltonian_spec = circuit_spec, qf = 0),
        [0.0, 0.1];
        hamiltonian_spec = circuit_spec,
        flux_controls = [FluxControl(:bad_flux, :qf, (p, t) -> 0.0)],
    )

    @test_throws ArgumentError evolve(
        sys,
        ψ0,
        short_tlist;
        hamiltonian_spec = circuit_spec,
        flux_controls = [
            FluxControl(:dup_flux_a, :q, (p, t) -> 0.0),
            FluxControl(:dup_flux_b, :q, (p, t) -> 0.01),
        ],
    )

    other_tunable = TunableCoupler(:c; EJmax = 15.0, EC = 0.30, flux = 0.05, asymmetry = 0.15, ng = 0.0, ncut = 5)
    uncoupled_sys = CompositeSystem(tq, other_tunable)
    uncoupled_ψ0 = basis_state(uncoupled_sys; hamiltonian_spec = circuit_spec, q = 0, c = 0)
    @test_throws ArgumentError evolve(
        uncoupled_sys,
        uncoupled_ψ0,
        [0.0, 0.1];
        hamiltonian_spec = circuit_spec,
        flux_controls = [
            FluxControl(:dup_label, :q, (p, t) -> 0.0),
            FluxControl(:dup_label, :c, (p, t) -> 0.0),
        ],
    )

    @test_throws ArgumentError evolve(
        sys,
        basis_state(sys; q = 0),
        short_tlist;
        flux_controls = [FluxControl(:effective_bad, :q, (p, t) -> 0.01)],
    )
end

@testset "Coupled circuit capacitive dynamics and flux control" begin
    drive_spec = CircuitHamiltonianSpec(charge_cutoff = 2)
    q = Transmon(:q; EJ = 20.0, EC = 0.25, ng = 0.05, ncut = 5)
    r = Resonator(:r; ω = 6.2, dim = 3)
    drive_sys = CompositeSystem(q, r, CircuitCapacitiveCoupling(:q, :r; G = 0.03))
    ψ0_drive = basis_state(drive_sys; hamiltonian_spec = drive_spec, q = 0, r = 0)
    drive_tlist = collect(range(0.0, 8.0; length = 161))
    drive_observables = [ObservableSpec(:nr, :r, :n), ObservableSpec(:charge_q, :q, :charge)]
    baseline = evolve(
        drive_sys,
        ψ0_drive,
        drive_tlist;
        hamiltonian_spec = drive_spec,
        observables = drive_observables,
    )
    drive = SubsystemDrive(
        :r_drive,
        :r,
        :x,
        (p, t) -> p.Ω * cos(p.ωd * t),
    )
    driven = evolve(
        drive_sys,
        ψ0_drive,
        drive_tlist;
        hamiltonian_spec = drive_spec,
        drives = [drive],
        observables = drive_observables,
        params = (; Ω = 0.18, ωd = 6.2),
    )
    @test maximum(abs.(real.(observable_trace(driven, :charge_q).values) .- real.(observable_trace(baseline, :charge_q).values))) > 1e-2
    @test maximum(real.(observable_trace(driven, :nr).values)) > 1e-2
    @test maximum(population_trace(driven, :q, 1).values) > 1e-2

    flux_spec = CircuitHamiltonianSpec(charge_cutoff = 3)
    tq = TunableTransmon(:q; EJmax = 20.0, EC = 0.25, flux = 0.12, asymmetry = 0.20, ng = 0.05, ncut = 6)
    flux_resonator = Resonator(:r; ω = 6.0, dim = 2)
    flux_sys = CompositeSystem(tq, flux_resonator, CircuitCapacitiveCoupling(:q, :r; G = 0.025))
    ψ0_flux = basis_state(flux_sys; hamiltonian_spec = flux_spec, q = 0, r = 0)
    short_tlist = collect(range(0.0, 2.0; length = 21))
    flux_observables = [ObservableSpec(:charge_q, :q, :charge), ObservableSpec(:nr, :r, :n)]
    zero_control = FluxControl(:zero_flux, :q, (p, t) -> 0.0)

    flux_baseline = evolve(
        flux_sys,
        ψ0_flux,
        short_tlist;
        hamiltonian_spec = flux_spec,
        observables = flux_observables,
    )
    zero_result = evolve(
        flux_sys,
        ψ0_flux,
        short_tlist;
        hamiltonian_spec = flux_spec,
        observables = flux_observables,
        flux_controls = [zero_control],
    )
    @test observable_trace(zero_result, :charge_q).values ≈ observable_trace(flux_baseline, :charge_q).values atol = 1e-4
    @test observable_trace(zero_result, :nr).values ≈ observable_trace(flux_baseline, :nr).values atol = 1e-4

    constant_shift = 0.05
    constant_control = FluxControl(:constant_flux, :q, (p, t) -> constant_shift)
    shifted_sys = with_subsystem_parameter(flux_sys, :q, :flux, tq.flux + constant_shift)
    shifted_result = evolve(
        shifted_sys,
        basis_state(shifted_sys; hamiltonian_spec = flux_spec, q = 0, r = 0),
        short_tlist;
        hamiltonian_spec = flux_spec,
        observables = flux_observables,
    )
    constant_result = evolve(
        flux_sys,
        ψ0_flux,
        short_tlist;
        hamiltonian_spec = flux_spec,
        observables = flux_observables,
        flux_controls = [constant_control],
    )
    @test observable_trace(constant_result, :charge_q).values ≈ observable_trace(shifted_result, :charge_q).values atol = 1e-4
    @test observable_trace(constant_result, :nr).values ≈ observable_trace(shifted_result, :nr).values atol = 1e-4

    flux_pulse = FluxControl(
        :flux_pulse,
        :q,
        (p, t) -> p.δ * sin(p.ω * t);
        derivative = (p, t) -> p.δ * p.ω * cos(p.ω * t),
    )
    long_tlist = collect(range(0.0, 12.0; length = 161))
    pulsed_result = evolve(
        flux_sys,
        ψ0_flux,
        long_tlist;
        hamiltonian_spec = flux_spec,
        observables = flux_observables,
        flux_controls = [flux_pulse],
        params = (; δ = 0.08, ω = 5.6),
    )
    @test maximum(abs.(real.(observable_trace(pulsed_result, :charge_q).values))) > 1e-4
    @test maximum(real.(observable_trace(pulsed_result, :nr).values)) > 1e-4
    @test maximum(population_trace(pulsed_result, :q, 1).values) > 1e-4
end

@testset "Circuit sweeps and dynamics" begin
    circuit_spec = CircuitHamiltonianSpec(charge_cutoff = 2)
    tq = TunableTransmon(:q; EJmax = 20.0, EC = 0.25, flux = 0.0, asymmetry = 0.0, ng = 0.0, ncut = 6)
    q = Transmon(:q; EJ = 20.0, EC = 0.25, ng = 0.0, ncut = 6)

    flux_result = simulate_sweep(CompositeSystem(tq), SweepSpec(:q, :flux, [0.0, 0.15]; levels = 4); hamiltonian_spec = circuit_spec)
    ej_result = simulate_sweep(CompositeSystem(q), SweepSpec(:q, :EJ, [18.0, 22.0]; levels = 4); hamiltonian_spec = circuit_spec)
    ec_result = simulate_sweep(CompositeSystem(q), SweepSpec(:q, :EC, [0.20, 0.30]; levels = 4); hamiltonian_spec = circuit_spec)
    ng_result = simulate_sweep(CompositeSystem(q), SweepSpec(:q, :ng, [0.0, 0.25]; levels = 4); hamiltonian_spec = circuit_spec)

    @test flux_result.values == [0.0, 0.15]
    @test length(flux_result.spectra) == 2
    @test transition_frequencies(flux_result.spectra[2])[1] < transition_frequencies(flux_result.spectra[1])[1]
    @test ej_result.spectra[1].energies != ej_result.spectra[2].energies
    @test ec_result.spectra[1].energies != ec_result.spectra[2].energies
    @test ng_result.spectra[1].energies != ng_result.spectra[2].energies

    ψ0 = basis_state(CompositeSystem(tq); hamiltonian_spec = circuit_spec, q = 2)
    tlist = collect(range(0.0, 1.0; length = 11))
    drive = SubsystemDrive(
        :q_drive,
        :q,
        :sinphi,
        (p, t) -> p.Ω * cos(p.ωd * t),
    )

    result = evolve(
        CompositeSystem(tq),
        ψ0,
        tlist;
        hamiltonian_spec = circuit_spec,
        drives = [drive],
        observables = [ObservableSpec(:charge, :q, :charge), ObservableSpec(:phi, :q, :cosphi)],
        params = (; Ω = 0.2, ωd = 1.0),
    )

    charge_trace = observable_trace(result, :charge)
    phi_trace = observable_trace(result, :phi)

    @test result.model.hamiltonian_spec isa CircuitHamiltonianSpec
    @test result.times == tlist
    @test length(result.states) == length(tlist)
    @test charge_trace.times == tlist
    @test phi_trace.times == tlist
    @test maximum(abs.(real.(charge_trace.values))) > 1e-4
    @test maximum(abs.(real.(phi_trace.values))) > 1e-2
    @test all(value -> isapprox(imag(value), 0.0; atol = 1e-8), charge_trace.values)
    @test_throws ArgumentError evolve(CompositeSystem(tq), ψ0, tlist; hamiltonian_spec = circuit_spec, observables = [ObservableSpec(:bad, :q, :x)])
end

@testset "Sweep utilities" begin
    tq = TunableTransmon(:tq; EJmax = 20.0, EC = 0.25, flux = 0.0, asymmetry = 0.0, ncut = 6)
    tunable_sys = CompositeSystem(tq)

    flux_sweep = SweepSpec(:tq, :flux, [0.0, 0.15, 0.30]; levels = 4)
    flux_result = simulate_sweep(tunable_sys, flux_sweep)

    @test flux_result.values == [0.0, 0.15, 0.30]
    @test length(flux_result.spectra) == 3
    @test transition_frequencies(flux_result.spectra[1])[1] > transition_frequencies(flux_result.spectra[end])[1]
    @test flux_result.spectra[2].model.system.subsystems[1].flux == 0.15

    q1 = Transmon(:q1; EJ = 20.0, EC = 0.25, ncut = 5)
    r1 = Resonator(:r1; ω = 6.8, dim = 4)
    g1 = CapacitiveCoupling(:q1, :r1; g = 0.02)
    coupled_sys = CompositeSystem(q1, r1, g1)

    coupling_sweep = SweepSpec(:q1, :r1, :g, [0.02, 0.08, 0.14]; levels = 5)
    coupling_result = simulate_sweep(coupled_sys, coupling_sweep)

    @test length(coupling_result.spectra) == 3
    @test coupling_result.spectra[3].model.system.couplings[1].g == 0.14
    @test coupling_result.spectra[1].energies != coupling_result.spectra[3].energies

    @test_throws ArgumentError SweepSpec(:tq, :flux, Float64[]; levels = 4)
    @test_throws ArgumentError simulate_sweep(tunable_sys, SweepSpec(:missing, :flux, [0.0, 0.1]))
end

@testset "Sweep analysis" begin
    tq = TunableTransmon(:tq; EJmax = 20.0, EC = 0.25, flux = 0.0, asymmetry = 0.0, ncut = 6)
    result = simulate_sweep(CompositeSystem(tq), SweepSpec(:tq, :flux, [0.0, 0.15, 0.30]; levels = 4))

    transition01 = transition_curve(result)
    @test transition01.values == [0.0, 0.15, 0.30]
    @test transition01.label == :transition_01
    @test transition01.data[1] ≈ transition_frequencies(result.spectra[1])[1] atol = 1e-10
    @test transition01.data[end] < transition01.data[1]

    anh_curve = anharmonicity_curve(result)
    @test anh_curve.label == :anharmonicity
    @test all(value -> isapprox(value, -tq.EC; atol = 1e-10), anh_curve.data)

    summary = sweep_summary(result)
    @test summary.values == transition01.values
    @test summary.transition_01 == transition01.data
    @test summary.anharmonicity == anh_curve.data

    gap = minimum_gap(result)
    @test gap.gap == minimum(transition01.data)
    @test gap.sweep_value == result.values[gap.index]
    @test gap.level_pair == (1, 2)

    @test_throws ArgumentError transition_curve(result; transition = (2, 2))
    @test_throws ArgumentError minimum_gap(result; level_pair = (1, 5))
end

@testset "Eigensystem simulation" begin
    q = Transmon(:q; EJ = 20.0, EC = 0.25, ncut = 5)
    system = CompositeSystem(q)
    esys = eigensystem(system; levels = 4)
    spec = spectrum(system; levels = 4)

    @test length(esys.energies) == 4
    @test length(esys.states) == 4
    @test issorted(esys.energies)
    @test esys.energies ≈ spec.energies atol = 1e-10

    tq = TunableTransmon(:tq; EJmax = 20.0, EC = 0.25, flux = 0.0, asymmetry = 0.0, ncut = 6)
    sweep_spec = SweepSpec(:tq, :flux, [0.0, 0.15, 0.30]; levels = 4)
    energy_sweep = simulate_sweep(CompositeSystem(tq), sweep_spec)
    eigensystem_sweep = simulate_eigensystem_sweep(CompositeSystem(tq), sweep_spec)

    @test eigensystem_sweep.values == energy_sweep.values
    @test length(eigensystem_sweep.spectra) == length(energy_sweep.spectra)
    @test length(eigensystem_sweep.spectra[1].states) == 4
    @test transition_curve(eigensystem_sweep).data ≈ transition_curve(energy_sweep).data atol = 1e-10
    @test minimum_gap(eigensystem_sweep).gap ≈ minimum_gap(energy_sweep).gap atol = 1e-10
    @test minimum_gap(eigensystem_sweep).sweep_value == minimum_gap(energy_sweep).sweep_value
end

@testset "Bare-label lookup" begin
    q = Transmon(:q; EJ = 20.0, EC = 0.25, ncut = 4)
    r = Resonator(:r; ω = 6.05, dim = 3)
    coupled_sys = CompositeSystem(q, r, CapacitiveCoupling(:q, :r; g = 0.02))
    sweep = simulate_eigensystem_sweep(
        coupled_sys,
        SweepSpec(:q, :r, :g, [0.02, 0.60]; levels = 6),
    )
    basis = bare_product_basis(coupled_sys; subsystem_levels = Dict(:q => 3, :r => 2))
    labeled = label_sweep(sweep, basis)

    @test length(basis.subsystems) == 2
    @test length(basis.labels) == 6
    @test length(basis.states) == 6
    @test length(labeled.label_maps) == 2

    assigned = dressed_index(labeled, (1, 0), 1)
    @test !ismissing(assigned)
    @test bare_label(labeled, assigned, 1) == (1, 0)
    @test energy_by_bare_label(labeled, (1, 0), 1) ≈ sweep.spectra[1].energies[assigned] atol = 1e-10

    curve = energy_curve(labeled, (1, 0))
    @test curve.values == sweep.values
    @test curve.label == :bare_1_0
    @test length(curve.data) == length(sweep.values)
    @test curve.data[1] ≈ energy_by_bare_label(labeled, (1, 0), 1) atol = 1e-10

    components = dressed_state_components(labeled, (1, 0), 1)
    probabilities = last.(components)
    @test !isempty(components)
    @test issorted(probabilities; rev = true)
    @test sum(probabilities) ≈ 1.0 atol = 1e-8
    @test length(dressed_state_components(labeled, (1, 0), 1; top = 2)) == 2

    strict_labeled = label_sweep(sweep, basis; overlap_threshold = 0.99)
    @test count(ismissing, strict_labeled.label_maps[end].bare_to_dressed) > 0
    missing_position = findfirst(ismissing, strict_labeled.label_maps[end].bare_to_dressed)
    @test !isnothing(missing_position)
    missing_curve = energy_curve(strict_labeled, strict_labeled.basis.labels[missing_position])
    @test isnan(missing_curve.data[end])
end

@testset "README quickstart examples" begin
    q1 = Transmon(:q1; EJ = 20.0, EC = 0.25, ncut = 6)
    static_sys = CompositeSystem(q1)
    static_spec = spectrum(static_sys; levels = 4)
    static_energies = static_spec.energies
    static_freqs = transition_frequencies(static_spec)
    static_alpha = anharmonicity(static_spec)
    static_view = (; energies = static_energies, ω01 = static_freqs[1], α = static_alpha)

    @test length(static_view.energies) == 4
    @test static_view.ω01 ≈ static_freqs[1] atol = 1e-10
    @test static_view.α ≈ -q1.EC atol = 1e-10

    tq = TunableTransmon(:tq; EJmax = 20.0, EC = 0.25, flux = 0.0, asymmetry = 0.0, ncut = 6)
    tunable_sys = CompositeSystem(tq)
    tunable_sweep = SweepSpec(:tq, :flux, [0.0, 0.15, 0.30]; levels = 4)
    tunable_result = simulate_sweep(tunable_sys, tunable_sweep)
    ω01_curve = transition_curve(tunable_result)
    alpha_curve = anharmonicity_curve(tunable_result)
    tunable_gap = minimum_gap(tunable_result)
    tunable_summary = sweep_summary(tunable_result)
    tunable_view = (; flux_values = tunable_summary.values, ω01 = ω01_curve.data, α = alpha_curve.data, minimum_gap = tunable_gap)

    @test tunable_view.flux_values == [0.0, 0.15, 0.30]
    @test length(tunable_view.ω01) == 3
    @test tunable_view.minimum_gap.gap == minimum(tunable_view.ω01)

    q_fixed = Transmon(:qf; EJ = 19.0, EC = 0.24, ncut = 5)
    coupler = TunableCoupler(:c1; EJmax = 15.0, EC = 0.30, flux = 0.0, asymmetry = 0.1, ncut = 4)
    resonator = Resonator(:r1; ω = 6.7, dim = 3)
    mixed_sys = CompositeSystem(
        q_fixed,
        coupler,
        resonator,
        CapacitiveCoupling(:qf, :c1; g = 0.07),
        CapacitiveCoupling(:c1, :r1; g = 0.05),
    )
    mixed_spec = spectrum(mixed_sys; levels = 5)
    mixed_result = simulate_sweep(mixed_sys, SweepSpec(:c1, :flux, [0.0, 0.10, 0.20]; levels = 5))
    mixed_summary = sweep_summary(mixed_result)
    mixed_view = (; low_lying_energies = mixed_spec.energies, coupler_flux_values = mixed_summary.values, transition_01 = mixed_summary.transition_01)

    @test length(mixed_view.low_lying_energies) == 5
    @test mixed_view.coupler_flux_values == [0.0, 0.10, 0.20]
    @test length(mixed_view.transition_01) == 3

    q_static = Transmon(:q1; EJ = 20.0, EC = 0.25, ncut = 5)
    r_static = Resonator(:r1; ω = 6.8, dim = 4)
    coupled_sys = CompositeSystem(
        q_static,
        r_static,
        CapacitiveCoupling(:q1, :r1; g = 0.02),
    )
    coupling_result = simulate_sweep(coupled_sys, SweepSpec(:q1, :r1, :g, [0.02, 0.08, 0.14]; levels = 5))
    coupling_curve = transition_curve(coupling_result)
    coupling_summary = sweep_summary(coupling_result)
    coupling_view = (; coupling_values = coupling_summary.values, ω01 = coupling_curve.data, weakest_coupling_spectrum = coupling_result.spectra[1].energies)

    @test coupling_view.coupling_values == [0.02, 0.08, 0.14]
    @test length(coupling_view.ω01) == 3
    @test length(coupling_view.weakest_coupling_spectrum) == 5

    resonator = Resonator(:r1; ω = 1.0, dim = 2)
    dynamics_sys = CompositeSystem(resonator)
    ψ0 = basis_state(dynamics_sys; r1 = 0)
    tlist = collect(range(0.0, 6.0; length = 61))
    drive = SubsystemDrive(
        :r1_x_drive,
        :r1,
        :x,
        (p, t) -> p.Ω * cos(p.ωd * t),
    )
    dynamics_result = evolve(
        dynamics_sys,
        ψ0,
        tlist;
        drives = [drive],
        observables = [ObservableSpec(:nr, :r1, :n)],
        params = (; Ω = 0.35, ωd = 1.0),
    )
    dynamics_trace = observable_trace(dynamics_result, :nr)
    ground_population = population_trace(dynamics_result, :r1, 0)
    excited_population = population_trace(dynamics_result, :r1, 1)
    dynamics_view = (
        ;
        times = dynamics_trace.times,
        driven_population = real.(dynamics_trace.values),
        ground_population = ground_population.values,
        excited_population = excited_population.values,
        final_state = final_state(dynamics_result),
    )

    @test dynamics_view.times == tlist
    @test length(dynamics_view.driven_population) == length(tlist)
    @test length(dynamics_view.ground_population) == length(tlist)
    @test length(dynamics_view.excited_population) == length(tlist)
    @test maximum(abs.(dynamics_view.driven_population)) > 1e-3
    @test maximum(dynamics_view.excited_population) > 1e-3
end

@testset "Phase 3A basis states" begin
    q1 = Transmon(:q1; EJ = 20.0, EC = 0.25, ncut = 3)
    r1 = Resonator(:r1; ω = 6.8, dim = 2)
    sys = CompositeSystem(q1, r1, CapacitiveCoupling(:q1, :r1; g = 0.02))
    model = build_model(sys)

    ψ0 = basis_state(sys; q1 = 1)
    ψ1 = basis_state(model; q1 = 1, r1 = 1)

    @test length(ψ0.data) == 6
    @test length(ψ1.data) == 6
    @test ψ0.data[3] ≈ 1.0 + 0.0im atol = 1e-10
    @test ψ1.data[4] ≈ 1.0 + 0.0im atol = 1e-10

    @test_throws ArgumentError basis_state(sys; missing = 1)
    @test_throws ArgumentError basis_state(sys; q1 = 3)
end

@testset "Phase 3A closed-system dynamics" begin
    resonator = Resonator(:r1; ω = 6.0, dim = 2)
    sys = CompositeSystem(resonator)
    model = build_model(sys)
    ψ0 = basis_state(sys; r1 = 1)
    tlist = collect(range(0.0, 1.0; length = 11))

    result = evolve(
        sys,
        ψ0,
        tlist;
        observables = [ObservableSpec(:nr, :r1, :n)],
    )

    nr_trace = observable_trace(result, :nr)
    p0_trace = population_trace(result, :r1, 0)
    p1_trace = population_trace(result, :r1, 1)
    @test result.model.dimensions[:r1] == 2
    @test result.times == tlist
    @test length(result.states) == length(tlist)
    @test nr_trace.times == tlist
    @test length(nr_trace.values) == length(tlist)
    @test p0_trace.times == tlist
    @test p1_trace.times == tlist
    @test all(value -> isapprox(real(value), 1.0; atol = 1e-6), nr_trace.values)
    @test all(value -> isapprox(imag(value), 0.0; atol = 1e-6), nr_trace.values)
    @test all(value -> isapprox(value, 0.0; atol = 1e-6), p0_trace.values)
    @test all(value -> isapprox(value, 1.0; atol = 1e-6), p1_trace.values)
    @test final_state(result) == result.states[end]
    @test_throws ArgumentError observable_trace(result, :missing)
    @test_throws ArgumentError population_trace(result, :missing, 0)
    @test_throws ArgumentError population_trace(result, :r1, 2)
    @test_throws ArgumentError evolve(sys, ψ0, tlist; observables = [ObservableSpec(:nr, :r1, :n)], saveat = [0.0, 1.0])

    manual_expectation = [expect(number_operator(model, :r1), state) for state in result.states]
    @test nr_trace.values ≈ manual_expectation atol = 1e-6
end

@testset "Phase 3A driven closed-system dynamics" begin
    resonator = Resonator(:r1; ω = 1.0, dim = 2)
    sys = CompositeSystem(resonator)
    ψ0 = basis_state(sys; r1 = 0)
    tlist = collect(range(0.0, 6.0; length = 61))

    drive = SubsystemDrive(
        :r1_x_drive,
        :r1,
        :x,
        (p, t) -> p.Ω * cos(p.ωd * t),
    )

    result = evolve(
        sys,
        ψ0,
        tlist;
        drives = [drive],
        observables = [ObservableSpec(:nr, :r1, :n)],
        params = (; Ω = 0.35, ωd = 1.0),
    )

    nr_trace = observable_trace(result, :nr)
    p0_trace = population_trace(result, :r1, 0)
    p1_trace = population_trace(result, :r1, 1)
    @test length(result.states) == length(tlist)
    @test length(nr_trace.values) == length(tlist)
    @test length(p0_trace.values) == length(tlist)
    @test length(p1_trace.values) == length(tlist)
    @test maximum(abs.(real.(nr_trace.values))) > 1e-3
    @test maximum(abs.(imag.(nr_trace.values))) < 1e-7
    @test maximum(p1_trace.values) > 1e-3
    @test all(value -> 0.0 <= value <= 1.0 + 1e-6, p0_trace.values)
    @test all(value -> 0.0 <= value <= 1.0 + 1e-6, p1_trace.values)
    @test p0_trace.values .+ p1_trace.values ≈ ones(length(tlist)) atol = 1e-6
    @test real(nr_trace.values[end]) != real(nr_trace.values[1])
end

@testset "Nonadiabatic effective flux dynamics" begin
    nonadiabatic_spec = EffectiveHamiltonianSpec(NonadiabaticDuffingEffectiveMethod())
    tq = TunableTransmon(:q; EJmax = 20.0, EC = 0.25, flux = 0.15, asymmetry = 0.25, ncut = 4)
    sys = CompositeSystem(tq)
    ψ0 = basis_state(sys; hamiltonian_spec = nonadiabatic_spec, q = 0)
    tlist = collect(range(0.0, 24.0; length = 181))
    flux_pulse = FluxControl(
        :flux_drive,
        :q,
        (p, t) -> p.δ * sin(p.ωf * t);
        derivative = (p, t) -> p.δ * p.ωf * cos(p.ωf * t),
    )
    drive = SubsystemDrive(
        :q_x_drive,
        :q,
        :x,
        (p, t) -> p.Ω * cos(p.ωd * t),
    )

    flux_only = evolve(
        sys,
        ψ0,
        tlist;
        hamiltonian_spec = nonadiabatic_spec,
        flux_controls = [flux_pulse],
        observables = [ObservableSpec(:nq, :q, :n)],
        params = (; δ = 0.06, ωf = 6.0),
    )
    mixed_result = evolve(
        sys,
        ψ0,
        tlist;
        hamiltonian_spec = nonadiabatic_spec,
        drives = [drive],
        flux_controls = [flux_pulse],
        observables = [ObservableSpec(:nq, :q, :n)],
        params = (; δ = 0.06, ωf = 6.0, Ω = 0.02, ωd = 6.0),
    )
    flux_population = population_trace(flux_only, :q, 1)
    mixed_population = population_trace(mixed_result, :q, 1)
    @test maximum(flux_population.values) > 1e-4
    @test maximum(abs.(real.(observable_trace(mixed_result, :nq).values) .- real.(observable_trace(flux_only, :nq).values))) > 1e-4

    constant_shift = 0.04
    constant_control = FluxControl(:constant_flux, :q, (p, t) -> constant_shift; derivative = (p, t) -> 0.0)
    shifted_sys = with_subsystem_parameter(sys, :q, :flux, tq.flux + constant_shift)
    shifted_result = evolve(
        shifted_sys,
        basis_state(shifted_sys; hamiltonian_spec = nonadiabatic_spec, q = 0),
        tlist;
        hamiltonian_spec = nonadiabatic_spec,
        observables = [ObservableSpec(:nq, :q, :n)],
    )
    constant_result = evolve(
        sys,
        ψ0,
        tlist;
        hamiltonian_spec = nonadiabatic_spec,
        flux_controls = [constant_control],
        observables = [ObservableSpec(:nq, :q, :n)],
    )
    @test observable_trace(constant_result, :nq).values ≈ observable_trace(shifted_result, :nq).values atol = 1e-7

    @test_throws ArgumentError evolve(
        sys,
        ψ0,
        [0.0, 0.1];
        hamiltonian_spec = nonadiabatic_spec,
        flux_controls = [FluxControl(:missing_derivative, :q, (p, t) -> 0.01 * sin(t))],
    )

    resonator = Resonator(:r1; ω = 6.0, dim = 2)
    coupled_sys = CompositeSystem(tq, resonator, CapacitiveCoupling(:q, :r1; g = 0.02))
    coupled_ψ0 = basis_state(coupled_sys; hamiltonian_spec = nonadiabatic_spec, q = 0, r1 = 0)
    @test_throws ArgumentError evolve(
        coupled_sys,
        coupled_ψ0,
        [0.0, 0.1];
        hamiltonian_spec = nonadiabatic_spec,
        flux_controls = [flux_pulse],
        params = (; δ = 0.01, ωf = 1.0),
    )

    other_tunable = TunableCoupler(:c; EJmax = 15.0, EC = 0.30, flux = 0.05, asymmetry = 0.10, ncut = 3)
    uncoupled_multi_sys = CompositeSystem(tq, other_tunable)
    uncoupled_multi_ψ0 = basis_state(uncoupled_multi_sys; hamiltonian_spec = nonadiabatic_spec, q = 0, c = 0)
    @test_throws ArgumentError evolve(
        uncoupled_multi_sys,
        uncoupled_multi_ψ0,
        [0.0, 0.1];
        hamiltonian_spec = nonadiabatic_spec,
        flux_controls = [flux_pulse],
        params = (; δ = 0.01, ωf = 1.0),
    )
end

@testset "Phase 3A tunable-coupler two-qubit dynamics" begin
    tq1 = TunableTransmon(:tq1; EJmax = 20.0, EC = 0.25, flux = 0.0, asymmetry = 0.05, ncut = 3)
    c1 = TunableCoupler(:c1; EJmax = 15.0, EC = 0.30, flux = 0.0, asymmetry = 0.10, ncut = 3)
    q2 = Transmon(:q2; EJ = 19.0, EC = 0.24, ncut = 3)
    sys = CompositeSystem(
        tq1,
        c1,
        q2,
        CapacitiveCoupling(:tq1, :c1; g = 0.10),
        CapacitiveCoupling(:c1, :q2; g = 0.10),
    )

    spec = spectrum(sys; levels = 5)
    flux_sweep = simulate_sweep(sys, SweepSpec(:c1, :flux, [0.0, 0.12]; levels = 5))

    ψ0 = basis_state(sys; tq1 = 0, c1 = 0, q2 = 0)
    tlist = collect(range(0.0, 16.0; length = 121))
    drive = SubsystemDrive(
        :tq1_x_drive,
        :tq1,
        :x,
        (p, t) -> p.Ω * cos(p.ωd * t),
    )

    baseline = evolve(
        sys,
        ψ0,
        tlist;
        observables = [ObservableSpec(:ntq1, :tq1, :n), ObservableSpec(:nq2, :q2, :n)],
    )
    driven = evolve(
        sys,
        ψ0,
        tlist;
        drives = [drive],
        observables = [ObservableSpec(:ntq1, :tq1, :n), ObservableSpec(:nq2, :q2, :n)],
        params = (; Ω = 0.25, ωd = 5.95),
    )

    tq1_baseline = population_trace(baseline, :tq1, 1)
    tq1_population = population_trace(driven, :tq1, 1)
    c1_population = population_trace(driven, :c1, 1)
    q2_population = population_trace(driven, :q2, 1)
    tq1_number = observable_trace(driven, :ntq1)
    q2_number = observable_trace(driven, :nq2)

    @test length(spec.energies) == 5
    @test flux_sweep.values == [0.0, 0.12]
    @test length(flux_sweep.spectra) == 2
    @test tq1_population.times == tlist
    @test c1_population.times == tlist
    @test q2_population.times == tlist
    @test length(tq1_number.values) == length(tlist)
    @test length(q2_number.values) == length(tlist)
    @test maximum(tq1_baseline.values) ≈ 0.0 atol = 1e-8
    @test maximum(tq1_population.values) > 1e-2
    @test maximum(c1_population.values) > 1e-3
    @test maximum(q2_population.values) > 1e-3
    @test maximum(real.(tq1_number.values)) > 1e-2
    @test maximum(real.(q2_number.values)) > 1e-3
end

@testset "Projected unitary and subspace analysis" begin
    identity4 = ComplexF64[
        1 0 0 0
        0 1 0 0
        0 0 1 0
        0 0 0 1
    ]

    phase_sys = CompositeSystem(
        Resonator(:r1; ω = 1.0, dim = 2),
        Resonator(:r2; ω = 1.5, dim = 2),
    )
    phase_spec = subspace_spec(phase_sys; subsystem_levels = (r1 = [0, 1], r2 = [0, 1]))
    phase_trace = projected_unitary(phase_sys, phase_spec, [0.0, 1.0])
    local_phase = ComplexF64[
        exp(0.4im) 0 0 0
        0 exp(0.1im) 0 0
        0 0 exp(-0.2im) 0
        0 0 0 exp(-0.5im)
    ]

    @test phase_spec.labels == ["|00>", "|01>", "|10>", "|11>"]
    @test phase_trace.unitaries[1] ≈ identity4 atol = 1e-12
    @test maximum(phase_trace.leakages) ≈ 0.0 atol = 1e-12
    @test conditional_phase(phase_trace.unitaries[end]) ≈ 0.0 atol = 1e-6
    @test strip_local_z_phases(phase_trace.unitaries[end]) ≈ identity4 atol = 1e-6
    @test conditional_phase(phase_trace.unitaries[end] * local_phase) ≈ conditional_phase(phase_trace.unitaries[end]) atol = 1e-6

    dressed_sys = CompositeSystem(
        Resonator(:a; ω = 1.0, dim = 2),
        Resonator(:b; ω = 1.2, dim = 2),
        CapacitiveCoupling(:a, :b; g = 0.02),
    )
    dressed_spec = subspace_spec(dressed_sys; subsystem_levels = (a = [0, 1], b = [0, 1]), basis = :dressed_static)
    dressed_trace = projected_unitary(dressed_sys, dressed_spec, [0.0, 1.0])

    @test minimum(dressed_spec.reference_overlaps) > 0.98
    @test dressed_trace.unitaries[1] ≈ identity4 atol = 1e-12

    leakage_sys = CompositeSystem(Resonator(:r1; ω = 1.0, dim = 3))
    leakage_spec = subspace_spec(leakage_sys; subsystem_levels = (r1 = [0, 1],))
    leakage_drive = SubsystemDrive(
        :r1_x_drive,
        :r1,
        :x,
        (p, t) -> p.Ω * cos(p.ωd * t),
    )
    leakage_trace = projected_unitary(
        leakage_sys,
        leakage_spec,
        collect(range(0.0, 4.0; length = 41));
        drives = [leakage_drive],
        params = (; Ω = 0.35, ωd = 1.0),
    )
    leakage_move = best_move(leakage_trace; source_index = 1, target_index = 2)

    @test maximum(leakage_trace.leakages) > 1e-2
    @test leakage_move.transfer_probability > 0.3
    @test leakage_move.leakage > 1e-2
end

@testset "Spectrum comparison helpers" begin
    snapshot = load_renger2026_snapshot()
    pair = renger2026_model_pair(snapshot)
    comparison = compare_model_spectra(
        pair.circuit.system,
        pair.effective.system;
        reference_hamiltonian_spec = pair.circuit.hamiltonian_spec,
        candidate_hamiltonian_spec = pair.effective.hamiltonian_spec,
        levels = 6,
    )

    @test length(comparison.level_deltas) == 6
    @test length(comparison.transition_deltas) == 5
    @test all(isfinite, comparison.level_deltas)
    @test all(isfinite, comparison.transition_deltas)

    q = TunableTransmon(:q; EJmax = 2.0, EC = 0.2, flux = 0.0, asymmetry = 0.1, ncut = 4)
    r = Resonator(:r; ω = 1.2, dim = 2)
    values = collect(range(0.0, 0.5; length = 41))
    effective_sweep = simulate_sweep(
        CompositeSystem(q, r, CapacitiveCoupling(:q, :r; g = 0.03)),
        SweepSpec(:q, :flux, values; levels = 4),
    )
    circuit_sweep = simulate_sweep(
        CompositeSystem(q, r, CircuitCapacitiveCoupling(:q, :r; G = 0.03)),
        SweepSpec(:q, :flux, values; levels = 4);
        hamiltonian_spec = CircuitHamiltonianSpec(charge_cutoff = 2),
    )
    gap_comparison = compare_minimum_gap(circuit_sweep, effective_sweep)

    @test gap_comparison.reference.sweep_value == values[end]
    @test gap_comparison.candidate.sweep_value == values[end]
    @test gap_comparison.sweep_value_delta == 0.0
    @test isfinite(gap_comparison.gap_delta)
end

@testset "Renger 2026 workflow builders" begin
    snapshot = load_renger2026_snapshot()
    qr_stage = renger2026_stage1_qr_system(snapshot)
    qcr_stage = renger2026_stage1_qcr_system(snapshot)
    pair = renger2026_model_pair(snapshot)

    @test length(subsystems(qr_stage.system)) == 2
    @test length(couplings(qr_stage.system)) == 1
    @test length(subsystems(qcr_stage.system)) == 3
    @test length(couplings(qcr_stage.system)) == 2
    @test pair.effective.hamiltonian_spec isa EffectiveHamiltonianSpec
    @test pair.circuit.hamiltonian_spec isa CircuitHamiltonianSpec
    @test snapshot.targets["cr_f01_ghz"] == 4.3
    @test subsystems(qr_stage.system)[2].ω == 4.3
    @test subsystems(qcr_stage.system)[3].ω == 4.3
    @test snapshot.devices[:QB1]["fit_classification"] == "dressed_anchor_backout"
    @test snapshot.devices[:QB2]["fit_classification"] == "dressed_anchor_backout"
    @test snapshot.devices[:TC1]["flux"] == 0.0
    @test snapshot.devices[:TC2]["flux"] == 0.0
    @test snapshot.devices[:TC1]["asymmetry"] == 0.1
    @test snapshot.devices[:TC2]["asymmetry"] == 0.1
    @test isapprox(snapshot.devices[:QB1]["f01_ghz"], 4.6679; atol = 0.003)
    @test isapprox(snapshot.devices[:QB2]["f01_ghz"], 4.47033; atol = 0.003)
    @test isapprox(snapshot.devices[:TC1]["f01_ghz"], 6.5; atol = 0.02)
    @test isapprox(snapshot.devices[:TC2]["f01_ghz"], 6.5; atol = 0.02)
    @test snapshot.devices[:TC1]["EJmax"] > 50.0
    @test snapshot.devices[:TC2]["EJmax"] > 50.0
    @test isapprox(
        qr_stage.system.couplings[1].g,
        snapshot.targets["beta_qr"] * sqrt(snapshot.devices[:QB1]["f01_ghz"] * snapshot.targets["cr_f01_ghz"]);
        atol = 1e-4,
    )
    @test isapprox(
        qcr_stage.system.couplings[1].g,
        snapshot.targets["beta_qc_qb1"] * sqrt(snapshot.devices[:QB1]["f01_ghz"] * snapshot.devices[:TC1]["f01_ghz"]);
        atol = 1e-4,
    )
    @test isapprox(
        pair.effective.system.couplings[4].g,
        snapshot.targets["beta_qc_qb2"] * sqrt(snapshot.devices[:QB2]["f01_ghz"] * snapshot.devices[:TC2]["f01_ghz"]);
        atol = 1e-4,
    )
    @test isapprox(
        qcr_stage.system.couplings[2].g,
        snapshot.targets["beta_cr"] * sqrt(snapshot.devices[:TC1]["f01_ghz"] * snapshot.targets["cr_f01_ghz"]);
        atol = 1e-4,
    )
    @test isapprox(
        pair.effective.system.couplings[3].g,
        snapshot.targets["beta_cr"] * sqrt(snapshot.devices[:TC2]["f01_ghz"] * snapshot.targets["cr_f01_ghz"]);
        atol = 1e-4,
    )

    effective_model = build_model(pair.effective.system)
    circuit_model = build_model(pair.circuit.system; hamiltonian_spec = pair.circuit.hamiltonian_spec)

    @test size(hamiltonian(effective_model).data) == (1200, 1200)
    @test size(hamiltonian(circuit_model).data) == (1875, 1875)
end

@testset "Renger 2026 explicit effective coupling seeds" begin
    snapshot = load_renger2026_snapshot()
    shared_snapshot = deepcopy(snapshot)
    shared_snapshot.devices[:TC1]["f01_ghz"] = 9.9
    shared_snapshot.devices[:TC2]["f01_ghz"] = 6.1
    shared_snapshot.targets["beta_cr"] = 0.031
    shared_snapshot.targets["g_qr_qb1"] = 0.0091
    shared_snapshot.targets["g_qr_qb2"] = 0.0087
    shared_snapshot.targets["g_qc_qb1_tc1"] = 0.082
    shared_snapshot.targets["g_qc_qb2_tc2"] = 0.107
    shared_snapshot.targets["g_cr_tc1"] = 0.123
    shared_snapshot.targets["g_cr_tc2"] = 0.124

    qr_stage = renger2026_stage1_qr_system(shared_snapshot)
    qcr_stage = renger2026_stage1_qcr_system(shared_snapshot; qubit = :QB1, coupler = :TC1)
    effective = renger2026_reduced_system(shared_snapshot; model = :effective)
    circuit = renger2026_reduced_system(shared_snapshot; model = :circuit)

    @test qr_stage.system.couplings[1].g == 0.0091
    @test qcr_stage.system.couplings[1].g == 0.082
    @test qcr_stage.system.couplings[2].g == 0.123
    @test effective.system.couplings[1].g == 0.082
    @test effective.system.couplings[2].g == 0.123
    @test effective.system.couplings[3].g == 0.124
    @test effective.system.couplings[4].g == 0.107
    @test circuit.system.couplings[2].G == 0.031
    @test circuit.system.couplings[3].G == 0.031
end

@testset "Renger 2026 snapshot overlay merge" begin
    baseline = load_renger2026_snapshot()
    mktemp() do path, io
        write(
            io,
            """
            [devices.TC1]
            EJmax = 38.5
            asymmetry = 0.14
            f01_ghz = 5.9

            [targets]
            beta_qc_qb1 = 0.019
            beta_cr = 0.028
            fig2_move_t_gate_ns = 90.0
            """,
        )
        close(io)

        snapshot = load_renger2026_snapshot(; overlay_path = path)

        @test snapshot.devices[:TC1]["EJmax"] == 38.5
        @test snapshot.devices[:TC1]["asymmetry"] == 0.14
        @test snapshot.devices[:TC1]["f01_ghz"] == 5.9
        @test snapshot.devices[:TC1]["EC"] == baseline.devices[:TC1]["EC"]
        @test snapshot.devices[:TC2]["EJmax"] == baseline.devices[:TC2]["EJmax"]
        @test snapshot.targets["beta_qc_qb1"] == 0.019
        @test snapshot.targets["beta_qc_qb2"] == 0.022
        @test snapshot.targets["beta_cr"] == 0.028
        @test snapshot.targets["fig2_move_t_gate_ns"] == 90.0
        @test occursin(path, snapshot.path)
    end
end

@testset "Renger 2026 working overlay schema" begin
    overlay = TOML.parsefile(joinpath(dirname(@__DIR__), "output", "renger2026", "fig2_ef_retune_working.toml"))

    @test overlay["source_snapshot"] == "output/renger2026/paper_local_priors.toml"
    @test Set(keys(overlay["devices"])) == Set(["TC1", "TC2"])

    allowed_device_keys = Set(["EJmax", "EC", "asymmetry", "f01_ghz", "anharmonicity_ghz"])
    for device_name in ("TC1", "TC2")
        @test Set(keys(overlay["devices"][device_name])) ⊆ allowed_device_keys
        @test !haskey(overlay["devices"][device_name], "flux")
    end

    allowed_target_keys = Set(["beta_qc_qb1", "beta_qc_qb2", "beta_cr", "fig2_move_t_gate_ns", "fig2_cz_t_gate_ns"])
    @test Set(keys(overlay["targets"])) ⊆ allowed_target_keys
end

@testset "MOVE workflow helpers" begin
    qr_sys = CompositeSystem(
        Transmon(:q; EJ = 0.9, EC = 0.2, ncut = 3),
        Resonator(:r; ω = 1.0, dim = 2),
        CapacitiveCoupling(:q, :r; g = 0.08),
    )
    qr_spec = subspace_spec(qr_sys; subsystem_levels = (q = [0, 1], r = [0, 1]))
    qr_trace = projected_unitary(qr_sys, qr_spec, collect(range(0.0, 40.0; length = 401)))
    qr_move = best_move(qr_trace; source_index = 3, target_index = 2)

    @test qr_move.transfer_probability > 0.99
    @test qr_move.leakage < 1e-6

    qcr_sys = CompositeSystem(
        Transmon(:q; EJ = 0.9, EC = 0.2, ncut = 3),
        TunableCoupler(:c; EJmax = 2.2, EC = 0.15, flux = 0.0, asymmetry = 0.1, ncut = 3),
        Resonator(:r; ω = 1.0, dim = 2),
        CapacitiveCoupling(:q, :c; g = 0.07),
        CapacitiveCoupling(:c, :r; g = 0.07),
    )
    qcr_tlist = collect(range(0.0, 80.0; length = 401))
    qcr_result = evolve(qcr_sys, basis_state(qcr_sys; q = 1, c = 0, r = 0), qcr_tlist)
    qcr_trace = projected_unitary(qcr_sys, subspace_spec(qcr_sys; subsystem_levels = (q = [0, 1], r = [0, 1])), qcr_tlist)
    qcr_move = best_move(qcr_trace; source_index = 3, target_index = 2)

    @test maximum(population_trace(qcr_result, :c, 1).values) > 0.05
    @test maximum(population_trace(qcr_result, :r, 1).values) > 0.45
    @test qcr_move.transfer_probability > 0.45
    @test qcr_move.leakage < 0.1

    snapshot = load_renger2026_snapshot()
    full_pair = renger2026_model_pair(snapshot)
    full_sys = full_pair.effective.system
    full_tlist = collect(range(0.0, 5.0; length = 51))
    full_result = evolve(
        full_sys,
        basis_state(full_sys; QB1 = 1, TC1 = 0, CR = 0, TC2 = 0, QB2 = 0),
        full_tlist;
        observables = [ObservableSpec(:ncr, :CR, :n)],
    )
    detuned_sys = with_subsystem_parameter(full_sys, :CR, :ω, 6.0)
    detuned_result = evolve(
        detuned_sys,
        basis_state(detuned_sys; QB1 = 1, TC1 = 0, CR = 0, TC2 = 0, QB2 = 0),
        full_tlist;
        observables = [ObservableSpec(:ncr, :CR, :n)],
    )

    full_ncr_peak = maximum(real.(observable_trace(full_result, :ncr).values))
    detuned_ncr_peak = maximum(real.(observable_trace(detuned_result, :ncr).values))
    full_qb2_peak = maximum(population_trace(full_result, :QB2, 1).values)
    detuned_qb2_peak = maximum(population_trace(detuned_result, :QB2, 1).values)

    @test abs(full_ncr_peak - detuned_ncr_peak) > 1e-5
    @test abs(full_qb2_peak - detuned_qb2_peak) > 1e-7
end

@testset "CZ workflow helpers" begin
    q1 = Transmon(:q1; EJ = 0.9, EC = 0.2, ng = 0.0, ncut = 3)
    q2 = Transmon(:q2; EJ = 0.95, EC = 0.2, ng = 0.0, ncut = 3)
    c = TunableCoupler(:c; EJmax = 2.6, EC = 0.12, flux = 0.0, asymmetry = 0.1, ng = 0.0, ncut = 3)
    cz_sys = CompositeSystem(
        q1,
        c,
        q2,
        CircuitCapacitiveCoupling(:q1, :c; G = 0.03),
        CircuitCapacitiveCoupling(:c, :q2; G = 0.03),
    )
    circuit_spec = CircuitHamiltonianSpec(charge_cutoff = 2)
    cz_subspace = subspace_spec(
        cz_sys;
        hamiltonian_spec = circuit_spec,
        subsystem_levels = (q1 = [0, 1], q2 = [0, 1]),
        basis = :dressed_static,
    )
    flux = FluxControl(:cz_flux, :c, (p, t) -> p.δ)
    cz_trace = projected_unitary(
        cz_sys,
        cz_subspace,
        collect(range(0.0, 100.0; length = 201));
        hamiltonian_spec = circuit_spec,
        flux_controls = [flux],
        params = (; δ = 0.04),
    )
    cz_point = best_cz(cz_trace)

    @test wrapped_phase_distance(cz_point.conditional_phase, π) < 0.01
    @test cz_point.leakage < 0.01
end
