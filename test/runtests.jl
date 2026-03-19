using Test
using QuantumCircuit
using QuantumToolbox: expect

effective_EJ(EJmax, flux, asymmetry) = EJmax * sqrt(cospi(flux)^2 + asymmetry^2 * sinpi(flux)^2)

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
    @test_throws ArgumentError TunableTransmon(:tq1; EJmax = 20.0, EC = 0.25, asymmetry = 1.2)
    @test_throws ArgumentError TunableCoupler(:c1; EJmax = -15.0, EC = 0.30, flux = 0.0)

    @test TunableTransmon(:tq1; EJmax = 20.0, EC = 0.25, ng = 0.2).ng == 0.2
    @test TunableCoupler(:c1; EJmax = 15.0, EC = 0.30).ng == 0.0
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

@testset "Circuit Hamiltonian and operators" begin
    circuit_spec = CircuitHamiltonianSpec(charge_cutoff = 2)
    q_ng0 = Transmon(:q; EJ = 20.0, EC = 0.25, ng = 0.0, ncut = 6)
    q_ng1 = Transmon(:q; EJ = 20.0, EC = 0.25, ng = 0.25, ncut = 6)

    effective_spec0 = spectrum(CompositeSystem(q_ng0); levels = 4)
    effective_spec1 = spectrum(CompositeSystem(q_ng1); levels = 4)
    circuit_spec0 = spectrum(CompositeSystem(q_ng0); levels = 4, hamiltonian_spec = circuit_spec)
    circuit_spec1 = spectrum(CompositeSystem(q_ng1); levels = 4, hamiltonian_spec = circuit_spec)

    @test effective_spec0.energies ≈ effective_spec1.energies atol = 1e-10
    @test circuit_spec0.energies != circuit_spec1.energies

    tq_flux0 = TunableTransmon(:tq; EJmax = 20.0, EC = 0.25, flux = 0.0, asymmetry = 0.0, ng = 0.0, ncut = 6)
    tq_flux1 = TunableTransmon(:tq; EJmax = 20.0, EC = 0.25, flux = 0.25, asymmetry = 0.0, ng = 0.0, ncut = 6)
    @test transition_frequencies(spectrum(CompositeSystem(tq_flux1); levels = 4, hamiltonian_spec = circuit_spec))[1] <
        transition_frequencies(spectrum(CompositeSystem(tq_flux0); levels = 4, hamiltonian_spec = circuit_spec))[1]

    resonator = Resonator(:r1; ω = 1.25, dim = 4)
    @test spectrum(CompositeSystem(resonator); levels = 4, hamiltonian_spec = circuit_spec).energies ≈
        spectrum(CompositeSystem(resonator); levels = 4).energies atol = 1e-10

    circuit_model = build_model(CompositeSystem(q_ng0); hamiltonian_spec = circuit_spec)
    effective_model = build_model(CompositeSystem(q_ng0))
    resonator_model = build_model(CompositeSystem(resonator); hamiltonian_spec = circuit_spec)

    @test size(charge_operator(circuit_model, :q).data) == (5, 5)
    @test size(cosphi_operator(circuit_model, :q).data) == (5, 5)
    @test size(sinphi_operator(circuit_model, :q).data) == (5, 5)
    @test_throws ArgumentError annihilation_operator(circuit_model, :q)
    @test_throws ArgumentError number_operator(circuit_model, :q)
    @test_throws ArgumentError charge_operator(effective_model, :q)
    @test size(number_operator(resonator_model, :r1).data) == (4, 4)
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
