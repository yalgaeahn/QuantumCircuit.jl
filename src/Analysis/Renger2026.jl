using TOML
using ..Architecture:
    CapacitiveCoupling,
    CircuitCapacitiveCoupling,
    CompositeSystem,
    Resonator,
    TunableCoupler,
    TunableTransmon
using ..Model: CircuitHamiltonianSpec, EffectiveHamiltonianSpec

struct Renger2026Snapshot
    path::String
    devices::Dict{Symbol, Dict{String, Any}}
    targets::Dict{String, Any}
end

function load_renger2026_snapshot(path::AbstractString = _renger2026_snapshot_path())
    parsed = TOML.parsefile(path)
    devices = Dict(Symbol(name) => Dict{String, Any}(data) for (name, data) in parsed["devices"])
    targets = Dict{String, Any}(parsed["targets"])
    return Renger2026Snapshot(String(path), devices, targets)
end

function renger2026_stage1_qr_system(
    snapshot::Renger2026Snapshot = load_renger2026_snapshot();
    qubit::Symbol = :QB1,
    qubit_ncut::Integer = 5,
    resonator_dim::Integer = 3,
)
    targets = snapshot.targets
    qubit_device = _renger2026_transmon(snapshot.devices, qubit; ncut = qubit_ncut)
    resonator = Resonator(:CR; ω = Float64(targets["cr_f01_ghz"]), dim = resonator_dim)
    coupling = CapacitiveCoupling(qubit, :CR; g = _seed_effective_coupling(Float64(targets["beta_qr"]), _device_frequency(snapshot.devices[qubit]), resonator.ω))
    return (; system = CompositeSystem(qubit_device, resonator, coupling), hamiltonian_spec = EffectiveHamiltonianSpec())
end

function renger2026_stage1_qcr_system(
    snapshot::Renger2026Snapshot = load_renger2026_snapshot();
    qubit::Symbol = :QB1,
    coupler::Symbol = :TC1,
    qubit_ncut::Integer = 5,
    coupler_ncut::Integer = 4,
    resonator_dim::Integer = 3,
)
    targets = snapshot.targets
    qubit_device = _renger2026_transmon(snapshot.devices, qubit; ncut = qubit_ncut)
    coupler_device = _renger2026_coupler(snapshot.devices, coupler; ncut = coupler_ncut)
    resonator = Resonator(:CR; ω = Float64(targets["cr_f01_ghz"]), dim = resonator_dim)
    g_qc = _seed_effective_coupling(Float64(_beta_qc(targets, qubit)), _device_frequency(snapshot.devices[qubit]), _device_frequency(snapshot.devices[coupler]))
    g_cr = _seed_effective_coupling(Float64(targets["beta_cr"]), _device_frequency(snapshot.devices[coupler]), resonator.ω)

    system = CompositeSystem(
        qubit_device,
        coupler_device,
        resonator,
        CapacitiveCoupling(qubit, coupler; g = g_qc),
        CapacitiveCoupling(coupler, :CR; g = g_cr),
    )
    return (; system, hamiltonian_spec = EffectiveHamiltonianSpec())
end

function renger2026_reduced_system(
    snapshot::Renger2026Snapshot = load_renger2026_snapshot();
    model::Symbol = :effective,
    qubit_ncut::Integer = 5,
    coupler_ncut::Integer = 4,
    resonator_dim::Integer = 3,
    charge_cutoff::Integer = 2,
)
    q1 = _renger2026_transmon(snapshot.devices, :QB1; ncut = qubit_ncut)
    q2 = _renger2026_transmon(snapshot.devices, :QB2; ncut = qubit_ncut)
    c1 = _renger2026_coupler(snapshot.devices, :TC1; ncut = coupler_ncut)
    c2 = _renger2026_coupler(snapshot.devices, :TC2; ncut = coupler_ncut)
    resonator = Resonator(:CR; ω = Float64(snapshot.targets["cr_f01_ghz"]), dim = resonator_dim)

    if model == :effective
        g_q1c1 = _seed_effective_coupling(Float64(_beta_qc(snapshot.targets, :QB1)), _device_frequency(snapshot.devices[:QB1]), _device_frequency(snapshot.devices[:TC1]))
        g_q2c2 = _seed_effective_coupling(Float64(_beta_qc(snapshot.targets, :QB2)), _device_frequency(snapshot.devices[:QB2]), _device_frequency(snapshot.devices[:TC2]))
        g_c1r = _seed_effective_coupling(Float64(snapshot.targets["beta_cr"]), _device_frequency(snapshot.devices[:TC1]), resonator.ω)
        g_c2r = _seed_effective_coupling(Float64(snapshot.targets["beta_cr"]), _device_frequency(snapshot.devices[:TC2]), resonator.ω)

        system = CompositeSystem(
            q1,
            c1,
            resonator,
            c2,
            q2,
            CapacitiveCoupling(:QB1, :TC1; g = g_q1c1),
            CapacitiveCoupling(:TC1, :CR; g = g_c1r),
            CapacitiveCoupling(:CR, :TC2; g = g_c2r),
            CapacitiveCoupling(:TC2, :QB2; g = g_q2c2),
        )
        return (; system, hamiltonian_spec = EffectiveHamiltonianSpec())
    elseif model == :circuit
        system = CompositeSystem(
            q1,
            c1,
            resonator,
            c2,
            q2,
            CircuitCapacitiveCoupling(:QB1, :TC1; G = Float64(_beta_qc(snapshot.targets, :QB1))),
            CircuitCapacitiveCoupling(:TC1, :CR; G = Float64(snapshot.targets["beta_cr"])),
            CircuitCapacitiveCoupling(:CR, :TC2; G = Float64(snapshot.targets["beta_cr"])),
            CircuitCapacitiveCoupling(:TC2, :QB2; G = Float64(_beta_qc(snapshot.targets, :QB2))),
        )
        return (; system, hamiltonian_spec = CircuitHamiltonianSpec(charge_cutoff = charge_cutoff))
    end

    throw(ArgumentError("model must be :effective or :circuit, got $model."))
end

function renger2026_model_pair(
    snapshot::Renger2026Snapshot = load_renger2026_snapshot();
    qubit_ncut::Integer = 5,
    coupler_ncut::Integer = 4,
    resonator_dim::Integer = 3,
    charge_cutoff::Integer = 2,
)
    return (
        effective = renger2026_reduced_system(
            snapshot;
            model = :effective,
            qubit_ncut = qubit_ncut,
            coupler_ncut = coupler_ncut,
            resonator_dim = resonator_dim,
            charge_cutoff = charge_cutoff,
        ),
        circuit = renger2026_reduced_system(
            snapshot;
            model = :circuit,
            qubit_ncut = qubit_ncut,
            coupler_ncut = coupler_ncut,
            resonator_dim = resonator_dim,
            charge_cutoff = charge_cutoff,
        ),
    )
end

_renger2026_snapshot_path() = joinpath(dirname(dirname(@__DIR__)), "output", "renger2026", "paper_local_priors.toml")

function _renger2026_transmon(devices::Dict{Symbol, Dict{String, Any}}, name::Symbol; ncut::Integer)
    device = devices[name]
    return TunableTransmon(
        name;
        EJmax = Float64(device["EJmax"]),
        EC = Float64(device["EC"]),
        flux = Float64(device["flux"]),
        asymmetry = Float64(device["asymmetry"]),
        ng = Float64(device["ng"]),
        ncut = Int(ncut),
    )
end

function _renger2026_coupler(devices::Dict{Symbol, Dict{String, Any}}, name::Symbol; ncut::Integer)
    device = devices[name]
    return TunableCoupler(
        name;
        EJmax = Float64(device["EJmax"]),
        EC = Float64(device["EC"]),
        flux = Float64(device["flux"]),
        asymmetry = Float64(device["asymmetry"]),
        ng = Float64(device["ng"]),
        ncut = Int(ncut),
    )
end

_device_frequency(device::Dict{String, Any}) = Float64(device["f01_ghz"])
_seed_effective_coupling(beta::Float64, ωa::Float64, ωb::Float64) = beta * sqrt(ωa * ωb)
_beta_qc(targets::Dict{String, Any}, qubit::Symbol) = qubit == :QB1 ? targets["beta_qc_qb1"] : targets["beta_qc_qb2"]
