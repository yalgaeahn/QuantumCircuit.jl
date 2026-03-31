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

function load_renger2026_snapshot(
    path::AbstractString = _renger2026_snapshot_path();
    overlay_path::Union{Nothing, AbstractString} = nothing,
)
    parsed = TOML.parsefile(path)
    if overlay_path !== nothing
        overlay = TOML.parsefile(overlay_path)
        parsed = _deep_merge_dict(parsed, overlay)
    end
    devices = Dict(Symbol(name) => Dict{String, Any}(data) for (name, data) in parsed["devices"])
    targets = Dict{String, Any}(parsed["targets"])
    merged_path = overlay_path === nothing ? String(path) : string(path, " + ", overlay_path)
    return Renger2026Snapshot(merged_path, devices, targets)
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
    coupling = CapacitiveCoupling(qubit, :CR; g = _g_qr(targets, qubit))
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
    g_qc = _g_qc(targets, qubit, coupler)
    g_cr = _g_cr(targets, coupler)

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
    # WARNING: charge_cutoff=2 is the default for structural/dimensional testing only.
    # For physically accurate circuit-model spectra, charge_cutoff must be large enough
    # to converge the local energy levels (verified at flux=0, the hardest point):
    #   QB1/QB2 (EJmax~17 GHz):            cc=5  → converged (<1 MHz from cc=10)
    #   TC retune (EJmax~32–36 GHz):       cc=10 → converged (<0.1 MHz from cc=13)
    #   TC baseline (EJmax~51 GHz):        cc=10 → converged (<0.1 MHz from cc=13)
    # At cc=2: QB1 appears at ~7.2 GHz (target 4.67 GHz), TC1 at ~19 GHz (target 6.5 GHz).
    # At cc=7: TC1 retune still 16 MHz too high — not converged.
    # For the 3-body QCR branch models used in Fig. 2, use charge_cutoff=10
    # (Hilbert space 21²×3 = 1323 states — fast and fully converged).
)
    q1 = _renger2026_transmon(snapshot.devices, :QB1; ncut = qubit_ncut)
    q2 = _renger2026_transmon(snapshot.devices, :QB2; ncut = qubit_ncut)
    c1 = _renger2026_coupler(snapshot.devices, :TC1; ncut = coupler_ncut)
    c2 = _renger2026_coupler(snapshot.devices, :TC2; ncut = coupler_ncut)
    resonator = Resonator(:CR; ω = Float64(snapshot.targets["cr_f01_ghz"]), dim = resonator_dim)

    if model == :effective
        g_q1c1 = _g_qc(snapshot.targets, :QB1, :TC1)
        g_q2c2 = _g_qc(snapshot.targets, :QB2, :TC2)
        g_c1r = _g_cr(snapshot.targets, :TC1)
        g_c2r = _g_cr(snapshot.targets, :TC2)

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
            CircuitCapacitiveCoupling(:TC1, :CR; G = Float64(_beta_cr(snapshot.targets))),
            CircuitCapacitiveCoupling(:CR, :TC2; G = Float64(_beta_cr(snapshot.targets))),
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

_beta_qc(targets::Dict{String, Any}, qubit::Symbol) = qubit == :QB1 ? targets["beta_qc_qb1"] : targets["beta_qc_qb2"]
_beta_cr(targets::Dict{String, Any}) = targets["beta_cr"]
_g_qr(targets::Dict{String, Any}, qubit::Symbol) = qubit == :QB1 ? Float64(targets["g_qr_qb1"]) : qubit == :QB2 ? Float64(targets["g_qr_qb2"]) : throw(ArgumentError("Unsupported qubit $qubit for QR seed coupling."))
_g_qc(targets::Dict{String, Any}, qubit::Symbol, coupler::Symbol) = qubit == :QB1 && coupler == :TC1 ? Float64(targets["g_qc_qb1_tc1"]) : qubit == :QB2 && coupler == :TC2 ? Float64(targets["g_qc_qb2_tc2"]) : throw(ArgumentError("Unsupported qubit-coupler pair ($qubit, $coupler) for QC seed coupling."))
_g_cr(targets::Dict{String, Any}, coupler::Symbol) = coupler == :TC1 ? Float64(targets["g_cr_tc1"]) : coupler == :TC2 ? Float64(targets["g_cr_tc2"]) : throw(ArgumentError("Unsupported coupler $coupler for CR seed coupling."))

function _deep_merge_dict(base::Dict, overlay::Dict)
    merged = deepcopy(base)
    for (key, value) in overlay
        if haskey(merged, key) && merged[key] isa Dict && value isa Dict
            merged[key] = _deep_merge_dict(merged[key], value)
        else
            merged[key] = deepcopy(value)
        end
    end
    return merged
end
