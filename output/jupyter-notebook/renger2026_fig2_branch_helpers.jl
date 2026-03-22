round4(x::Real) = round(Float64(x); digits = 4)
compact_value(x) = x isa Real ? round4(x) : x
compact_namedtuple(nt::NamedTuple) = (; (key => compact_value(value) for (key, value) in pairs(nt))...)

const FIG2_COLORS = (
    diabatic_a = :gray40,
    diabatic_b = :gray65,
    dressed_a = :royalblue3,
    dressed_b = :firebrick3,
)

const FIG2_BRANCHES = (
    move = (label = "QB1-TC1-CR", qb_key = :QB1, tc_key = :TC1, beta_qc_key = "beta_qc_qb1"),
    cz = (label = "QB2-TC2-CR", qb_key = :QB2, tc_key = :TC2, beta_qc_key = "beta_qc_qb2"),
)

branch_spec(branch::Symbol) = getproperty(FIG2_BRANCHES, branch)

function branch_flux_ranges(branch::Symbol)
    branch == :move && return ((0.18, 0.26), (0.12, 0.40))
    return ((0.12, 0.24), (0.12, 0.40))
end

function branch_flux_axes(branch::Symbol; q_points, tc_points)
    (q_lo, q_hi), (tc_lo, tc_hi) = branch_flux_ranges(branch)
    return (
        q_fluxes = collect(range(q_lo, q_hi; length = q_points)),
        tc_fluxes = collect(range(tc_lo, tc_hi; length = tc_points)),
    )
end

function fig2_gate_times(snapshot; move_default = 100.0, cz_default = 120.0)
    return (
        move = Float64(get(snapshot.targets, "fig2_move_t_gate_ns", move_default)),
        cz = Float64(get(snapshot.targets, "fig2_cz_t_gate_ns", cz_default)),
    )
end

function parameter_rows(snapshot; move_t_gate_ns, cz_t_gate_ns, charge_cutoff)
    return [
        (
            device = "QB1",
            fit = snapshot.devices[:QB1]["fit_classification"],
            f01_ghz = round4(snapshot.devices[:QB1]["f01_ghz"]),
            EC_ghz = round4(snapshot.devices[:QB1]["EC"]),
            EJmax_ghz = round4(snapshot.devices[:QB1]["EJmax"]),
            flux = round4(snapshot.devices[:QB1]["flux"]),
            asymmetry = round4(snapshot.devices[:QB1]["asymmetry"]),
        ),
        (
            device = "TC1",
            fit = snapshot.devices[:TC1]["fit_classification"],
            f01_ghz = round4(snapshot.devices[:TC1]["f01_ghz"]),
            EC_ghz = round4(snapshot.devices[:TC1]["EC"]),
            EJmax_ghz = round4(snapshot.devices[:TC1]["EJmax"]),
            flux = round4(snapshot.devices[:TC1]["flux"]),
            asymmetry = round4(snapshot.devices[:TC1]["asymmetry"]),
        ),
        (
            device = "QB2",
            fit = snapshot.devices[:QB2]["fit_classification"],
            f01_ghz = round4(snapshot.devices[:QB2]["f01_ghz"]),
            EC_ghz = round4(snapshot.devices[:QB2]["EC"]),
            EJmax_ghz = round4(snapshot.devices[:QB2]["EJmax"]),
            flux = round4(snapshot.devices[:QB2]["flux"]),
            asymmetry = round4(snapshot.devices[:QB2]["asymmetry"]),
        ),
        (
            device = "TC2",
            fit = snapshot.devices[:TC2]["fit_classification"],
            f01_ghz = round4(snapshot.devices[:TC2]["f01_ghz"]),
            EC_ghz = round4(snapshot.devices[:TC2]["EC"]),
            EJmax_ghz = round4(snapshot.devices[:TC2]["EJmax"]),
            flux = round4(snapshot.devices[:TC2]["flux"]),
            asymmetry = round4(snapshot.devices[:TC2]["asymmetry"]),
        ),
        (
            device = "CR",
            fit = "paper target",
            f01_ghz = round4(snapshot.targets["cr_f01_ghz"]),
            EC_ghz = missing,
            EJmax_ghz = missing,
            flux = missing,
            asymmetry = missing,
        ),
    ], (
        beta_qc_qb1 = round4(snapshot.targets["beta_qc_qb1"]),
        beta_qc_qb2 = round4(snapshot.targets["beta_qc_qb2"]),
        beta_cr = round4(snapshot.targets["beta_cr"]),
        cr_f01_ghz = round4(snapshot.targets["cr_f01_ghz"]),
        move_t_gate_ns = round4(move_t_gate_ns),
        cz_t_gate_ns = round4(cz_t_gate_ns),
        charge_cutoff = charge_cutoff,
    )
end

function subsystem_by_name(system::CompositeSystem, target::Symbol)
    for subsystem in subsystems(system)
        name(subsystem) == target && return subsystem
    end
    error("No subsystem named $(target) exists in the system.")
end

function snapshot_with_overrides(snapshot; device_overrides = Dict{Symbol, Any}(), target_overrides = Dict{String, Any}())
    devices = Dict(name => deepcopy(data) for (name, data) in snapshot.devices)
    for (device_name, overrides) in device_overrides
        device = get!(devices, device_name, Dict{String, Any}())
        for (key, value) in pairs(overrides)
            device[String(key)] = value
        end
    end

    targets = Dict{String, Any}(snapshot.targets)
    for (key, value) in pairs(target_overrides)
        targets[String(key)] = value
    end

    return typeof(snapshot)(snapshot.path, devices, targets)
end

branch_tc_flux_offset(snapshot, branch::Symbol) = Float64(get(snapshot.devices[branch_spec(branch).tc_key], "flux", 0.0))
effective_tc_flux(snapshot, tc_flux; branch = :move) = tc_flux + branch_tc_flux_offset(snapshot, branch)

function local_qb(snapshot, flux; branch = :move, ncut = 5)
    spec = branch_spec(branch)
    p = snapshot.devices[spec.qb_key]
    return TunableTransmon(
        :QB;
        EJmax = p["EJmax"],
        EC = p["EC"],
        flux = flux,
        asymmetry = p["asymmetry"],
        ng = p["ng"],
        ncut = ncut,
    )
end

function local_tc(snapshot, flux; branch = :move, ncut = 4)
    spec = branch_spec(branch)
    p = snapshot.devices[spec.tc_key]
    return TunableCoupler(
        :TC;
        EJmax = p["EJmax"],
        EC = p["EC"],
        flux = effective_tc_flux(snapshot, flux; branch),
        asymmetry = p["asymmetry"],
        ng = p["ng"],
        ncut = ncut,
    )
end

local_cr(snapshot; dim = 3) = Resonator(:CR; ω = snapshot.targets["cr_f01_ghz"], dim = dim)
local_spec(; charge_cutoff) = CircuitHamiltonianSpec(charge_cutoff = charge_cutoff)

function local_exact_system(snapshot; branch = :move, q_flux = 0.20, tc_flux = 0.34, qb_ncut = 5, tc_ncut = 4, cr_dim = 3)
    spec = branch_spec(branch)
    return CompositeSystem(
        local_qb(snapshot, q_flux; branch, ncut = qb_ncut),
        local_tc(snapshot, tc_flux; branch, ncut = tc_ncut),
        local_cr(snapshot; dim = cr_dim),
        CircuitCapacitiveCoupling(:QB, :TC; G = snapshot.targets[spec.beta_qc_key]),
        CircuitCapacitiveCoupling(:TC, :CR; G = snapshot.targets["beta_cr"]),
    )
end

function retuned_system(base_system::CompositeSystem, snapshot; branch = :move, q_flux, tc_flux)
    step1 = with_subsystem_parameter(base_system, :QB, :flux, q_flux)
    return with_subsystem_parameter(step1, :TC, :flux, effective_tc_flux(snapshot, tc_flux; branch))
end

function local_eigenbundle(subsystem; charge_cutoff)
    model = build_model(CompositeSystem(subsystem); hamiltonian_spec = CircuitHamiltonianSpec(charge_cutoff = charge_cutoff))
    es = eigenstates(hamiltonian(model))
    values = Float64.(real.(es.values))
    dims = basis(model.dimensions[name(subsystem)], 0).dims
    states = [QuantumObject(es.vectors[:, i]; type = Ket(), dims = dims) for i in axes(es.vectors, 2)]
    return values, states
end

q_bundle(snapshot, flux; branch = :move, charge_cutoff, qb_ncut = 5) = local_eigenbundle(local_qb(snapshot, flux; branch, ncut = qb_ncut); charge_cutoff)
tc_bundle(snapshot, flux; branch = :move, charge_cutoff, tc_ncut = 4) = local_eigenbundle(local_tc(snapshot, flux; branch, ncut = tc_ncut); charge_cutoff)
cr_bundle(snapshot; resonator_dim = 3) = ([snapshot.targets["cr_f01_ghz"] * n for n in 0:(resonator_dim - 1)], [basis(resonator_dim, n) for n in 0:(resonator_dim - 1)])

function qubit_frequency(snapshot, flux; branch = :move, charge_cutoff, qb_ncut = 5)
    values, _ = q_bundle(snapshot, flux; branch, charge_cutoff, qb_ncut)
    return values[2] - values[1]
end

function coupler_frequency(snapshot, flux; branch = :move, charge_cutoff, tc_ncut = 4)
    values, _ = tc_bundle(snapshot, flux; branch, charge_cutoff, tc_ncut)
    return values[2] - values[1]
end

function diabatic_state(snapshot, q_flux, tc_flux, qlevel, tclevel, crlevel; branch = :move, charge_cutoff, qb_ncut = 5, tc_ncut = 4, resonator_dim = 3)
    _, qstates = q_bundle(snapshot, q_flux; branch, charge_cutoff, qb_ncut)
    _, tcstates = tc_bundle(snapshot, tc_flux; branch, charge_cutoff, tc_ncut)
    _, crstates = cr_bundle(snapshot; resonator_dim)
    return kron(qstates[qlevel + 1], tcstates[tclevel + 1], crstates[crlevel + 1])
end

function adiabatic_fan_rows(base_system, snapshot, q_fluxes, tc_flux; branch_pair_a, branch_pair_b, which, branch, charge_cutoff, fan_levels = 8, qb_ncut = 5, tc_ncut = 4, resonator_dim = 3)
    spec = local_spec(; charge_cutoff)
    rows = NamedTuple[]

    for q_flux in q_fluxes
        system = retuned_system(base_system, snapshot; branch, q_flux, tc_flux)
        model = build_model(system; hamiltonian_spec = spec)
        es = eigenstates(hamiltonian(model))
        raw_energies = Float64.(real.(es.values))
        order = sortperm(raw_energies)
        energies = raw_energies[order]
        vectors = es.vectors[:, order]
        ground = energies[1]
        dims = diabatic_state(snapshot, q_flux, tc_flux, 0, 0, 0; branch, charge_cutoff, qb_ncut, tc_ncut, resonator_dim).dims
        coupled_states = [QuantumObject(vectors[:, i]; type = Ket(), dims = dims) for i in axes(vectors, 2)]
        q_values, _ = q_bundle(snapshot, q_flux; branch, charge_cutoff, qb_ncut)
        tc_values, _ = tc_bundle(snapshot, tc_flux; branch, charge_cutoff, tc_ncut)

        target_a = diabatic_state(snapshot, q_flux, tc_flux, branch_pair_a...; branch, charge_cutoff, qb_ncut, tc_ncut, resonator_dim)
        target_b = diabatic_state(snapshot, q_flux, tc_flux, branch_pair_b...; branch, charge_cutoff, qb_ncut, tc_ncut, resonator_dim)
        overlaps_a = Float64[abs2(dag(target_a) * ψ) for ψ in coupled_states]
        overlaps_b = Float64[abs2(dag(target_b) * ψ) for ψ in coupled_states]

        last_level = min(length(energies), fan_levels + 1)
        fan_abs_indices = collect(2:last_level)
        fan_ghz = [energies[idx] - ground for idx in fan_abs_indices]
        fan_overlaps_a = [overlaps_a[idx] for idx in fan_abs_indices]
        fan_overlaps_b = [overlaps_b[idx] for idx in fan_abs_indices]
        fan_total_support = fan_overlaps_a .+ fan_overlaps_b
        issorted(fan_ghz) || error("Adiabatic fan is not sorted at q_flux=$(q_flux), tc_flux=$(tc_flux).")

        push!(rows, (
            q_flux = q_flux,
            tc_flux = tc_flux,
            q_ghz = q_values[2] - q_values[1],
            tc_ghz = tc_values[2] - tc_values[1],
            fan_abs_indices = fan_abs_indices,
            fan_ghz = fan_ghz,
            fan_overlaps_a = fan_overlaps_a,
            fan_overlaps_b = fan_overlaps_b,
            fan_total_support = fan_total_support,
            which = which,
            branch = branch,
        ))
    end

    return rows
end

function adjacent_pair_metrics(row; gap_bias = 0.05)
    metrics = NamedTuple[]
    length(row.fan_ghz) >= 2 || error("Need at least two excited levels to define an avoided crossing pair.")

    for lower_rel in 1:(length(row.fan_ghz) - 1)
        upper_rel = lower_rel + 1
        gap_ghz = row.fan_ghz[upper_rel] - row.fan_ghz[lower_rel]
        a_pair_support = row.fan_overlaps_a[lower_rel] + row.fan_overlaps_a[upper_rel]
        b_pair_support = row.fan_overlaps_b[lower_rel] + row.fan_overlaps_b[upper_rel]
        lower_branch_support = row.fan_total_support[lower_rel]
        upper_branch_support = row.fan_total_support[upper_rel]
        total_support = a_pair_support + b_pair_support
        balance_bonus = 0.1 * min(a_pair_support, b_pair_support)
        order_support = max(
            row.fan_overlaps_a[lower_rel] + row.fan_overlaps_b[upper_rel],
            row.fan_overlaps_b[lower_rel] + row.fan_overlaps_a[upper_rel],
        )
        score = total_support + balance_bonus + 0.05 * order_support - gap_bias * gap_ghz

        push!(metrics, (
            lower_rel = lower_rel,
            upper_rel = upper_rel,
            pair_rel = (lower_rel, upper_rel),
            pair_abs = (row.fan_abs_indices[lower_rel], row.fan_abs_indices[upper_rel]),
            gap_ghz = gap_ghz,
            a_pair_support = a_pair_support,
            b_pair_support = b_pair_support,
            lower_branch_support = lower_branch_support,
            upper_branch_support = upper_branch_support,
            order_support = order_support,
            score = score,
        ))
    end

    return metrics
end

function choose_highlighted_pair(rows; gap_bias = 0.05)
    best = nothing

    for (row_index, row) in enumerate(rows)
        for metric in adjacent_pair_metrics(row; gap_bias)
            candidate = merge(metric, (
                row_index = row_index,
                q_flux = row.q_flux,
                tc_flux = row.tc_flux,
                q_ghz = row.q_ghz,
                tc_ghz = row.tc_ghz,
            ))
            if best === nothing || candidate.score > best.score
                best = candidate
            end
        end
    end

    best === nothing && error("Could not identify an adjacent adiabatic pair to highlight.")
    return best
end

function adiabatic_fan_curves(rows; ylims = nothing)
    q_values = [row.q_ghz for row in rows]
    fan_count = maximum(length(row.fan_ghz) for row in rows)
    curves = NamedTuple[]

    for idx in 1:fan_count
        y_values = Float64[]
        for row in rows
            push!(y_values, idx <= length(row.fan_ghz) ? row.fan_ghz[idx] : NaN)
        end
        finite_values = collect(filter(isfinite, y_values))
        isempty(finite_values) && continue
        if ylims !== nothing
            ymin, ymax = ylims
            any(ymin <= value <= ymax for value in finite_values) || continue
        end
        push!(curves, (x = q_values, y = y_values, label = nothing, color = (:gray70, 0.9), linewidth = 1.75))
    end

    return curves
end

fan_y_values(curves) = reduce(vcat, [collect(filter(isfinite, curve.y)) for curve in curves]; init = Float64[])

function highlighted_pair_curves(rows, pair_choice)
    lower_rel, upper_rel = pair_choice.pair_rel
    lower = [row.fan_ghz[lower_rel] for row in rows]
    upper = [row.fan_ghz[upper_rel] for row in rows]
    return lower, upper
end

function third_mode_support_summary(row, pair_choice)
    lower_rel, upper_rel = pair_choice.pair_rel
    other_rels = [idx for idx in eachindex(row.fan_total_support) if idx != lower_rel && idx != upper_rel]

    if isempty(other_rels)
        return (
            rel_index = missing,
            abs_index = missing,
            support = 0.0,
            comparable = false,
        )
    end

    supports = [row.fan_total_support[idx] for idx in other_rels]
    best_position = argmax(supports)
    best_rel = other_rels[best_position]
    weaker_highlighted = min(row.fan_total_support[lower_rel], row.fan_total_support[upper_rel])
    third_support = supports[best_position]
    comparable = weaker_highlighted > 0 ? third_support >= 0.8 * weaker_highlighted : third_support > 0.1

    return (
        rel_index = best_rel,
        abs_index = row.fan_abs_indices[best_rel],
        support = third_support,
        comparable = comparable,
    )
end

function adiabatic_gap_summary(rows, pair_choice; gap_bias = 0.05)
    lower_rel, upper_rel = pair_choice.pair_rel
    gaps = [row.fan_ghz[upper_rel] - row.fan_ghz[lower_rel] for row in rows]
    best_index = argmin(gaps)
    best_row = rows[best_index]
    gap_ghz = gaps[best_index]
    a_pair_support = best_row.fan_overlaps_a[lower_rel] + best_row.fan_overlaps_a[upper_rel]
    b_pair_support = best_row.fan_overlaps_b[lower_rel] + best_row.fan_overlaps_b[upper_rel]
    order_support = max(
        best_row.fan_overlaps_a[lower_rel] + best_row.fan_overlaps_b[upper_rel],
        best_row.fan_overlaps_b[lower_rel] + best_row.fan_overlaps_a[upper_rel],
    )
    support_score = a_pair_support + b_pair_support + 0.1 * min(a_pair_support, b_pair_support) + 0.05 * order_support - gap_bias * gap_ghz
    third_mode = third_mode_support_summary(best_row, pair_choice)

    return (
        q_flux = best_row.q_flux,
        q_ghz = best_row.q_ghz,
        tc_flux = best_row.tc_flux,
        tc_ghz = best_row.tc_ghz,
        pair_abs = pair_choice.pair_abs,
        gap_ghz = gap_ghz,
        lower_a_overlap = best_row.fan_overlaps_a[lower_rel],
        lower_b_overlap = best_row.fan_overlaps_b[lower_rel],
        upper_a_overlap = best_row.fan_overlaps_a[upper_rel],
        upper_b_overlap = best_row.fan_overlaps_b[upper_rel],
        a_pair_support = a_pair_support,
        b_pair_support = b_pair_support,
        lower_branch_support = best_row.fan_total_support[lower_rel],
        upper_branch_support = best_row.fan_total_support[upper_rel],
        support_score = support_score,
        third_mode_abs_index = third_mode.abs_index,
        third_mode_support = third_mode.support,
        third_mode_comparable = third_mode.comparable,
    )
end

function adiabatic_validation_summary(rows, pair_choice, gap_summary)
    lower_rel, upper_rel = pair_choice.pair_rel
    gaps = [row.fan_ghz[upper_rel] - row.fan_ghz[lower_rel] for row in rows]
    best_index = argmin(gaps)
    expected_gap = rows[best_index].fan_ghz[upper_rel] - rows[best_index].fan_ghz[lower_rel]
    return (
        sorted_energies_monotone = all(issorted(row.fan_ghz) for row in rows),
        highlighted_pair_adjacent = pair_choice.pair_abs[2] == pair_choice.pair_abs[1] + 1,
        pair_indices_constant = all(
            length(row.fan_ghz) >= upper_rel &&
            row.fan_abs_indices[lower_rel] == pair_choice.pair_abs[1] &&
            row.fan_abs_indices[upper_rel] == pair_choice.pair_abs[2]
            for row in rows
        ),
        highlighted_noncrossing = all(row.fan_ghz[upper_rel] >= row.fan_ghz[lower_rel] for row in rows),
        reported_gap_matches = isapprox(gap_summary.gap_ghz, expected_gap; atol = 1e-12),
    )
end

function local_projector_qlevel(snapshot, q_flux, tc_flux, qlevel; branch = :move, charge_cutoff, qb_ncut = 5, tc_ncut = 4, resonator_dim = 3)
    _, qstates = q_bundle(snapshot, q_flux; branch, charge_cutoff, qb_ncut)
    _, tcstates = tc_bundle(snapshot, tc_flux; branch, charge_cutoff, tc_ncut)
    _, crstates = cr_bundle(snapshot; resonator_dim)
    Pq = qstates[qlevel + 1] * dag(qstates[qlevel + 1])
    Itc = sum(ψ * dag(ψ) for ψ in tcstates)
    Icr = sum(ψ * dag(ψ) for ψ in crstates)
    return kron(Pq, Itc, Icr)
end

function qubit_excited_scan(base_system, snapshot, q_fluxes, tc_fluxes, init_levels, t_gate; branch = :move, charge_cutoff, qb_ncut = 5, tc_ncut = 4, resonator_dim = 3)
    x_values = [qubit_frequency(snapshot, q_flux; branch, charge_cutoff, qb_ncut) for q_flux in q_fluxes]
    y_values = [coupler_frequency(snapshot, tc_flux; branch, charge_cutoff, tc_ncut) for tc_flux in tc_fluxes]
    matrix = zeros(length(tc_fluxes), length(q_fluxes))
    spec = local_spec(; charge_cutoff)

    for (j, tc_flux) in enumerate(tc_fluxes), (i, q_flux) in enumerate(q_fluxes)
        system = retuned_system(base_system, snapshot; branch, q_flux, tc_flux)
        model = build_model(system; hamiltonian_spec = spec)
        ψ0 = diabatic_state(snapshot, q_flux, tc_flux, init_levels...; branch, charge_cutoff, qb_ncut, tc_ncut, resonator_dim)
        result = evolve(model, ψ0, [0.0, t_gate])
        ψf = result.states[end]
        Pq = local_projector_qlevel(snapshot, q_flux, tc_flux, 1; branch, charge_cutoff, qb_ncut, tc_ncut, resonator_dim)
        matrix[j, i] = real(dag(ψf) * Pq * ψf)
    end

    return x_values, y_values, matrix
end

function ordered_heatmap_axes(x_values, y_values, matrix)
    x_order = sortperm(x_values)
    y_order = sortperm(y_values)
    return x_values[x_order], y_values[y_order], matrix[y_order, x_order]
end

function best_scan_point(snapshot, q_fluxes, tc_fluxes, matrix; branch = :move, objective = :min, charge_cutoff, qb_ncut = 5, tc_ncut = 4)
    linear_index = objective == :max ? argmax(matrix) : argmin(matrix)
    index = CartesianIndices(matrix)[linear_index]
    row, col = Tuple(index)
    q_flux = q_fluxes[col]
    tc_flux = tc_fluxes[row]
    return (
        q_flux = q_flux,
        tc_flux = tc_flux,
        q_ghz = qubit_frequency(snapshot, q_flux; branch, charge_cutoff, qb_ncut),
        tc_ghz = coupler_frequency(snapshot, tc_flux; branch, charge_cutoff, tc_ncut),
        value = matrix[row, col],
        q_edge = col == 1 ? :low : (col == length(q_fluxes) ? :high : :interior),
        tc_edge = row == 1 ? :low : (row == length(tc_fluxes) ? :high : :interior),
    )
end

function shift_flux_range(values, direction::Symbol; fraction = 0.6, min_flux = 0.0, max_flux = 0.5)
    direction == :interior && return values
    lo = first(values)
    hi = last(values)
    width = hi - lo
    signed_shift = direction == :high ? fraction * width : -fraction * width
    shifted_lo = clamp(lo + signed_shift, min_flux, max_flux - width)
    shifted_hi = shifted_lo + width
    return collect(range(shifted_lo, shifted_hi; length = length(values)))
end

function adaptive_qubit_excited_scan(base_system, snapshot, q_fluxes, tc_fluxes, init_levels, t_gate; branch = :move, objective = :min, charge_cutoff, max_rounds = 4, qb_ncut = 5, tc_ncut = 4, resonator_dim = 3)
    current_q_fluxes = q_fluxes
    current_tc_fluxes = tc_fluxes
    result = nothing

    for round in 1:max_rounds
        x_values, y_values, matrix = qubit_excited_scan(
            base_system,
            snapshot,
            current_q_fluxes,
            current_tc_fluxes,
            init_levels,
            t_gate;
            branch,
            charge_cutoff,
            qb_ncut,
            tc_ncut,
            resonator_dim,
        )
        best_point = best_scan_point(snapshot, current_q_fluxes, current_tc_fluxes, matrix; branch, objective, charge_cutoff, qb_ncut, tc_ncut)
        result = (
            q_fluxes = current_q_fluxes,
            tc_fluxes = current_tc_fluxes,
            x_values = x_values,
            y_values = y_values,
            matrix = matrix,
            best_point = best_point,
            rounds = round,
        )
        (best_point.q_edge == :interior && best_point.tc_edge == :interior) && return result
        round == max_rounds && return result
        current_q_fluxes = shift_flux_range(current_q_fluxes, best_point.q_edge)
        current_tc_fluxes = shift_flux_range(current_tc_fluxes, best_point.tc_edge)
    end

    return result
end

function quench_qubit_excited_population(base_system, snapshot, q_flux, tc_flux, init_levels, t_gate; branch = :move, charge_cutoff, qb_ncut = 5, tc_ncut = 4, resonator_dim = 3)
    system = retuned_system(base_system, snapshot; branch, q_flux, tc_flux)
    model = build_model(system; hamiltonian_spec = local_spec(; charge_cutoff))
    ψ0 = diabatic_state(snapshot, q_flux, tc_flux, init_levels...; branch, charge_cutoff, qb_ncut, tc_ncut, resonator_dim)
    result = evolve(model, ψ0, [0.0, t_gate])
    ψf = result.states[end]
    Pq = local_projector_qlevel(snapshot, q_flux, tc_flux, 1; branch, charge_cutoff, qb_ncut, tc_ncut, resonator_dim)
    return Float64(real(dag(ψf) * Pq * ψf))
end

final_qubit_excited_population(args...; kwargs...) = quench_qubit_excited_population(args...; kwargs...)

function computational_phase_summary(base_system, snapshot, q_flux, tc_flux, t_gate; branch = :cz, charge_cutoff)
    system = retuned_system(base_system, snapshot; branch, q_flux, tc_flux)
    spec = subspace_spec(
        system;
        hamiltonian_spec = local_spec(; charge_cutoff),
        subsystem_levels = (QB = [0, 1], CR = [0, 1]),
        basis = :dressed_static,
    )
    trace = projected_unitary(system, spec, [0.0, t_gate]; hamiltonian_spec = local_spec(; charge_cutoff))
    return (
        q_flux = q_flux,
        tc_flux = tc_flux,
        q_ghz = qubit_frequency(snapshot, q_flux; branch, charge_cutoff),
        tc_ghz = coupler_frequency(snapshot, tc_flux; branch, charge_cutoff),
        conditional_phase_pi = conditional_phase(trace.unitaries[end]) / pi,
        max_leakage = maximum(trace.leakages),
    )
end

function centered_span_window(best_point, x_values, y_values; delta_q_ghz, delta_tc_ghz)
    xlo = maximum((minimum(x_values), best_point.q_ghz - delta_q_ghz / 2))
    xhi = minimum((maximum(x_values), best_point.q_ghz + delta_q_ghz / 2))
    ylo = maximum((minimum(y_values), best_point.tc_ghz - delta_tc_ghz / 2))
    yhi = minimum((maximum(y_values), best_point.tc_ghz + delta_tc_ghz / 2))
    return (xlims = (xlo, xhi), ylims = (ylo, yhi))
end
