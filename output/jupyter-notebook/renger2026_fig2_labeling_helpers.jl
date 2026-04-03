fig2_round4(x::Real) = round(Float64(x); digits = 4)
fig2_compact_value(x) = x isa Real ? fig2_round4(x) : x
fig2_compact_namedtuple(nt::NamedTuple) = (; (key => fig2_compact_value(value) for (key, value) in pairs(nt))...)

const FIG2_LABELING_COLORS = (
    diabatic_a = :gray40,
    diabatic_b = :gray65,
    dressed_upper = :royalblue3,
    dressed_lower = :firebrick3,
)

const FIG2_LABELING_BRANCHES = (
    move = (
        label = "QB1-TC1-CR",
        qb_key = :QB1,
        tc_key = :TC1,
        beta_qc_key = "beta_qc_qb1",
        target_labels = ((1, 0, 0), (0, 0, 1)),
        target_names = ("|100> branch", "|001> branch"),
        q_flux_range = (0.21, 0.25),
        tc_flux = 0.34,
        qubit_ncut = 5,
        coupler_ncut = 4,
        resonator_dim = 3,
        levels = 6,
    ),
    cz = (
        label = "QB2-TC2-CR",
        qb_key = :QB2,
        tc_key = :TC2,
        beta_qc_key = "beta_qc_qb2",
        target_labels = ((1, 0, 1), (2, 0, 0)),
        target_names = ("|101> branch", "|200> branch"),
        q_flux_range = (0.16, 0.24),
        tc_flux = 0.34,
        qubit_ncut = 5,
        coupler_ncut = 4,
        resonator_dim = 3,
        levels = 6,
    ),
)

function fig2_branch_config(branch::Symbol)
    branch == :move && return FIG2_LABELING_BRANCHES.move
    branch == :cz && return FIG2_LABELING_BRANCHES.cz
    throw(ArgumentError("Unsupported Fig. 2 branch $branch. Expected :move or :cz."))
end

function fig2_snapshot_paths(repo_root::AbstractString)
    baseline_path = joinpath(repo_root, "output", "renger2026", "paper_local_priors.toml")
    overlay_candidate = joinpath(repo_root, "output", "renger2026", "fig2_ef_retune_working.toml")
    overlay_path = isfile(overlay_candidate) ? overlay_candidate : nothing
    return (; baseline_path, overlay_path)
end

function load_fig2_labeling_snapshot(repo_root::AbstractString)
    paths = fig2_snapshot_paths(repo_root)
    return load_renger2026_snapshot(paths.baseline_path; overlay_path = paths.overlay_path)
end

function fig2_branch_defaults(branch::Symbol = :move; q_points::Integer = 17, charge_cutoff::Integer = 10)
    q_points > 1 || throw(ArgumentError("q_points must be greater than 1."))
    charge_cutoff > 0 || throw(ArgumentError("charge_cutoff must be positive."))
    config = fig2_branch_config(branch)
    q_fluxes = collect(range(config.q_flux_range[1], config.q_flux_range[2]; length = q_points))
    return (
        branch = branch,
        branch_label = config.label,
        q_fluxes = q_fluxes,
        tc_flux = config.tc_flux,
        target_labels = config.target_labels,
        target_names = config.target_names,
        levels = config.levels,
        subsystem_levels = Dict(:QB => 3, :TC => 2, :CR => 2),
        charge_cutoff = charge_cutoff,
        qubit_ncut = config.qubit_ncut,
        coupler_ncut = config.coupler_ncut,
        resonator_dim = config.resonator_dim,
    )
end

function fig2_branch_summary(snapshot; branch::Symbol = :move, defaults = fig2_branch_defaults(branch))
    config = fig2_branch_config(branch)
    qb = snapshot.devices[config.qb_key]
    tc = snapshot.devices[config.tc_key]
    return (
        branch = branch,
        branch_label = config.label,
        qb_device = String(config.qb_key),
        tc_device = String(config.tc_key),
        target_labels = defaults.target_labels,
        q_flux_window = (first(defaults.q_fluxes), last(defaults.q_fluxes)),
        tc_flux = defaults.tc_flux,
        qb_f01_ghz = Float64(qb["f01_ghz"]),
        tc_f01_ghz = Float64(tc["f01_ghz"]),
        cr_f01_ghz = Float64(snapshot.targets["cr_f01_ghz"]),
        beta_qc = Float64(snapshot.targets[config.beta_qc_key]),
        beta_cr = Float64(snapshot.targets["beta_cr"]),
        charge_cutoff = defaults.charge_cutoff,
    )
end

function fig2_branch_system(
    snapshot;
    branch::Symbol = :move,
    q_flux = nothing,
    tc_flux = nothing,
    qubit_ncut::Integer = fig2_branch_config(branch).qubit_ncut,
    coupler_ncut::Integer = fig2_branch_config(branch).coupler_ncut,
    resonator_dim::Integer = fig2_branch_config(branch).resonator_dim,
)
    config = fig2_branch_config(branch)
    qb_data = snapshot.devices[config.qb_key]
    tc_data = snapshot.devices[config.tc_key]
    q_flux_value = isnothing(q_flux) ? Float64(qb_data["flux"]) : Float64(q_flux)
    tc_flux_value = isnothing(tc_flux) ? Float64(tc_data["flux"]) : Float64(tc_flux)

    qb = TunableTransmon(
        :QB;
        EJmax = Float64(qb_data["EJmax"]),
        EC = Float64(qb_data["EC"]),
        flux = q_flux_value,
        asymmetry = Float64(qb_data["asymmetry"]),
        ng = Float64(qb_data["ng"]),
        ncut = Int(qubit_ncut),
    )
    tc = TunableCoupler(
        :TC;
        EJmax = Float64(tc_data["EJmax"]),
        EC = Float64(tc_data["EC"]),
        flux = tc_flux_value,
        asymmetry = Float64(tc_data["asymmetry"]),
        ng = Float64(tc_data["ng"]),
        ncut = Int(coupler_ncut),
    )
    cr = Resonator(:CR; ω = Float64(snapshot.targets["cr_f01_ghz"]), dim = Int(resonator_dim))

    return CompositeSystem(
        qb,
        tc,
        cr,
        CircuitCapacitiveCoupling(:QB, :TC; G = Float64(snapshot.targets[config.beta_qc_key])),
        CircuitCapacitiveCoupling(:TC, :CR; G = Float64(snapshot.targets["beta_cr"])),
    )
end

function fig2_qubit_frequency_ghz(snapshot, q_flux; branch::Symbol = :move, charge_cutoff::Integer = 10, qubit_ncut::Integer = 5)
    config = fig2_branch_config(branch)
    qb_data = snapshot.devices[config.qb_key]
    qb = TunableTransmon(
        :QB;
        EJmax = Float64(qb_data["EJmax"]),
        EC = Float64(qb_data["EC"]),
        flux = Float64(q_flux),
        asymmetry = Float64(qb_data["asymmetry"]),
        ng = Float64(qb_data["ng"]),
        ncut = Int(qubit_ncut),
    )
    esys = eigensystem(CompositeSystem(qb); levels = 2, hamiltonian_spec = CircuitHamiltonianSpec(charge_cutoff = charge_cutoff))
    return esys.energies[2] - esys.energies[1]
end

function fig2_qubit_frequency_axis(snapshot, q_fluxes; branch::Symbol = :move, charge_cutoff::Integer = 10, qubit_ncut::Integer = 5)
    return [fig2_qubit_frequency_ghz(snapshot, q_flux; branch, charge_cutoff, qubit_ncut) for q_flux in q_fluxes]
end

function fig2_labeled_crossing_curves(
    x_values,
    labeled::LabeledEigensystemSweepResult;
    branch::Symbol = :move,
    bare_labels = fig2_branch_config(branch).target_labels,
    bare_label_names = fig2_branch_config(branch).target_names,
    subtract_ground::Bool = true,
)
    length(x_values) == length(labeled.sweep.values) ||
        throw(ArgumentError("x_values must align with the labeled sweep values."))

    curves = NamedTuple[]
    palette = (FIG2_LABELING_COLORS.diabatic_a, FIG2_LABELING_COLORS.diabatic_b)

    for index in eachindex(bare_labels)
        push!(
            curves,
            (
                x = x_values,
                y = energy_curve(labeled, bare_labels[index]; subtract_ground = subtract_ground).data,
                label = bare_label_names[index],
                color = palette[index],
                linestyle = :dash,
                linewidth = 2,
            ),
        )
    end

    return curves
end

function fig2_labeled_crossing_figure(
    x_values,
    labeled::LabeledEigensystemSweepResult;
    branch::Symbol = :move,
    bare_labels = fig2_branch_config(branch).target_labels,
    bare_label_names = fig2_branch_config(branch).target_names,
    title::AbstractString = "Fig. 2 move branch: bare-labeled dressed curves",
    xlabel::AbstractString = "omega_QB,01 (GHz)",
    ylabel::AbstractString = "E - E_ground (GHz)",
    subtract_ground::Bool = true,
)
    curves = fig2_labeled_crossing_curves(
        x_values,
        labeled;
        branch,
        bare_labels,
        bare_label_names,
        subtract_ground,
    )
    return line_figure(curves; title, xlabel, ylabel, size = NOTEBOOK_WIDE, legend_position = :rb)
end

function fig2_adiabatic_pair_curves(
    labeled::LabeledEigensystemSweepResult;
    pair_abs,
    subtract_ground::Bool = true,
)
    length(pair_abs) == 2 || throw(ArgumentError("pair_abs must contain exactly two dressed-state indices."))
    lower_index, upper_index = sort(Tuple(Int(index) for index in pair_abs))

    lower_curve = Float64[]
    upper_curve = Float64[]
    for eigensystem_result in labeled.sweep.spectra
        1 <= lower_index <= length(eigensystem_result.energies) ||
            throw(ArgumentError("lower_index $lower_index is out of bounds for the available eigensystem levels."))
        1 <= upper_index <= length(eigensystem_result.energies) ||
            throw(ArgumentError("upper_index $upper_index is out of bounds for the available eigensystem levels."))

        ground = subtract_ground ? eigensystem_result.energies[1] : 0.0
        push!(lower_curve, eigensystem_result.energies[lower_index] - ground)
        push!(upper_curve, eigensystem_result.energies[upper_index] - ground)
    end

    return (
        lower = lower_curve,
        upper = upper_curve,
        pair_abs = (lower_index, upper_index),
    )
end

function fig2_adiabatic_crossing_figure(
    x_values,
    labeled::LabeledEigensystemSweepResult;
    pair_abs,
    branch::Symbol = :move,
    bare_labels = fig2_branch_config(branch).target_labels,
    bare_label_names = fig2_branch_config(branch).target_names,
    title::AbstractString = "Fig. 2 move branch via adiabatic avoided-crossing pair",
    xlabel::AbstractString = "omega_QB,01 (GHz)",
    ylabel::AbstractString = "E - E_ground (GHz)",
    subtract_ground::Bool = true,
    include_bare_guides::Bool = true,
    x_center = nothing,
    x_half_width = nothing,
)
    length(x_values) == length(labeled.sweep.values) ||
        throw(ArgumentError("x_values must align with the labeled sweep values."))

    adiabatic = fig2_adiabatic_pair_curves(labeled; pair_abs, subtract_ground)
    curves = NamedTuple[]

    if include_bare_guides
        append!(
            curves,
            fig2_labeled_crossing_curves(
                x_values,
                labeled;
                branch,
                bare_labels,
                bare_label_names,
                subtract_ground,
            ),
        )
    end

    push!(
        curves,
        (
            x = x_values,
            y = adiabatic.upper,
            label = "upper dressed branch",
            color = FIG2_LABELING_COLORS.dressed_upper,
        ),
    )
    push!(
        curves,
        (
            x = x_values,
            y = adiabatic.lower,
            label = "lower dressed branch",
            color = FIG2_LABELING_COLORS.dressed_lower,
        ),
    )

    xlims = nothing
    ylims = nothing
    if !isnothing(x_center) && !isnothing(x_half_width)
        x_center = Float64(x_center)
        x_half_width = Float64(x_half_width)
        xlims = (x_center - x_half_width, x_center + x_half_width)
        window_indices = findall(index -> xlims[1] <= x_values[index] <= xlims[2], eachindex(x_values))
        if !isempty(window_indices)
            window_y = vcat(adiabatic.lower[window_indices], adiabatic.upper[window_indices])
            y_pad = max(0.02, 0.15 * (maximum(window_y) - minimum(window_y)))
            ylims = (minimum(window_y) - y_pad, maximum(window_y) + y_pad)
        end
    end

    return line_figure(curves; title, xlabel, ylabel, size = NOTEBOOK_WIDE, legend_position = :rb, xlims, ylims)
end

function fig2_valid_labeled_indices(
    labeled::LabeledEigensystemSweepResult;
    bare_labels,
    subtract_ground::Bool = true,
)
    curves = [energy_curve(labeled, bare_label; subtract_ground = subtract_ground).data for bare_label in bare_labels]
    isempty(curves) && throw(ArgumentError("bare_labels must contain at least one label."))
    valid = Int[]
    for index in eachindex(curves[1])
        all(isfinite(curve[index]) for curve in curves) && push!(valid, index)
    end
    return valid
end

function fig2_best_labeled_crossing(
    labeled::LabeledEigensystemSweepResult,
    q_fluxes,
    q_frequency_axis;
    bare_labels,
    subtract_ground::Bool = true,
)
    length(q_fluxes) == length(q_frequency_axis) == length(labeled.sweep.values) ||
        throw(ArgumentError("q_fluxes, q_frequency_axis, and labeled sweep values must have the same length."))
    length(bare_labels) == 2 || throw(ArgumentError("fig2_best_labeled_crossing expects exactly two bare labels."))

    curve_a = energy_curve(labeled, bare_labels[1]; subtract_ground = subtract_ground).data
    curve_b = energy_curve(labeled, bare_labels[2]; subtract_ground = subtract_ground).data
    valid_indices = fig2_valid_labeled_indices(labeled; bare_labels, subtract_ground)
    isempty(valid_indices) &&
        throw(ArgumentError("No sweep points have all requested bare labels assigned under the current labeling threshold."))

    gap_values = [abs(curve_a[index] - curve_b[index]) for index in valid_indices]
    best_offset = argmin(gap_values)
    best_index = valid_indices[best_offset]
    dressed_index_a = dressed_index(labeled, bare_labels[1], best_index)
    dressed_index_b = dressed_index(labeled, bare_labels[2], best_index)
    if ismissing(dressed_index_a) || ismissing(dressed_index_b)
        throw(ArgumentError("Internal labeling inconsistency at best_index $best_index."))
    end
    pair_abs = Tuple(sort((Int(dressed_index_a), Int(dressed_index_b))))

    return (
        index = best_index,
        valid_indices = valid_indices,
        valid_count = length(valid_indices),
        dressed_index_a = Int(dressed_index_a),
        dressed_index_b = Int(dressed_index_b),
        pair_abs = pair_abs,
        q_flux = q_fluxes[best_index],
        q_ghz = q_frequency_axis[best_index],
        bare_energy_a = curve_a[best_index],
        bare_energy_b = curve_b[best_index],
        labeled_gap_ghz = gap_values[best_offset],
    )
end
