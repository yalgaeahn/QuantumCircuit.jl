using Pkg

function find_repo_root(start::AbstractString)
    candidates = [
        normpath(start),
        normpath(joinpath(start, "..")),
        normpath(joinpath(start, "..", "..")),
    ]

    for candidate in unique(candidates)
        project_toml = joinpath(candidate, "Project.toml")
        if isfile(project_toml)
            content = read(project_toml, String)
            occursin("QuantumCircuit", content) && return candidate
        end
    end

    error("Could not find the QuantumCircuit.jl repository root from $(start).")
end

repo_root = find_repo_root(pwd())
Pkg.activate(repo_root)
Pkg.instantiate()

using Printf
using TOML
using CairoMakie
using QuantumCircuit
using QuantumToolbox: basis, dag, eigenstates, kron, Ket, QuantumObject

include(joinpath(repo_root, "output", "jupyter-notebook", "makie_helpers.jl"))
activate_notebook_theme!()

const FIG2_CHARGE_CUTOFF = 3
const FIG2_MOVE_PLOT_TC_FLUX = 0.34
const FIG2_CZ_PLOT_TC_FLUX = 0.34
const FIG2_MOVE_T_GATE_DEFAULT = 100.0
const FIG2_CZ_T_GATE_DEFAULT = 120.0
const FIG2_MOVE_RAMP_NS = 4.0
const FIG2_CZ_RAMP_NS = 4.0
const FIG2_MOVE_OSC_TLIST = collect(0.0:2.0:200.0)
const FIG2_CZ_OSC_TLIST = collect(0.0:2.0:260.0)
const FIG2_CZ_RETURN_F_THRESHOLD = 0.15
const FIG2_OSC_T_MIN_NS = 10.0
const FIG2_ZOOM_Q_SPAN_GHZ = 0.10
const FIG2_MOVE_ZOOM_TC_SPAN_GHZ = 1.60
const FIG2_CZ_ZOOM_TC_SPAN_GHZ = 1.90
const FIG2_ZOOM_Q_POINTS = 61
const FIG2_ZOOM_TC_POINTS = 101
const RUN_TRUNCATION_CHECKS = false
figure_exports = Dict{Symbol, NamedTuple}()

nothing


overlay_candidate = joinpath(repo_root, "output", "renger2026", "fig2_ef_retune_working.toml")
const FIG2_OVERLAY_PATH = isfile(overlay_candidate) ? overlay_candidate : nothing
snapshot = load_renger2026_snapshot(; overlay_path = FIG2_OVERLAY_PATH)
include(joinpath(repo_root, "output", "jupyter-notebook", "renger2026_fig2_branch_helpers.jl"))

gate_times = fig2_gate_times(snapshot; move_default = FIG2_MOVE_T_GATE_DEFAULT, cz_default = FIG2_CZ_T_GATE_DEFAULT)
move_system = local_exact_system(
    snapshot;
    branch = :move,
    q_flux = branch_q_park_flux(snapshot, :move),
    tc_flux = branch_tc_park_flux(snapshot, :move),
)
cz_system = local_exact_system(
    snapshot;
    branch = :cz,
    q_flux = branch_q_park_flux(snapshot, :cz),
    tc_flux = branch_tc_park_flux(snapshot, :cz),
)

nothing


parameter_summary, coupling_summary = parameter_rows(
    snapshot;
    move_t_gate_ns = gate_times.move,
    cz_t_gate_ns = gate_times.cz,
    charge_cutoff = FIG2_CHARGE_CUTOFF,
)
branch_summary = (
    move_branch = branch_spec(:move).label,
    cz_branch = branch_spec(:cz).label,
    overlay_path = FIG2_OVERLAY_PATH,
)

(; branch_summary, parameter_summary, coupling_summary)


move_plot_tc_flux = FIG2_MOVE_PLOT_TC_FLUX
move_q_fluxes = collect(range(0.21, 0.25; length = 17))
move_rows = adiabatic_fan_rows(
    move_system,
    snapshot,
    move_q_fluxes,
    move_plot_tc_flux;
    branch_pair_a = (1, 0, 0),
    branch_pair_b = (0, 0, 1),
    which = :move,
    branch = :move,
    charge_cutoff = FIG2_CHARGE_CUTOFF,
)
move_pair = choose_highlighted_pair(move_rows)

move_x = [row.q_ghz for row in move_rows]
move_lower, move_upper = highlighted_pair_curves(move_rows, move_pair)
move_pair_span = maximum(move_upper) - minimum(move_lower)
move_y_pad = max(0.02, 0.15 * move_pair_span)
move_ylims = (minimum(move_lower) - move_y_pad, maximum(move_upper) + move_y_pad)
move_fan_curves = adiabatic_fan_curves(move_rows; ylims = move_ylims)
move_background_y = fan_y_values(move_fan_curves)
move_y_values = isempty(move_background_y) ? vcat(move_lower, move_upper) : vcat(move_lower, move_upper, move_background_y)
move_xlims = (minimum(move_x), maximum(move_x))
move_ylims = (minimum(move_y_values) - move_y_pad, maximum(move_y_values) + move_y_pad)

move_fig = line_figure(
    vcat(
        move_fan_curves,
        [
            (x = move_x, y = move_upper, label = "upper dressed branch", color = FIG2_COLORS.dressed_a),
            (x = move_x, y = move_lower, label = "lower dressed branch", color = FIG2_COLORS.dressed_b),
        ],
    );
    title = "Fig. 2(c) analog: MOVE avoided crossing",
    xlabel = "omega_QB,01 (GHz)",
    ylabel = "E - E_ground (GHz)",
    size = NOTEBOOK_WIDE,
    legend_position = :rb,
    xlims = move_xlims,
    ylims = move_ylims,
)
move_saved = save_figure(move_fig, repo_root, "figure2c_move_avoided_crossing")
figure_exports[:figure2c] = move_saved

display(move_fig)

move_gap = adiabatic_gap_summary(move_rows, move_pair)
move_validation = adiabatic_validation_summary(move_rows, move_pair, move_gap)
move_reference_third_mode = third_mode_support_summary(move_rows[move_pair.row_index], move_pair)

(
    saved = move_saved,
    plotted_tc_flux = round4(move_plot_tc_flux),
    highlighted_pair = move_pair.pair_abs,
    move_reference = compact_namedtuple((
        q_flux = move_pair.q_flux,
        q_ghz = move_pair.q_ghz,
        tc_flux = move_pair.tc_flux,
        tc_ghz = move_pair.tc_ghz,
        pair_abs = move_pair.pair_abs,
        gap_ghz = move_pair.gap_ghz,
        a_pair_support = move_pair.a_pair_support,
        b_pair_support = move_pair.b_pair_support,
        score = move_pair.score,
    )),
    move_reference_third_mode = compact_namedtuple(move_reference_third_mode),
    move_gap = compact_namedtuple(move_gap),
    move_validation = move_validation,
)


cz_plot_tc_flux = FIG2_CZ_PLOT_TC_FLUX
cz_q_fluxes = collect(range(0.16, 0.24; length = 17))
cz_rows = adiabatic_fan_rows(
    cz_system,
    snapshot,
    cz_q_fluxes,
    cz_plot_tc_flux;
    branch_pair_a = (1, 0, 1),
    branch_pair_b = (2, 0, 0),
    which = :cz,
    branch = :cz,
    charge_cutoff = FIG2_CHARGE_CUTOFF,
)
cz_pair = choose_highlighted_pair(cz_rows)

cz_x = [row.q_ghz for row in cz_rows]
cz_lower, cz_upper = highlighted_pair_curves(cz_rows, cz_pair)
cz_pair_span = maximum(cz_upper) - minimum(cz_lower)
cz_y_pad = max(0.02, 0.15 * cz_pair_span)
cz_ylims = (minimum(cz_lower) - cz_y_pad, maximum(cz_upper) + cz_y_pad)
cz_fan_curves = adiabatic_fan_curves(cz_rows; ylims = cz_ylims)
cz_background_y = fan_y_values(cz_fan_curves)
cz_y_values = isempty(cz_background_y) ? vcat(cz_lower, cz_upper) : vcat(cz_lower, cz_upper, cz_background_y)
cz_xlims = (minimum(cz_x), maximum(cz_x))
cz_ylims = (minimum(cz_y_values) - cz_y_pad, maximum(cz_y_values) + cz_y_pad)

cz_fig = line_figure(
    vcat(
        cz_fan_curves,
        [
            (x = cz_x, y = cz_upper, label = "upper dressed branch", color = FIG2_COLORS.dressed_a),
            (x = cz_x, y = cz_lower, label = "lower dressed branch", color = FIG2_COLORS.dressed_b),
        ],
    );
    title = "Fig. 2(d) analog: CZ-side avoided crossing",
    xlabel = "omega_QB,01 (GHz)",
    ylabel = "E - E_ground (GHz)",
    size = NOTEBOOK_WIDE,
    legend_position = :rb,
    xlims = cz_xlims,
    ylims = cz_ylims,
)
cz_saved = save_figure(cz_fig, repo_root, "figure2d_cz_avoided_crossing")
figure_exports[:figure2d] = cz_saved

display(cz_fig)

cz_gap = adiabatic_gap_summary(cz_rows, cz_pair)
cz_validation = adiabatic_validation_summary(cz_rows, cz_pair, cz_gap)
cz_reference_third_mode = third_mode_support_summary(cz_rows[cz_pair.row_index], cz_pair)

(
    saved = cz_saved,
    plotted_tc_flux = round4(cz_plot_tc_flux),
    highlighted_pair = cz_pair.pair_abs,
    cz_reference = compact_namedtuple((
        q_flux = cz_pair.q_flux,
        q_ghz = cz_pair.q_ghz,
        tc_flux = cz_pair.tc_flux,
        tc_ghz = cz_pair.tc_ghz,
        pair_abs = cz_pair.pair_abs,
        gap_ghz = cz_pair.gap_ghz,
        a_pair_support = cz_pair.a_pair_support,
        b_pair_support = cz_pair.b_pair_support,
        score = cz_pair.score,
    )),
    cz_reference_third_mode = compact_namedtuple(cz_reference_third_mode),
    cz_gap = compact_namedtuple(cz_gap),
    cz_validation = cz_validation,
)


move_scan_axes = branch_flux_axes(:move; q_points = 21, tc_points = 21)
move_scan_q_fluxes = move_scan_axes.q_fluxes
move_scan_tc_fluxes = move_scan_axes.tc_fluxes
move_t_gate = gate_times.move

move_osc_scan = adaptive_oscillation_metric_scan(
    move_system,
    snapshot,
    move_scan_q_fluxes,
    move_scan_tc_fluxes,
    (1, 0, 0),
    FIG2_MOVE_OSC_TLIST;
    branch = :move,
    metric = :move,
    charge_cutoff = FIG2_CHARGE_CUTOFF,
    t_min_ns = FIG2_OSC_T_MIN_NS,
)
move_osc_x, move_osc_y, move_osc_matrix = ordered_heatmap_axes(
    move_osc_scan.x_values,
    move_osc_scan.y_values,
    move_osc_scan.matrix,
)
move_osc_fig = heatmap_figure(
    move_osc_x,
    move_osc_y,
    move_osc_matrix;
    title = "Fig. 2(e) diagnostic: MOVE oscillation-envelope proxy",
    xlabel = "omega_QB,01 (GHz)",
    ylabel = "omega_TC (GHz)",
    colorlabel = "min_t P_e (MOVE qubit)",
    colorrange = (0.0, 1.0),
    xlims = (minimum(move_osc_x), maximum(move_osc_x)),
    ylims = (minimum(move_osc_y), maximum(move_osc_y)),
)
move_osc_saved = save_figure(move_osc_fig, repo_root, "figure2e_move_population_map_oscillation_proxy")
figure_exports[:figure2e_oscillation_proxy] = move_osc_saved

move_scan = adaptive_pulse_population_scan(
    move_system,
    snapshot,
    move_scan_q_fluxes,
    move_scan_tc_fluxes,
    (1, 0, 0),
    move_t_gate;
    branch = :move,
    metric = :move,
    objective = :max,
    ramp_ns = FIG2_MOVE_RAMP_NS,
    charge_cutoff = FIG2_CHARGE_CUTOFF,
)

x_move = move_scan.x_values
y_move = move_scan.y_values
matrix_move = move_scan.matrix
x_move_sorted, y_move_sorted, matrix_move_sorted = ordered_heatmap_axes(x_move, y_move, matrix_move)
best_move_point = move_scan.best_point
move_display = (
    xlims = (minimum(x_move_sorted), maximum(x_move_sorted)),
    ylims = (minimum(y_move_sorted), maximum(y_move_sorted)),
)

move_heatmap_fig = Figure(size = NOTEBOOK_WIDE)
ax_move = Axis(
    move_heatmap_fig[1, 1];
    title = "Fig. 2(e) analog: MOVE final qubit-state probability after pulse",
    xlabel = "omega_QB,01 (GHz)",
    ylabel = "omega_TC (GHz)",
)
move_heatmap = heatmap!(
    ax_move,
    x_move_sorted,
    y_move_sorted,
    matrix_move_sorted';
    colormap = :Reds,
    colorrange = (0.0, 1.0),
)
xlims!(ax_move, move_display.xlims...)
ylims!(ax_move, move_display.ylims...)
scatter!(
    ax_move,
    [best_move_point.q_ghz],
    [best_move_point.tc_ghz];
    color = :white,
    strokecolor = :black,
    markersize = 12,
)
Colorbar(
    move_heatmap_fig[1, 2],
    move_heatmap;
    label = "QB state probability",
    ticks = ([0.0, 1.0], ["|g⟩", "|e⟩"]),
)
move_heatmap_saved = save_figure(move_heatmap_fig, repo_root, "figure2e_move_population_map")
figure_exports[:figure2e] = move_heatmap_saved

move_zoom = dense_best_centered_pulse_population_scan(
    move_system,
    snapshot,
    move_scan,
    (1, 0, 0),
    move_t_gate;
    branch = :move,
    metric = :move,
    objective = :max,
    ramp_ns = FIG2_MOVE_RAMP_NS,
    charge_cutoff = FIG2_CHARGE_CUTOFF,
    q_span_ghz = FIG2_ZOOM_Q_SPAN_GHZ,
    tc_span_ghz = FIG2_MOVE_ZOOM_TC_SPAN_GHZ,
    q_points = FIG2_ZOOM_Q_POINTS,
    tc_points = FIG2_ZOOM_TC_POINTS,
)
move_zoom_fig = Figure(size = NOTEBOOK_WIDE)
ax_move_zoom = Axis(
    move_zoom_fig[1, 1];
    title = "Fig. 2(e) analog: dense best-centered MOVE pulse zoom",
    xlabel = "omega_QB,01 (GHz)",
    ylabel = "omega_TC (GHz)",
)
move_zoom_hm = heatmap!(
    ax_move_zoom,
    move_zoom.x_sorted,
    move_zoom.y_sorted,
    move_zoom.matrix_sorted';
    colormap = :Reds,
    colorrange = (0.0, 1.0),
)
xlims!(ax_move_zoom, move_zoom.display.xlims...)
ylims!(ax_move_zoom, move_zoom.display.ylims...)
scatter!(
    ax_move_zoom,
    [move_zoom.best_point.q_ghz],
    [move_zoom.best_point.tc_ghz];
    color = :white,
    strokecolor = :black,
    markersize = 12,
)
Colorbar(
    move_zoom_fig[1, 2],
    move_zoom_hm;
    label = "QB state probability",
    ticks = ([0.0, 1.0], ["|g⟩", "|e⟩"]),
)
move_zoom_saved = save_figure(move_zoom_fig, repo_root, "figure2e_move_population_map_best_zoom")
figure_exports[:figure2e_best_zoom] = move_zoom_saved

move_best_trace = pulse_population_trace(
    move_system,
    snapshot,
    best_move_point.q_flux,
    best_move_point.tc_flux,
    (1, 0, 0),
    move_t_gate;
    branch = :move,
    ramp_ns = FIG2_MOVE_RAMP_NS,
    charge_cutoff = FIG2_CHARGE_CUTOFF,
)
move_trace_fig = Figure(size = NOTEBOOK_TALL)
ax_move_trace = Axis(
    move_trace_fig[1, 1];
    title = "MOVE best-point pulse trace",
    xlabel = "time (ns)",
    ylabel = "probability",
)
lines!(ax_move_trace, move_best_trace.times, move_best_trace.p_e; color = :firebrick3, linewidth = 3, label = "P_e")
lines!(ax_move_trace, move_best_trace.times, 1 .- move_best_trace.p_e; color = :royalblue3, linewidth = 3, label = "P_g = 1 - P_e")
vlines!(ax_move_trace, [move_best_trace.final_time_ns]; color = :black, linestyle = :dot)
axislegend(ax_move_trace; position = :rb)
move_trace_saved = save_figure(move_trace_fig, repo_root, "figure2e_move_best_trace")
figure_exports[:figure2e_best_trace] = move_trace_saved

display(move_heatmap_fig)
display(move_zoom_fig)
display(move_trace_fig)

(
    oscillation_proxy_saved = move_osc_saved,
    main_saved = move_heatmap_saved,
    best_zoom_saved = move_zoom_saved,
    best_trace_saved = move_trace_saved,
    coarse_best_point = compact_namedtuple(best_move_point),
    dense_best_point = compact_namedtuple(move_zoom.best_point),
    selected_best_time_ns = round4(best_move_point.best_time_ns),
    final_p_g = round4(best_move_point.aux),
    requested_zoom = compact_namedtuple(move_zoom.requested_span_ghz),
    realized_zoom = compact_namedtuple(move_zoom.realized_span_ghz),
    dense_grid = move_zoom.grid_size,
)


cz_scan_axes = branch_flux_axes(:cz; q_points = 21, tc_points = 21)
cz_scan_q_fluxes = cz_scan_axes.q_fluxes
cz_scan_tc_fluxes = cz_scan_axes.tc_fluxes
cz_t_gate = gate_times.cz

cz_osc_scan = adaptive_oscillation_metric_scan(
    cz_system,
    snapshot,
    cz_scan_q_fluxes,
    cz_scan_tc_fluxes,
    (1, 0, 1),
    FIG2_CZ_OSC_TLIST;
    branch = :cz,
    metric = :cz,
    charge_cutoff = FIG2_CHARGE_CUTOFF,
    t_min_ns = FIG2_OSC_T_MIN_NS,
    f_threshold = FIG2_CZ_RETURN_F_THRESHOLD,
)
cz_osc_x, cz_osc_y, cz_osc_matrix = ordered_heatmap_axes(
    cz_osc_scan.x_values,
    cz_osc_scan.y_values,
    cz_osc_scan.matrix,
)
cz_osc_fig = heatmap_figure(
    cz_osc_x,
    cz_osc_y,
    cz_osc_matrix;
    title = "Fig. 2(f) diagnostic: CZ oscillation-envelope proxy",
    xlabel = "omega_QB,01 (GHz)",
    ylabel = "omega_TC (GHz)",
    colorlabel = "1 - max return score",
    colorrange = (0.0, 1.0),
    xlims = (minimum(cz_osc_x), maximum(cz_osc_x)),
    ylims = (minimum(cz_osc_y), maximum(cz_osc_y)),
)
cz_osc_saved = save_figure(cz_osc_fig, repo_root, "figure2f_cz_population_map_oscillation_proxy")
figure_exports[:figure2f_oscillation_proxy] = cz_osc_saved

cz_scan = adaptive_pulse_population_scan(
    cz_system,
    snapshot,
    cz_scan_q_fluxes,
    cz_scan_tc_fluxes,
    (1, 0, 1),
    cz_t_gate;
    branch = :cz,
    metric = :cz,
    objective = :max,
    ramp_ns = FIG2_CZ_RAMP_NS,
    charge_cutoff = FIG2_CHARGE_CUTOFF,
)

x_cz = cz_scan.x_values
y_cz = cz_scan.y_values
matrix_cz = cz_scan.matrix
x_cz_sorted, y_cz_sorted, matrix_cz_sorted = ordered_heatmap_axes(x_cz, y_cz, matrix_cz)
best_cz_point = cz_scan.best_point
cz_display = (
    xlims = (minimum(x_cz_sorted), maximum(x_cz_sorted)),
    ylims = (minimum(y_cz_sorted), maximum(y_cz_sorted)),
)

cz_heatmap_fig = Figure(size = NOTEBOOK_WIDE)
ax_cz = Axis(
    cz_heatmap_fig[1, 1];
    title = "Fig. 2(f) analog: CZ final qubit-state probability after pulse",
    xlabel = "omega_QB,01 (GHz)",
    ylabel = "omega_TC (GHz)",
)
cz_heatmap = heatmap!(
    ax_cz,
    x_cz_sorted,
    y_cz_sorted,
    matrix_cz_sorted';
    colormap = :Blues,
    colorrange = (0.0, 1.0),
)
xlims!(ax_cz, cz_display.xlims...)
ylims!(ax_cz, cz_display.ylims...)
scatter!(
    ax_cz,
    [best_cz_point.q_ghz],
    [best_cz_point.tc_ghz];
    color = :white,
    strokecolor = :black,
    markersize = 12,
)
Colorbar(
    cz_heatmap_fig[1, 2],
    cz_heatmap;
    label = "QB state probability",
    ticks = ([0.0, 1.0], ["|e⟩", "|f⟩"]),
)
cz_heatmap_saved = save_figure(cz_heatmap_fig, repo_root, "figure2f_cz_population_map")
figure_exports[:figure2f] = cz_heatmap_saved

cz_zoom = dense_best_centered_pulse_population_scan(
    cz_system,
    snapshot,
    cz_scan,
    (1, 0, 1),
    cz_t_gate;
    branch = :cz,
    metric = :cz,
    objective = :max,
    ramp_ns = FIG2_CZ_RAMP_NS,
    charge_cutoff = FIG2_CHARGE_CUTOFF,
    q_span_ghz = FIG2_ZOOM_Q_SPAN_GHZ,
    tc_span_ghz = FIG2_CZ_ZOOM_TC_SPAN_GHZ,
    q_points = FIG2_ZOOM_Q_POINTS,
    tc_points = FIG2_ZOOM_TC_POINTS,
)
cz_zoom_fig = Figure(size = NOTEBOOK_WIDE)
ax_cz_zoom = Axis(
    cz_zoom_fig[1, 1];
    title = "Fig. 2(f) analog: dense best-centered CZ pulse zoom",
    xlabel = "omega_QB,01 (GHz)",
    ylabel = "omega_TC (GHz)",
)
cz_zoom_hm = heatmap!(
    ax_cz_zoom,
    cz_zoom.x_sorted,
    cz_zoom.y_sorted,
    cz_zoom.matrix_sorted';
    colormap = :Blues,
    colorrange = (0.0, 1.0),
)
xlims!(ax_cz_zoom, cz_zoom.display.xlims...)
ylims!(ax_cz_zoom, cz_zoom.display.ylims...)
scatter!(
    ax_cz_zoom,
    [cz_zoom.best_point.q_ghz],
    [cz_zoom.best_point.tc_ghz];
    color = :white,
    strokecolor = :black,
    markersize = 12,
)
Colorbar(
    cz_zoom_fig[1, 2],
    cz_zoom_hm;
    label = "QB state probability",
    ticks = ([0.0, 1.0], ["|e⟩", "|f⟩"]),
)
cz_zoom_saved = save_figure(cz_zoom_fig, repo_root, "figure2f_cz_population_map_best_zoom")
figure_exports[:figure2f_best_zoom] = cz_zoom_saved

cz_best_trace = pulse_population_trace(
    cz_system,
    snapshot,
    best_cz_point.q_flux,
    best_cz_point.tc_flux,
    (1, 0, 1),
    cz_t_gate;
    branch = :cz,
    ramp_ns = FIG2_CZ_RAMP_NS,
    charge_cutoff = FIG2_CHARGE_CUTOFF,
)
cz_trace_fig = Figure(size = NOTEBOOK_TALL)
ax_cz_trace = Axis(
    cz_trace_fig[1, 1];
    title = "CZ best-point pulse trace",
    xlabel = "time (ns)",
    ylabel = "probability",
)
lines!(ax_cz_trace, cz_best_trace.times, cz_best_trace.p_e; color = :royalblue3, linewidth = 3, label = "P_e")
lines!(ax_cz_trace, cz_best_trace.times, cz_best_trace.p_f; color = :goldenrod3, linewidth = 3, label = "P_f")
vlines!(ax_cz_trace, [cz_best_trace.final_time_ns]; color = :black, linestyle = :dot)
axislegend(ax_cz_trace; position = :rb)
cz_trace_saved = save_figure(cz_trace_fig, repo_root, "figure2f_cz_best_trace")
figure_exports[:figure2f_best_trace] = cz_trace_saved

display(cz_heatmap_fig)
display(cz_zoom_fig)
display(cz_trace_fig)

(
    oscillation_proxy_saved = cz_osc_saved,
    main_saved = cz_heatmap_saved,
    best_zoom_saved = cz_zoom_saved,
    best_trace_saved = cz_trace_saved,
    coarse_best_point = compact_namedtuple(best_cz_point),
    dense_best_point = compact_namedtuple(cz_zoom.best_point),
    selected_best_time_ns = round4(best_cz_point.best_time_ns),
    final_p_e = round4(best_cz_point.aux),
    requested_zoom = compact_namedtuple(cz_zoom.requested_span_ghz),
    realized_zoom = compact_namedtuple(cz_zoom.realized_span_ghz),
    dense_grid = cz_zoom.grid_size,
)


move_pulse_summary = final_pulse_state_probabilities(
    move_system,
    snapshot,
    best_move_point.q_flux,
    best_move_point.tc_flux,
    (1, 0, 0),
    move_t_gate;
    branch = :move,
    ramp_ns = FIG2_MOVE_RAMP_NS,
    charge_cutoff = FIG2_CHARGE_CUTOFF,
)

cz_pulse_summary = final_pulse_state_probabilities(
    cz_system,
    snapshot,
    best_cz_point.q_flux,
    best_cz_point.tc_flux,
    (1, 0, 1),
    cz_t_gate;
    branch = :cz,
    ramp_ns = FIG2_CZ_RAMP_NS,
    charge_cutoff = FIG2_CHARGE_CUTOFF,
)

cz_phase_proxy = pulse_computational_phase_summary(
    cz_system,
    snapshot,
    best_cz_point.q_flux,
    best_cz_point.tc_flux,
    cz_t_gate;
    branch = :cz,
    ramp_ns = FIG2_CZ_RAMP_NS,
    charge_cutoff = FIG2_CHARGE_CUTOFF,
)

(
    move_selected = compact_namedtuple(best_move_point),
    move_pulse = compact_namedtuple(move_pulse_summary),
    cz_selected = compact_namedtuple(best_cz_point),
    cz_pulse = compact_namedtuple(cz_pulse_summary),
    cz_computational_phase = compact_namedtuple(cz_phase_proxy),
)


if RUN_TRUNCATION_CHECKS
    move_cutoff4 = final_pulse_state_probabilities(
        move_system,
        snapshot,
        best_move_point.q_flux,
        best_move_point.tc_flux,
        (1, 0, 0),
        move_t_gate;
        branch = :move,
        ramp_ns = FIG2_MOVE_RAMP_NS,
        charge_cutoff = 4,
    )
    cz_cutoff4 = final_pulse_state_probabilities(
        cz_system,
        snapshot,
        best_cz_point.q_flux,
        best_cz_point.tc_flux,
        (1, 0, 1),
        cz_t_gate;
        branch = :cz,
        ramp_ns = FIG2_CZ_RAMP_NS,
        charge_cutoff = 4,
    )

    (
        move_cutoff3_final_p_e = round4(best_move_point.value),
        move_cutoff4_final_p_e = round4(move_cutoff4.final_p_e),
        cz_cutoff3_final_p_f = round4(best_cz_point.value),
        cz_cutoff4_final_p_f = round4(cz_cutoff4.final_p_f),
    )
else
    "Set RUN_TRUNCATION_CHECKS = true to rerun the pulse-schedule best points at charge_cutoff = 4."
end


final_summary = (
    move_branch = branch_spec(:move).label,
    cz_branch = branch_spec(:cz).label,
    move_gap = compact_namedtuple(move_gap),
    cz_gap = compact_namedtuple(cz_gap),
    best_move_point = compact_namedtuple(best_move_point),
    best_cz_point = compact_namedtuple(best_cz_point),
    move_pulse = compact_namedtuple(move_pulse_summary),
    cz_pulse = compact_namedtuple(cz_pulse_summary),
    cz_computational_phase = compact_namedtuple(cz_phase_proxy),
    exported_figures = figure_exports,
)

final_summary






