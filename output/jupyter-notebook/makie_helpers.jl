using CairoMakie

const NOTEBOOK_PANEL = (920, 540)
const NOTEBOOK_WIDE = (1080, 620)
const NOTEBOOK_TALL = (920, 760)

function activate_notebook_theme!()
    CairoMakie.activate!(type = "svg")
    set_theme!(
        Theme(
            fontsize = 18,
            figure_padding = 12,
            Axis = (
                backgroundcolor = :white,
                topspinevisible = false,
                rightspinevisible = false,
                xgridvisible = false,
                ygridvisible = false,
                xlabelsize = 20,
                ylabelsize = 20,
                titlesize = 22,
                xticklabelsize = 16,
                yticklabelsize = 16,
            ),
            Legend = (
                framevisible = false,
                patchsize = (28, 14),
                labelsize = 16,
            ),
            Colorbar = (
                ticklabelsize = 14,
                labelsize = 16,
            ),
            Lines = (
                linewidth = 3,
            ),
        ),
    )
    return nothing
end

function line_figure(
    curves;
    title::AbstractString,
    xlabel::AbstractString,
    ylabel::AbstractString,
    size = NOTEBOOK_PANEL,
    legend_position = :rt,
    xlims = nothing,
    ylims = nothing,
)
    fig = Figure(size = size)
    ax = Axis(fig[1, 1]; title, xlabel, ylabel)
    has_labeled_curve = false

    for curve in curves
        label = get(curve, :label, nothing)
        has_labeled_curve |= !isnothing(label)
        if isnothing(label)
            lines!(
                ax,
                curve.x,
                curve.y;
                color = get(curve, :color, :steelblue3),
                linestyle = get(curve, :linestyle, :solid),
                linewidth = get(curve, :linewidth, 3),
            )
        else
            lines!(
                ax,
                curve.x,
                curve.y;
                label,
                color = get(curve, :color, :steelblue3),
                linestyle = get(curve, :linestyle, :solid),
                linewidth = get(curve, :linewidth, 3),
            )
        end
    end

    isnothing(xlims) || xlims!(ax, xlims...)
    isnothing(ylims) || ylims!(ax, ylims...)
    has_labeled_curve && axislegend(ax; position = legend_position)
    return fig
end

function heatmap_figure(
    x_values,
    y_values,
    matrix;
    title::AbstractString,
    xlabel::AbstractString,
    ylabel::AbstractString,
    colorlabel::AbstractString,
    size = NOTEBOOK_WIDE,
    colormap = :viridis,
    colorrange = nothing,
    xlims = nothing,
    ylims = nothing,
)
    fig = Figure(size = size)
    ax = Axis(fig[1, 1]; title, xlabel, ylabel)
    xspan = (minimum(x_values), maximum(x_values))
    yspan = (minimum(y_values), maximum(y_values))
    hm = isnothing(colorrange) ?
        image!(ax, xspan, yspan, matrix; colormap, interpolate = false) :
        image!(ax, xspan, yspan, matrix; colormap, colorrange, interpolate = false)
    isnothing(xlims) || xlims!(ax, xlims...)
    isnothing(ylims) || ylims!(ax, ylims...)
    Colorbar(fig[1, 2], hm; label = colorlabel)
    return fig
end

function save_figure(fig, project_root::AbstractString, stem::AbstractString; subdir::AbstractString = "renger2026")
    output_dir = joinpath(project_root, "output", "figures", subdir)
    mkpath(output_dir)
    svg_path = joinpath(output_dir, stem * ".svg")
    png_path = joinpath(output_dir, stem * ".png")
    save(svg_path, fig)
    save(png_path, fig)
    return (svg = svg_path, png = png_path)
end
