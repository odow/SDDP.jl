#  Copyright 2017-20, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public License,
#  v. 2.0. If a copy of the MPL was not distributed with this file, You can
#  obtain one at http://mozilla.org/MPL/2.0/.

const ASSET_DIR = dirname(@__FILE__)
const SPAGHETTI_HTML_FILE = joinpath(ASSET_DIR, "spaghetti_plot.html")
const SPAGHETTI_JS_FILE = joinpath(ASSET_DIR, "spaghetti_plot.js")
const D3_JS_FILE = joinpath(ASSET_DIR, "d3.v3.min.js")

const PLOT_DATA = Dict{String,Any}(
    "cumulative" => false,
    "title" => "",
    "ylabel" => "",
    "xlabel" => "Stages",
    "interpolate" => "linear",
    "ymin" => "",
    "ymax" => "",
)

"""
	SDDP.SpaghettiPlot(; stages, scenarios)

Initialize a new `SpaghettiPlot` with `stages` stages and `scenarios` number of
replications.
"""
struct SpaghettiPlot
    simulations::Vector{Vector{Dict{Symbol,Any}}}
    data::Vector{Dict{String,Any}}
    function SpaghettiPlot(simulations)
        return new(simulations, Dict{String,Any}[])
    end
end

function Base.show(io::IO, plt::SpaghettiPlot)
    print(
        io,
        "A spaghetti plot with ",
        length(plt.simulations),
        " scenarios ",
        "and ",
        length(plt.simulations[1]),
        " stages.",
    )
end

"""
	SDDP.add_spaghetti(data_function::Function, plt::SpaghettiPlot; kwargs...)

# Description

Add a new figure to the SpaghettiPlot `plt`, where the y-value of the
`scenario`th line when x = `stage` is given by
`data_function(plt.simulations[scenario][stage])`.

# Keyword arguments

 * `xlabel`: set the xaxis label
 * `ylabel`: set the yaxis label
 * `title`: set the title of the plot
 * `ymin`: set the minimum y value
 * `ymax`: set the maximum y value
 * `cumulative`: plot the additive accumulation of the value across the stages
 * `interpolate`: interpolation method for lines between stages.
 	Defaults to `"linear"` see [the d3 docs](https://github.com/d3/d3-3.x-api-reference/blob/master/SVG-Shapes.md#line_interpolate)
	for all options.

# Examples

	simulations = simulate(model, 10)
	plt = SDDP.spaghetti_plot(simulations)
	SDDP.add_spaghetti(plt; title = "Stage objective") do data
		return data[:stage_objective]
	end
"""
function add_spaghetti(
    data_function::Function,
    plt::SpaghettiPlot;
    xlabel = "Stages",
    ylabel = "",
    cumulative = false,
    title = "",
    interpolate = "linear",
    ymin = "",
    ymax = "",
)
    plot_dict = Dict{String,Any}(
        "xlabel" => xlabel,
        "ylabel" => ylabel,
        "title" => title,
        "cumulative" => cumulative,
        "interpolate" => interpolate,
        "ymin" => ymin,
        "ymax" => ymax,
    )
    plot_dict["data"] = Vector{Float64}[]
    for (i, scenario) in enumerate(plt.simulations)
        push!(plot_dict["data"], Float64[])
        series_value = 0.0
        for stage in scenario
            y_value = float(data_function(stage))
            if cumulative
                series_value += y_value
            else
                series_value = y_value
            end
            push!(plot_dict["data"][i], series_value)
        end
    end
    push!(plt.data, plot_dict)
    return
end

function fill_template(dest::String, args...; template::String, launch::Bool = false)
    s = read(template, String)
    for arg in args
        s = replace(s, arg)
    end
    write(dest, s)
    if launch
        launch_file(dest)
    end
    return
end

"""
	plot(plt::SpaghettiPlot[, filename::String]; open::Bool = true)

The SpaghettiPlot plot `plt` to `filename`. If `filename` is not given, it will
be saved to a temporary directory. If `open = true`, then a browser window will
be opened to display the resulting HTML file.
"""
function plot(
    plt::SpaghettiPlot,
    filename::String = joinpath(tempdir(), string(Random.randstring(), ".html"));
    open::Bool = true,
)
    fill_template(
        filename,
        "<!--DATA-->" => JSON.json(plt.data),
        "<!--D3.JS-->" => read(D3_JS_FILE, String),
        "<!--SPAGHETTI_PLOT.JS-->" => read(SPAGHETTI_JS_FILE, String);
        template = SPAGHETTI_HTML_FILE,
        launch = open,
    )
    return
end

function save(p::SpaghettiPlot, args...; kwargs...)
    Base.depwarn("`SDDP.save` is deprecated. Use `SDDP.plot` instead.", :save)
    return plot(p, args...; kwargs...)
end

function launch_file(filename)
    if Sys.iswindows()
        run(`$(ENV["COMSPEC"]) /c start $(filename)`)
    elseif Sys.isapple()
        run(`open $(filename)`)
    elseif Sys.islinux() || Sys.isbsd()
        run(`xdg-open $(filename)`)
    else
        error("Unable to show plot. Try opening the file $(filename) manually.")
    end
    return
end
