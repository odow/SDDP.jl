#  Copyright 2018, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public License,
#  v. 2.0. If a copy of the MPL was not distributed with this file, You can
#  obtain one at http://mozilla.org/MPL/2.0/.

const ASSET_DIR = dirname(@__FILE__)
const SPAGHETTI_HTML_FILE = joinpath(ASSET_DIR, "spaghetti_plot.html")
const SPAGHETTI_JS_FILE = joinpath(ASSET_DIR, "spaghetti_plot.js")
const D3_JS_FILE = joinpath(ASSET_DIR, "d3.v3.min.js")

const PLOT_DATA = Dict{String, Any}(
	"cumulative"  => false,
	"title"       => "",
	"ylabel"      => "",
	"xlabel"      => "Stages",
	"interpolate" => "linear",
	"ymin"        => "",
	"ymax"        => ""
)


"""
	Kokako.SpaghettiPlot(; stages, scenarios)

Initialize a new `SpaghettiPlot` with `stages` stages and `scenarios` number of
replications.
"""
struct SpaghettiPlot
	stages::Int
	scenarios::Int
	data::Vector{Dict{String, Any}}
	function SpaghettiPlot(; stages, scenarios)
		return new(stages, scenarios, Dict{String, Any}[])
	end
end

"""
	Kokako.add_spaghetti(data_function::Function, plt::SpaghettiPlot; kwargs...)

# Description

Add a new figure to the SpaghettiPlot `plt`, where the y-value is given by
`data_function(stage, scenario)` for all `stage`s in `1:plt.stages` and
`scenario`s in `1:plt.scenarios`.

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
	plt = Kokako.spaghetti_plot(stages = 10, scenarios = 10)
	Kokako.add_spaghetti(plt; title = "Stage objective") do scenario, stage
		return simulations[scenario][stage][:stageobjective]
	end
"""
function add_spaghetti(data_function::Function, plt::SpaghettiPlot;
		xlabel = "Stages", ylabel = "", cumulative = false, title = "",
		interpolate = "linear", ymin = "", ymax = "")
	plot_dict = Dict{String, Any}(
		"xlabel" => xlabel,
		"ylabel" => ylabel,
		"title" => title,
		"cumulative" => cumulative,
		"interpolate" => interpolate,
		"ymin" => ymin,
		"ymax" => ymax
	)
	plot_dict["data"] = Vector{Float64}[]
	for scenario in 1:plt.scenarios
		push!(plot_dict["data"], Float64[])
		series_value = 0.0
		for stage in 1:plt.stages
			if cumulative
				series_value += float(data_function(scenario, stage))
			else
				series_value = float(data_function(scenario, stage))
			end
			push!(plot_dict["data"][scenario], series_value)
		end
	end
	push!(plt.data, plot_dict)
	return
end

"""
	show(plt::SpaghettiPlot[, filename::String])

Launch a browser and render the SpaghettiPlot plot `plt`. If `filename` is
given, save the resulting HTML file to `filename`.
"""
function Base.show(plt::SpaghettiPlot, filename::String =
		joinpath(tempdir(), string(Random.randstring(), ".html")))
	prep_html(plt, filename)
	return launch_file(filename)
end

function prep_html(plt::SpaghettiPlot, filename::String)
	html_string = read(SPAGHETTI_HTML_FILE, String)
	for pair in ["<!--DATA-->" => JSON.json(plt.data),
				 "<!--D3.JS-->" => read(D3_JS_FILE, String),
				 "<!--SPAGHETTI_PLOT.JS-->" => read(SPAGHETTI_JS_FILE, String)]
		html_string = replace(html_string, pair)
	end
	open(filename, "w") do io
        write(io, html_string)
    end
	return
end

function launch_file(filename)
    if Sys.iswindows()
		run(`$(ENV["COMSPEC"]) /c start $(filename)`)
	elseif Sys.isapple()
        run(`open $(filename)`)
    elseif Sys.islinux() || Sys.isbsd()
        run(`xdg-open $(filename)`)
	else
		error("Unable to show spaghetti plot. Try opening the file " *
		      "$(filename) manually.")
    end
	return
end
