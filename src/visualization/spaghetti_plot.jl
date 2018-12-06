#  Copyright 2018, Oscar Dowson.
#  This Source Code Form is subject to the terms of the Mozilla Public License,
#  v. 2.0. If a copy of the MPL was not distributed with this file, You can
#  obtain one at http://mozilla.org/MPL/2.0/.

const ASSET_DIR = dirname(@__FILE__)
const SIMULATION_HTML_FILE = "spaghetti_plot.html"
const SIMULATION_ASSETS    = ["d3.v3.min.js", "spaghetti_plot.js"]
const PLOT_DATA = Dict{String, Any}(
	"cumulative"  => false,
	"title"       => "",
	"ylabel"      => "",
	"xlabel"      => "Stages",
	"interpolate" => "linear",
	"ymin"        => "",
	"ymax"        => ""
)

struct SpaghettiPlot
	stages::Int
	scenarios::Int
	data::Vector{Dict{String, Any}}
end

"""
	Kokako.spaghetti_plot()

# Description

Initialize a new `SpaghettiPlot`.
"""
function spaghetti_plot(; stages, scenarios)
	return SpaghettiPlot(stages, scenarios, Dict{String, Any}[])
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
	show(plt::SpaghettiPlot)

Launch a browser and render the SpaghettiPlot plot `p`.
"""
function Base.show(plt::SpaghettiPlot)
	Base.show(joinpath(tempdir(), string(Random.randstring(), ".html")), plt)
end

"""
	show(filename::String, plt::SpaghettiPlot)

Launch a browser and render the SpaghettiPlot plot `plt`. Save the resulting
HTML file to `filename`.
"""
function Base.show(filename::String, plt::SpaghettiPlot)
	launch_file(prep_html(plt), SIMULATION_ASSETS, filename)
end

function prep_html(plt::SpaghettiPlot)
	html = read(joinpath(ASSET_DIR, SIMULATION_HTML_FILE), String)
	html = replace(html, "<!--DATA-->" => JSON.json(plt.data))
	return html
end

function launch_file(html_string, assets, output_file)
    for asset in assets
        cp(joinpath(ASSET_DIR, asset), joinpath(dirname(output_file), asset),
			force = true)
    end
    open(output_file, "w") do io
        write(io, html_string)
    end
	if Sys.iswindows()
		run(`$(ENV["COMSPEC"]) /c start $(output_file)`)
	elseif Sys.isapple()
        run(`open $(output_file)`)
    elseif Sys.islinux() || Sys.isbsd()
        run(`xdg-open $(output_file)`)
	else
		error("Unable to show spaghetti plot. Try opening the file " *
		      "$(output_file) manually.")
    end
	return
end
