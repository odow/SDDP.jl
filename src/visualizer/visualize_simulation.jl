#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

const SIMULATION_HTML_FILE = "simulation.template.html"
const SIMULATION_ASSETS    = ["d3.v3.min.js", "simulation.js", "simulation.css"]
const PLOT_DATA = Dict{String, Any}(
	"cumulative"  => false,
	"title"       => "",
	"ylabel"      => "",
	"xlabel"      => "Stages",
	"interpolate" => "linear",
	"ymin"        => "",
	"ymax"        => ""
)

immutable SimulationPlot
	data::Vector{Dict{String, Any}}
end
"""
	SDDP.newplot()

# Description

Initialize a new `SimulationPlot`.
"""
newplot() = SimulationPlot(Dict{String, Any}[])

"""
	SDDP.addplot!(p::SimulationPlot, ivals::AbstractVector{Int}, tvals::AbstractVector{Int}, f::Function; kwargs...)

# Description

Add a new figure to the SimulationPlot `p`, where the y-value is given by `f(i, t)`
for all `i` in `ivals` (one for each series) and `t` in `tvals` (one for each stage).

# Keywords
 * `xlabel`: set the xaxis label;
 * `ylabel`: set the yaxis label;
 * `title`: set the title of the plot;
 * `ymin`: set the minimum y value;
 * `ymax`: set the maximum y value;
 * `cumulative`: plot the additive accumulation of the value across the stages;
 * `interpolate`: interpolation method for lines between stages. Defaults to `"linear"`
    see [the d3 docs](https://github.com/d3/d3-3.x-api-reference/blob/master/SVG-Shapes.md#line_interpolate)
	for all options.

# Examples

	results = simulate(m, 10)
	p = SDDP.newplot()
	SDDP.addplot!(p, 1:10, 1:3, (i,t)->results[i][:stageobjective][t])
"""
function addplot!(p::SimulationPlot, i::AbstractVector{Int}, t::AbstractVector{Int}, f::Function; xlabel="Stages", ylabel="",
	cumulative=false, title="", interpolate="linear", ymin="", ymax="")

	plot_dict = Dict{String, Any}(
		"xlabel" => xlabel,
		"ylabel" => ylabel,
		"title" => title,
		"cumulative" => cumulative,
		"interpolate" => interpolate ,
		"ymin" => ymin,
		"ymax" => ymax

	)
	plot_dict["data"] = Vector{Float64}[]
	for ii in i
		push!(plot_dict["data"], Float64[])
		y = 0.0
		for tt in t
			if cumulative
				y += float(f(ii, tt))
			else
				y = float(f(ii, tt))
			end
			push!(plot_dict["data"][ii], y)
		end
	end
	push!(p.data, plot_dict)
end
"""
	show(p::SimulationPlot)

# Description

Launch a browser and render the SimulationPlot plot `p`.
"""
function show(p::SimulationPlot)
	html = prephtml(p)
	launch_file(html, SIMULATION_ASSETS)
end

function show(filename::String, p::SimulationPlot)
	html = prephtml(p)
	launch_file(html, SIMULATION_ASSETS, filename)
end

prephtml(p::SimulationPlot) = prephtml(SIMULATION_HTML_FILE, ("<!--DATA-->", json(p.data)))
