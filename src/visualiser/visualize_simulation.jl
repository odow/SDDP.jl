#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

const SIMULATION_ASSET_DIR = dirname(@__FILE__)
const SIMULATION_HTML_FILE = joinpath(SIMULATION_ASSET_DIR, "simulation.template.html")
const SIMULATION_ASSETS    = ["d3.v3.min.js", "simulation.js", "simulation.css"]

immutable SimulationPlot
	data::Vector{Dict{String, Any}}
end
"""
	SDDP.newplot()

# Description

Initialize a new `SimulationPlot`.
"""
newplot() = SimulationPlot(Dict{String, Any}[])

JSON.json(p::SimulationPlot) = json(p.data)

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
function addplot!(p::SimulationPlot, i::AbstractVector{Int}, t::AbstractVector{Int}, f::Function; xlabel="", ylabel="",
	cumulative=false, title="", interpolate="linear", ymin="", ymax="")

	plot_dict = deepcopy(PLOT_DATA)
	plot_dict["data"] = Vector{Float64}[]
	for ii in i
		push!(plot_dict["data"], Float64[])
		y = 0.0
		for tt in t
			if cumulative
				y += f(ii, tt)
			else
				y = f(ii, tt)
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
	html_string = readstring(SIMULATION_HTML_FILE)
	html_string = replace(html_string, "<!--DATA-->", json(p))
	launch_file(html_string, SIMULATION_ASSETS, SIMULATION_ASSET_DIR)
end

function launch_file(html_string, assets, asset_dir)
	temporary_html_file = replace(tempname(), ".tmp", ".html")
    for asset in assets
        cp(joinpath(asset_dir, asset), joinpath(dirname(temporary_html_file), asset), remove_destination=true)
    end
    open(temporary_html_file, "w") do f
        write(f, html_string)
    end
    launch_plot(temporary_html_file)
end
