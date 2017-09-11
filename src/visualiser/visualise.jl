#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

const ASSET_DIR = dirname(@__FILE__)
const HTML_FILE = joinpath(ASSET_DIR, "simulation.template.html")
const ASSETS    = ["d3.v3.min.js", "simulation.js", "simulation.css"]
const PLOT_DATA = Dict{String, Any}(
	"cumulative"  => false,
	"title"       => "",
	"ylabel"      => "",
	"xlabel"      => "Stages",
	"interpolate" => "linear",
	"ymin"        => "",
	"ymax"        => ""
)


function add_to_output!(output::Dict{String, Any}, sym, value)
	if haskey(output, string(sym))
		output[string(sym)] = value
	else
		error("Keyword $(sym)=$value not recognised in @visualise.")
	end
end

function adddata!(plot_dict::Dict{String, Any}, results::Vector, func::Function)
	plot_dict["data"] = Vector{Float64}[]
	for replication=1:length(results)
		push!(plot_dict["data"], Float64[])
		y = 0.0
		for stage=1:length(results[replication][:stageobjective])
			if plot_dict["cumulative"]
				y += func(results, replication, stage)
			else
				y = func(results, replication, stage)
			end
			push!(plot_dict["data"][replication], y)
		end
	end
end
adddata!(plot_dict::Dict{String, Any}, results::Dict, func::Function) = adddata!(plot_dict, [results], func)

function launch_plot(html_file)
	if is_windows()
		run(`$(ENV["COMSPEC"]) /c start $(html_file)`)
	elseif is_apple()
        run(`open $(html_file)`)
    elseif is_linux() || is_bsd()
        run(`xdg-open $(html_file)`)
    end
end
"""
	@visualise(results, replication, stage, begin
		... plot definitions ...
	end)

# Description

Plot everything using interactive javascript. This will launch an HTML page
to explore.

# Usage

    @visualise(results, i, t, begin
        ... one line for each plot ...
    end)

where results is the vector of result dictionaries from `simulate()`, `i` is the
simulation index (`1:length(results)`), and `t` is the stage index (`1:T`).

Each plot line gets transformed into an anonymous function

	(results, i, t) -> ... plot line ...

so can be any valid Julia syntax that uses `results`, `i`, or `t` as an argument.

After the plot definition, keyword arguments can be used (in parenthesises):
 * `title`       - set the title of the plot
 * `ylabel`      - set the yaxis label
 * `xlabel`      - set the xaxis label
 * `interpolate` - interpolate lines between stages. Defaults to `"linear"`
    see [the d3 docs](https://github.com/d3/d3-3.x-api-reference/blob/master/SVG-Shapes.md#line_interpolate)
	for all options.

# Results Object

`results::Vector{Dict{Symbol, Any}}` is a vector of dictionaries where each
dictionary corresponds to one simulation (therefore there will be
`N = length(results)` lines plotted in each graph).

# Examples
	@visualise(results, i, t, begin
		results[i][:stageobjective][t], (title="Accumulated Obj.", ylabel="\$", cumulative=true)
		results[i][:stageobjective][t], (title="Stage Objective",  ylabel="\$")
		results[i][:x][t],              (title="State",            ylabel="Level")
		results[i][:u][t],              (title="Control")
		prices[t, results[i][:markov][t]], (ylabel="Price", interpolate="step-after")
	end)
"""
macro visualise(results, replication, stage, block)
	@assert block.head == :block || error("Invalid syntax for @visualise")
	kw = Expr(:tuple,)
	push!(kw.args, results)
	push!(kw.args, replication)
	push!(kw.args, stage)
    code = quote
		plot_list = Dict{String, Any}[]
	end
    for it in block.args
        Base.Meta.isexpr(it, :line) && continue
		plot_dict = deepcopy(PLOT_DATA)
        if it.head == :tuple
            if length(it.args) == 2
				if it.args[2].head == :tuple
					for arg in it.args[2].args
						if arg.head != :(=)
							error("Must be a keyword argument in @visualise: $(arg)")
						end
						add_to_output!(plot_dict, arg.args[1], arg.args[2])
					end
				elseif it.args[2].head == :(=)
					add_to_output!(plot_dict, it.args[2].args[1],it.args[2].args[2])
				end
				f = Expr(:->, kw, Expr(:block, it.args[1]))
			elseif length(it.args) == 1
				f = Expr(:->, kw, Expr(:block, it.args))
			else
                error("Unknown arguments in @visualise")
			end
        else
			f = Expr(:->, kw, Expr(:block, it))
        end
		push!(code.args, quote
			adddata!($plot_dict, $(esc(results)), $(esc(f)))
			push!(plot_list, $plot_dict)
		end)
    end

	push!(code.args, quote
		temporary_html_file = replace(tempname(), ".tmp", ".html")
		html_string = readstring($HTML_FILE)
		for asset in $ASSETS
			cp(joinpath(ASSET_DIR, asset), joinpath(dirname(temporary_html_file), asset), remove_destination=true)
		end
		json_data = json(plot_list)
		open(temporary_html_file, "w") do f
			write(f, replace(html_string, "<!--DATA-->", json_data))
		end
		launch_plot(temporary_html_file)
	end)
    return code
end
