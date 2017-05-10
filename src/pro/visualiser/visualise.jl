#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

using JSON

export @visualise

const ASSET_DIR = dirname(@__FILE__)
const HTML_FILE = joinpath(ASSET_DIR, "visualise.html")
const ASSETS    = ["d3.v3.min.js", "visualise.js", "visualise.css"]
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

function adddata!(plot_dict::Dict{String, Any}, results::Vector{Dict{Symbol, Any}}, func::Function)
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
adddata!(plot_dict::Dict{String, Any}, results::Dict{Symbol, Any}, func::Function) = adddata!(plot_dict, [results], func)

function launch_plot(html_file)
	if is_windows()
		run(`$(ENV["COMSPEC"]) /c start $(html_file)`)
	elseif is_apple()
        run(`open $(html_file)`)
    elseif is_linux() || is_bsd()
        run(`xdg-open $(html_file)`)
    end
end

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
		open(temporary_html_file, "w") do f
			write(f, replace(html_string, "<!--DATA-->", json(plot_list)))
		end
		launch_plot(temporary_html_file)
	end)
    return code
end
