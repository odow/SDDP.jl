#  Copyright 2017, Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

# functions shared by all visualizations

const ASSET_DIR = dirname(@__FILE__)

function launch_file(html_string, assets, outputfile=replace(tempname(), ".tmp", ".html"))
    for asset in assets
        cp(joinpath(ASSET_DIR, asset), joinpath(dirname(outputfile), asset), remove_destination=true)
    end
    open(outputfile, "w") do f
        write(f, html_string)
    end
	if is_windows()
		run(`$(ENV["COMSPEC"]) /c start $(outputfile)`)
	elseif is_apple()
        run(`open $(outputfile)`)
    elseif is_linux() || is_bsd()
        run(`xdg-open $(outputfile)`)
    end
end

function prephtml(html_template, args...)
	html = readstring(joinpath(ASSET_DIR, html_template))
	for (key, val) in args
		html = replace(html, key, val)
	end
	html
end

include("visualize_simulation.jl")
include("visualize_value_function.jl")
