#  Copyright (c) 2017-25, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

import Documenter
import Literate
import Random
import Test

"Call julia docs/make.jl --fix to rebuild the doctests."
const FIX_DOCTESTS = any(isequal("--fix"), ARGS)

const EXAMPLES_DIR = joinpath(@__DIR__, "src", "examples")

_sorted_files(dir, ext) = sort(filter(f -> endswith(f, ext), readdir(dir)))

function list_of_sorted_files(prefix, dir, ext = ".md")
    return Any["$(prefix)/$(file)" for file in _sorted_files(dir, ext)]
end

function _include_sandbox(filename)
    mod = @eval module $(gensym()) end
    return Base.include(mod, filename)
end

# Run the farmer's problem first to precompile a bunch of SDDP.jl functions.
# This is a little sneaky, but it avoids leaking long (6 sec) compilation times
# into the examples.
include(joinpath(EXAMPLES_DIR, "the_farmers_problem.jl"))

if FIX_DOCTESTS
    # doctest=:fix only works with `\n` line endings. Replace any `\r\n` ones.
    for dir in joinpath.(@__DIR__, "src", ("tutorial", "explanation"))
        for filename in list_of_sorted_files(dir, dir)
            code = read(filename, String)
            write(filename, replace(code, "\r\n" => "\n"))
        end
    end
end

function add_binder_links(filename, content)
    filename = replace(filename, ".jl" => ".ipynb")
    links = """


    #md # [![](https://mybinder.org/badge_logo.svg)](@__BINDER_ROOT_URL__/$filename)
    #md # [![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](@__NBVIEWER_ROOT_URL__/$filename)
    """
    m = match(r"(\# \# .+)", content)
    return replace(content, m[1] => m[1] * links)
end

function _link_example(content, filename)
    title_line = findfirst(r"\n# .+?\n", content)
    line = content[title_line]
    ipynb = filename[1:end-3] * ".ipynb"
    new_title = string(
        line,
        "\n",
        "_This tutorial was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl)._\n",
        "[_Download the source as a `.jl` file_]($filename).\n",
        "[_Download the source as a `.ipynb` file_]($ipynb).\n",
    )
    contennt = replace(content, "nothing #hide" => ""),
    return replace(content, line => new_title)
end

for dir in joinpath.(@__DIR__, "src", ("examples", "tutorial", "explanation"))
    for jl_filename in list_of_sorted_files(dir, dir, ".jl")
        Random.seed!(12345)
        # `include` the file to test it before `#src` lines are removed. It is
        # in a testset to isolate local variables between files.
        Test.@testset "$jl_filename" begin
            _include_sandbox(jl_filename)
        end
        Random.seed!(12345)
        filename = replace(jl_filename, dirname(jl_filename) * "/" => "")
        Literate.markdown(
            jl_filename,
            dir;
            documenter = true,
            postprocess = content -> _link_example(content, filename),
            # Turn off the footer. We manually add a modified one.
            credit = false,
        )
        Literate.notebook(jl_filename, dir; execute = false, credit = false)
    end
end

# ==============================================================================
#  Modify the release notes
# ==============================================================================

function fix_release_line(
    line::String,
    url::String = "https://github.com/odow/SDDP.jl",
)
    # (#XXXX) -> ([#XXXX](url/issue/XXXX))
    while (m = match(r"\(\#([0-9]+)\)", line)) !== nothing
        id = m.captures[1]
        line = replace(line, m.match => "([#$id]($url/issues/$id))")
    end
    # ## vX.Y.Z -> [vX.Y.Z](url/releases/tag/vX.Y.Z)
    while (m = match(r"\#\# (v[0-9]+.[0-9]+.[0-9]+)", line)) !== nothing
        tag = m.captures[1]
        line = replace(line, m.match => "## [$tag]($url/releases/tag/$tag)")
    end
    # (Thanks @XXX) -> (Thanks [@XXX](https://github.com/XXX))
    while (m = match(r"\(Thanks \@(.+)\)", line)) !== nothing
        tag = m.captures[1]
        line = replace(
            line,
            m.match => "(Thanks [@$tag](https://github.com/$tag))",
        )
    end
    return line
end

open(joinpath(@__DIR__, "src", "changelog.md"), "r") do in_io
    open(joinpath(@__DIR__, "src", "release_notes.md"), "w") do out_io
        for line in readlines(in_io; keep = true)
            write(out_io, fix_release_line(line))
        end
    end
end

Documenter.makedocs(;
    sitename = "SDDP.jl",
    authors = "Oscar Dowson",
    format = Documenter.HTML(;
        analytics = "G-HZQQDVMPZW",
        # See  https://github.com/JuliaDocs/Documenter.jl/issues/868
        prettyurls = get(ENV, "CI", nothing) == "true",
        collapselevel = 1,
        sidebar_sitename = false,
        size_threshold_ignore = [
            "apireference.md",
            "examples/objective_state_newsvendor.md",
        ],
    ),
    clean = true,
    doctest = FIX_DOCTESTS ? :fix : true,
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "tutorial/first_steps.md",
            "tutorial/objective_uncertainty.md",
            "tutorial/markov_uncertainty.md",
            "tutorial/plotting.md",
            "tutorial/warnings.md",
            "tutorial/arma.md",
            "tutorial/decision_hazard.md",
            "tutorial/objective_states.md",
            "tutorial/pglib_opf.md",
            "tutorial/mdps.md",
            "tutorial/example_newsvendor.md",
            "tutorial/example_reservoir.md",
            "tutorial/example_milk_producer.md",
            "tutorial/inventory.md",
        ],
        "How-to guides" => [
            "guides/access_previous_variables.md",
            "guides/add_a_multidimensional_state_variable.md",
            "guides/add_a_risk_measure.md",
            "guides/add_integrality.md",
            "guides/add_multidimensional_noise.md",
            "guides/add_noise_in_the_constraint_matrix.md",
            "guides/choose_a_stopping_rule.md",
            "guides/create_a_general_policy_graph.md",
            "guides/debug_a_model.md",
            "guides/improve_computational_performance.md",
            "guides/simulate_using_a_different_sampling_scheme.md",
            "guides/create_a_belief_state.md",
        ],
        "Explanation" => ["explanation/theory_intro.md", "explanation/risk.md"],
        "Examples" => list_of_sorted_files("examples", EXAMPLES_DIR),
        "API Reference" => "apireference.md",
        "Release notes" => "release_notes.md",
    ],
    doctestfilters = [r"[\s\-]?\d\.\d{6}e[\+\-]\d{2}"],
)

Documenter.deploydocs(;
    repo = "github.com/odow/SDDP.jl.git",
    push_preview = true,
)
