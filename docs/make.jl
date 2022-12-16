#  Copyright (c) 2017-22, Oscar Dowson and SDDP.jl contributors.
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

import Documenter
import Literate
import Random

"Call julia docs/make.jl --fix to rebuild the doctests."
const FIX_DOCTESTS = any(isequal("--fix"), ARGS)

const EXAMPLES_DIR = joinpath(@__DIR__, "src", "examples")
const TUTORIAL_BASIC_DIR = joinpath(@__DIR__, "src", "tutorial", "basic")
const TUTORIAL_ADVANCED_DIR = joinpath(@__DIR__, "src", "tutorial", "advanced")
const TUTORIAL_THEORY_DIR = joinpath(@__DIR__, "src", "tutorial", "theory")
const GUIDES_DIR = joinpath(@__DIR__, "src", "guides")

_sorted_files(dir, ext) = sort(filter(f -> endswith(f, ext), readdir(dir)))

function list_of_sorted_files(prefix, dir, ext = ".md")
    return Any["$(prefix)/$(file)" for file in _sorted_files(dir, ext)]
end

# Run the farmer's problem first to precompile a bunch of SDDP.jl functions.
# This is a little sneaky, but it avoids leaking long (6 sec) compilation times
# into the examples.
include(joinpath(EXAMPLES_DIR, "the_farmers_problem.jl"))

if FIX_DOCTESTS
    # doctest=:fix only works with `\n` line endings. Replace any `\r\n` ones.
    for dir in [TUTORIAL_BASIC_DIR, TUTORIAL_ADVANCED_DIR, TUTORIAL_THEORY_DIR]
        for filename in list_of_sorted_files(dir, dir)
            code = read(filename, String)
            write(filename, replace(code, "\r\n" => "\n"))
        end
    end
end

for dir in [
    EXAMPLES_DIR,
    TUTORIAL_BASIC_DIR,
    TUTORIAL_ADVANCED_DIR,
    TUTORIAL_THEORY_DIR,
]
    for jl_filename in list_of_sorted_files(dir, dir, ".jl")
        Random.seed!(12345)
        Literate.markdown(jl_filename, dir; documenter = true, execute = true)
        md_filename = jl_filename[1:(end-3)] * ".md"
        md = read(md_filename, String)
        write(md_filename, replace(md, "nothing #hide" => ""))
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
    # @XXX -> [@XXX](https://github.com/XXX)
    while (m = match(r"@([a-zA-Z0-9\-]+)", line)) !== nothing
        tag = m.captures[1]
        line = replace(line, m.match => "[@$tag](https://github.com/$tag")
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

Documenter.makedocs(
    sitename = "SDDP.jl",
    authors = "Oscar Dowson",
    format = Documenter.HTML(
        # See  https://github.com/JuliaDocs/Documenter.jl/issues/868
        prettyurls = get(ENV, "CI", nothing) == "true",
        collapselevel = 1,
    ),
    clean = true,
    doctest = FIX_DOCTESTS ? :fix : true,
    strict = true,
    pages = [
        "Home" => "index.md",
        "Tutorials" => Any[
            "Basic"=>list_of_sorted_files("tutorial/basic", TUTORIAL_BASIC_DIR),
            "Advanced"=>list_of_sorted_files(
                "tutorial/advanced",
                TUTORIAL_ADVANCED_DIR,
            ),
            "Theory"=>list_of_sorted_files(
                "tutorial/theory",
                TUTORIAL_THEORY_DIR,
            ),
        ],
        "How-to guides" => list_of_sorted_files("guides", GUIDES_DIR),
        "Examples" => list_of_sorted_files("examples", EXAMPLES_DIR),
        "API Reference" => "apireference.md",
        "Release notes" => "release_notes.md",
    ],
    doctestfilters = [r"[\s\-]?\d\.\d{6}e[\+\-]\d{2}"],
)

Documenter.deploydocs(repo = "github.com/odow/SDDP.jl.git", push_preview = true)
