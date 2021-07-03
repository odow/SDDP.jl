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
            "Basic" =>
                list_of_sorted_files("tutorial/basic", TUTORIAL_BASIC_DIR),
            "Advanced" =>
                list_of_sorted_files("tutorial/advanced", TUTORIAL_ADVANCED_DIR),
            "Theory" =>
                list_of_sorted_files("tutorial/theory", TUTORIAL_THEORY_DIR),
        ],
        "How-to guides" => list_of_sorted_files("guides", GUIDES_DIR),
        "Examples" => list_of_sorted_files("examples", EXAMPLES_DIR),
        "API Reference" => "apireference.md",
    ],
    doctestfilters = [r"[\s\-]?\d\.\d{6}e[\+\-]\d{2}"],
)

Documenter.deploydocs(repo = "github.com/odow/SDDP.jl.git", push_preview = true)
