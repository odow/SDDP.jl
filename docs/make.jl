import Documenter
import Literate
import Random

"Call julia docs/make.jl --fix to rebuild the doctests."
const FIX_DOCTESTS = any(isequal("--fix"), ARGS)

const EXAMPLES_DIR = joinpath(@__DIR__, "src", "examples")
const TUTORIAL_DIR = joinpath(@__DIR__, "src", "tutorial")
const GUIDES_DIR = joinpath(@__DIR__, "src", "guides")

sorted_files(dir, ext) = sort(filter(f -> endswith(f, ext), readdir(dir)))

# Run the farmer's problem first to precompile a bunch of SDDP.jl functions.
# This is a little sneaky, but it avoids leaking long (6 sec) compilation times
# into the examples.
include(joinpath(EXAMPLES_DIR, "the_farmers_problem.jl"))


if FIX_DOCTESTS
    # doctest=:fix only works with `\n` line endings. Replace any `\r\n` ones.
    for file in sorted_files(TUTORIAL_DIR, ".md")
        filename = joinpath(TUTORIAL_DIR, file)
        code = read(filename, String)
        write(filename, replace(code, "\r\n" => "\n"))
    end
end

for dir in [EXAMPLES_DIR, TUTORIAL_DIR]
    for file in sorted_files(dir, ".jl")
        Random.seed!(12345)
        jl_filename = joinpath(dir, file)
        Literate.markdown(
            jl_filename,
            dir;
            documenter = true,
            execute = true,
        )
        md_filename = jl_filename[1:(end - 3)] * ".md"
        md = read(md_filename, String)
        write(md_filename, replace(md, "nothing #hide" => ""))
    end
end

Documenter.makedocs(
    sitename = "SDDP.jl",
    authors  = "Oscar Dowson",
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
            "tutorial/$(file)" for file in sorted_files(TUTORIAL_DIR, ".md")
        ],
        "How-to guides" => Any[
            "guides/$(file)" for file in sorted_files(GUIDES_DIR, ".md")
        ],
        "Examples" => Any[
            "examples/$(file)" for file in sorted_files(EXAMPLES_DIR, ".md")
        ],
        "API Reference" => "apireference.md",
    ],
    doctestfilters = [r"[\s\-]?\d\.\d{6}e[\+\-]\d{2}"],
)

Documenter.deploydocs(
    repo = "github.com/odow/SDDP.jl.git",
    push_preview = true,
)
