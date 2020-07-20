using Documenter, SDDP, Literate

# Run the farmer's problem first to precompile a bunch of SDDP.jl functions.
# This is a little sneaky, but it avoids leaking long (6 sec) compilation times
# into the examples.
include(joinpath(@__DIR__, "src", "examples", "the_farmers_problem.jl"))

"Call julia docs/make.jl --fix to rebuild the doctests."
const FIX_DOCTESTS = length(ARGS) == 1 && ARGS[1] == "--fix"

# doctest=:fix only works with `\n` line endings, so let's replace any `\r\n`
# ones.
function convert_line_endings()
    for file in readdir(joinpath(@__DIR__, "src", "tutorial"))
        !endswith(file, ".md") && continue
        filename = joinpath(@__DIR__, "src", "tutorial", file)
        code = read(filename, String)
        open(filename, "w") do io
            write(io, replace(code, "\r\n" => "\n"))
        end
    end
end
FIX_DOCTESTS && convert_line_endings()

for dir in ["examples", "tutorial"]
    for file in sort(readdir(joinpath(@__DIR__, "src", dir)))
        !endswith(file, ".jl") && continue
        filename = joinpath(@__DIR__, "src", dir, file)
        Literate.markdown(filename, dirname(filename); documenter=true)
    end
end

const EXAMPLES = Any[
    "examples/$(file)" for file  in filter(
        f -> endswith(f, ".md"),
        sort(readdir(joinpath(@__DIR__, "src", "examples")))
    )
]

const TUTORIAL_PAGES = Any[
    "tutorial/$(file)" for file  in filter(
        f -> endswith(f, ".md"),
        sort(readdir(joinpath(@__DIR__, "src", "tutorial")))
    )
]

const GUIDE_PAGES = Any[
    "guides/$(page)"
    for page in sort(readdir(joinpath(@__DIR__, "src", "guides")))
]

const ASSETS = [
    "logo.ico"
]

makedocs(
    sitename = "SDDP.jl",
    authors  = "Oscar Dowson",
    clean = true,
    doctest = FIX_DOCTESTS ? :fix : true,
    format = Documenter.HTML(
        assets = ASSETS,
        # See  https://github.com/JuliaDocs/Documenter.jl/issues/868
        prettyurls = get(ENV, "CI", nothing) == "true"
    ),
    strict = true,
    pages = [
        "Home" => "index.md",
        "Tutorials" => TUTORIAL_PAGES,
        "How-to guides" => GUIDE_PAGES,
        "Examples" => EXAMPLES,
        "API Reference" => "apireference.md"
    ],
    doctestfilters = [r"[\s\-]?\d\.\d{6}e[\+\-]\d{2}"]
)

deploydocs(repo = "github.com/odow/SDDP.jl.git")
