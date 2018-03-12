using Documenter, SDDP

makedocs(
    clean = false,
    format = :html,
    sitename = "SDDP.jl",
    pages = [
        "Manual" => "index.md",
        "Examples" => "examples.md",
        "Readings" => "readings.md",
        "Reference" => "apireference.md"
    ],
    assets = ["assets/custom.css", "3d.gif"]
)

deploydocs(
    repo   = "github.com/odow/SDDP.jl.git",
    target = "build",
    osname = "linux",
    julia  = "0.6",
    deps   = nothing,
    make   = nothing,
)
