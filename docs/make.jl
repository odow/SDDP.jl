using Documenter, SDDP

makedocs(
    format = :html,
    sitename = "SDDP.jl",
    pages = [
        "Introduction" => "index.md",
        "Quick Start" => "quick.md"
    ]
)
