using Documenter, ACME

makedocs(
    modules = [ACME],
    format = :html,
    sitename = "ACME.jl",
    pages = Any[
        "Home" => "index.md",
        "gettingstarted.md",
        "ug.md",
        "elements.md",
    ],
)

deploydocs(
    target = "build",
    repo = "github.com/HSU-ANT/ACME.jl.git",
    latest = "develop",
    deps = nothing,
    make = nothing,
    julia = "1.0",
)
