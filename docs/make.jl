using Documenter, ACME

makedocs(
    modules = [ACME],
    format = Documenter.Formats.HTML,
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
)
