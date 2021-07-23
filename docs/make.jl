using Documenter, ACME

makedocs(
    modules = [ACME],
    sitename = "ACME.jl",
    pages = Any[
        "Home" => "index.md",
        "gettingstarted.md",
        "ug.md",
        "elements.md",
    ],
)

deploydocs(
    repo = "github.com/HSU-ANT/ACME.jl.git",
    devbranch = "main",
)
