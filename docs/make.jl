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

if v"1.0" ≤ VERSION < v"1.1" # deploy from 1.0
    deploydocs(
        repo = "github.com/HSU-ANT/ACME.jl.git",
        devbranch = "develop",
    )
end
