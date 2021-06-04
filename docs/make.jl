using Documenter, ACME

makedocs(
    modules = [ACME2],
    sitename = "ACME2jl",
    pages = Any[
        "Home" => "index.md",
        "gettingstarted.md",
        "ug.md",
        "elements.md",
    ],
)

if v"1.0" ≤ VERSION < v"1.1" # deploy from 1.0
    deploydocs(
        repo = "github.com/markbennett95/ACME.jl.git",
        devbranch = "main",
    )
end
