# Copyright 2018 Martin Holters
# See accompanying license file.
@static if VERSION < v"0.6"
    ACMEdir = joinpath(dirname(@__FILE__), "..")
else
    ACMEdir = joinpath(@__DIR__, "..")
end
println("Checking copyright headers...")
for dirname in ("src", "examples", "test")
    dirname = joinpath(ACMEdir, dirname)
    for name in readdir(dirname)
        if endswith(name, ".jl")
            name = joinpath(dirname, name)
            years = sort!(unique(parse.([Int], readlines(`git log --format=%cd --date=format:%Y -- $name`))))
            if isempty(years)
                continue
            end
            println(name)
            open(name, "r") do io
                l = readline(io)
                #println(l)
                m = match(r"#\s*Copyright\s+(([0-9]+(,\s*)?)*)", l)
                if m === nothing
                    error("Missing copyright header in $name")
                end
                headyears = parse.([Int], strip.(split(m.captures[1], ',')))
                d = setdiff(years, headyears)
                if !isempty(d)
                    error("Missing years in copyright header of $name: $d")
                end
            end
        end
    end
end
