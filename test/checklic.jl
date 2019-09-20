# Copyright 2018, 2019 Martin Holters
# See accompanying license file.
ACMEdir = joinpath(@__DIR__, "..")
println("Checking copyright headers...")
for dirname in ("src", "examples", "test")
    dirname = joinpath(ACMEdir, dirname)
    for name in readdir(dirname)
        if endswith(name, ".jl")
            name = joinpath(dirname, name)
            local years
            try
                years = sort!(unique(parse.([Int], readlines(`git log --format=%cd --date=format:%Y -- $name`))))
            catch e
                @warn e
                continue
            end
            if isempty(years)
                continue
            end
            println(name)
            open(name, "r") do io
                l = readline(io)
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
