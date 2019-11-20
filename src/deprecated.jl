# Copyright 2019 Martin Holters
# See accompanying license file.

Base.@deprecate(composite_element(circ::Circuit, pins::Vector{<:Pair}),
    composite_element(circ,  pinmap=Dict(pins...)))
