# Copyright 2015 Martin Holters
# See accompanying license file.

using ACME

function sallenkey(::Type{Circuit})
    r1 = resistor(10e3)
    r2 = resistor(10e3)
    c1 = capacitor(10e-9)
    c2 = capacitor(10e-9)
    j_in = voltagesource()
    u1 = opamp()
    #u1 = opamp(Val{:macak}, 1000, -4, 4)
    j_out = voltageprobe()

    circ = Circuit()
    connect!(circ, j_in[:-], :gnd)
    connect!(circ, j_in[:+], r1[1])
    connect!(circ, r1[2], r2[1], c1[1])
    connect!(circ, r2[2], c2[1], u1["in+"])
    connect!(circ, c2[2], u1["out-"], j_out[:-], :gnd)
    connect!(circ, u1["in-"], u1["out+"], j_out[:+])

    return circ
end

sallenkey(::Type{DiscreteModel}; fs=44100) =
    DiscreteModel(sallenkey(Circuit), 1/fs)
sallenkey(; fs=44100) = sallenkey(DiscreteModel, fs=fs)
