# Copyright 2015, 2016 Martin Holters
# See accompanying license file.

using ACME

function diodeclipper(::Type{Circuit})
    j_in = voltagesource()
    r1 = resistor(1e3)
    c1 = capacitor(47e-9)
    d1 = diode(is=1e-15)
    d2 = diode(is=1.8e-15)
    j_out = voltageprobe()

    circ = Circuit()

    connect!(circ, j_in[:+], r1[1])
    connect!(circ, j_in[:-], :gnd)
    connect!(circ, r1[2], c1[1], d1[:+], d2[:-], j_out[:+])
    connect!(circ, :gnd, c1[2], d1[:-], d2[:+], j_out[:-])

    return circ
end

diodeclipper(::Type{DiscreteModel}; fs=44100) =
    DiscreteModel(diodeclipper(Circuit), 1/fs)
diodeclipper(; fs=44100) = diodeclipper(DiscreteModel, fs=fs)
