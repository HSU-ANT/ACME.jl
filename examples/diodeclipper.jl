# Copyright 2015 Martin Holters
# See accompanying license file.

using ACME

function diodeclipper()
    j_in = voltagesource()
    r1 = resistor(1e3)
    c1 = capacitor(47e-9)
    d1 = diode(1e-15)
    d2 = diode(1.8e-15)
    j_out = voltageprobe()

    circ = Circuit()

    connect!(circ, j_in[:+], r1[1])
    connect!(circ, j_in[:-], :gnd)
    connect!(circ, r1[2], c1[1], d1[:+], d2[:-], j_out[:+])
    connect!(circ, :gnd, c1[2], d1[:-], d2[:+], j_out[:-])

    return DiscreteModel(circ, 1./44100)
end
