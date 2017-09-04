# Copyright 2017 Martin Holters
# See accompanying license file.

using ACME

# Model of "Das Plus v2", see
# http://diy.musikding.de/wp-content/uploads/2013/06/PlusV2schalt.pdf
# for schematics.
# Note: R8, Rled, and LED have been omitted from the model.

function dasplus(::Type{Circuit}; gain=nothing, level=nothing)
    r1 = resistor(1e6)
    r2 = resistor(10e3)
    r3 = resistor(4.7e3)
    r4 = resistor(1e6)
    r5 = resistor(1e6)
    r6 = resistor(1e6)
    r7 = resistor(1e6)
    r9 = resistor(10e3)
    if gain == nothing
        p2 = potentiometer(500e3)
    else
        p2 = potentiometer(500e3, gain)
    end
    if level == nothing
        p1 = potentiometer(100e3)
    else
        p1 = potentiometer(100e3, level)
    end

    c1 = capacitor(10e-9)
    c2 = capacitor(1e-9)
    c3 = capacitor(47e-9)
    c4 = capacitor(1e-6)
    c5 = capacitor(10e-12)
    c6 = capacitor(1e-6)
    c7 = capacitor(1e-9)
    c8 = capacitor(100e-6)

    d1 = diode(is=3.5e-9, η=1.95) # 1N4148
    d2 = diode(is=3.5e-9, η=1.95) # 1N4148
    d3 = diode(is=12e-9, η=2) # 1N4001

    ic1 = opamp(Val{:macak}, 1e5, 0.5, 8.5)

    j1 = voltagesource() # input
    j2 = voltageprobe() # output
    j3 = voltagesource(9) # 9V power supply

    circ = Circuit()

    # power supply
    connect!(circ, j3[:+], :vcc)
    connect!(circ, j3[:-], :gnd)
    connect!(circ, d3[:-], c8[1], :vcc)
    connect!(circ, d3[:+], c8[2], :gnd)
    connect!(circ, :vcc, r7[1])
    connect!(circ, r7[2], r5[1], c4[1], :vb)
    connect!(circ, :gnd, r5[2], c4[2])

    # singal path
    connect!(circ, j1[:+], r1[1], c2[1], c1[1])
    connect!(circ, j1[:-], r1[2], c2[2], :gnd)
    connect!(circ, c1[2], r2[1])
    connect!(circ, r2[2], r6[1], ic1["in+"])
    connect!(circ, r6[2], :vb)
    connect!(circ, ic1["out+"], c5[1], c6[1], r4[1])
    connect!(circ, ic1["in-"], c5[2], r4[2], c3[1])
    connect!(circ, c3[2], r3[1])
    connect!(circ, r3[2], p2[3])
    connect!(circ, #=p2[1],=# p2[2], :gnd)

    connect!(circ, ic1["out-"], :gnd)
    connect!(circ, c6[2], r9[1])
    connect!(circ, r9[2], d1[:+], d2[:-], c7[1], p1[3])
    connect!(circ, d1[:-], d2[:+], c7[2], p1[1], j2[:-], :gnd)
    connect!(circ, j2[:+], p1[2])

    return circ
end

dasplus{S}(::Type{S}; gain=nothing, level=nothing, fs=44100) =
    DiscreteModel(dasplus(Circuit, gain=gain, level=level), 1//fs, S)
dasplus(;gain=nothing, level=nothing, fs=44100) =
    DiscreteModel(dasplus(Circuit, gain=gain, level=level), 1//fs)
