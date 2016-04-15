# Copyright 2016 Martin Holters
# See accompanying license file.

using ACME

# Model of "Der Super Over", see
# http://diy.musikding.de/wp-content/uploads/2013/06/superoverschalt.pdf
# for schematics.
# Note: R19 and LED D5 have been omitted from the model.

function superover(::Type{Circuit}; drive=nothing, tone=nothing, level=nothing, sym::Bool=false)
    r1 = resistor(2.2e6)
    r2 = resistor(10e3)
    r3 = resistor(470e3)
    r4 = resistor(10e3)
    r5 = resistor(100e3)
    r6 = resistor(33e3)
    r7 = resistor(4.7e3)
    r8 = resistor(10e3)
    r10 = resistor(10e3)
    r11 = resistor(470)
    r12 = resistor(4.7e3)
    r13 = resistor(1e6)
    r14 = resistor(10e3)
    r15 = resistor(1e3)
    r16 = resistor(100e3)
    r17 = resistor(33e3)
    r18 = resistor(33e3)
    r20 = resistor(22e3)
    if drive == nothing
        p1 = potentiometer(1e6)
    else
        p1 = potentiometer(1e6, drive)
    end
    if tone == nothing
        p2 = potentiometer(20e3)
    else
        p2 = potentiometer(20e3, tone)
    end
    if level == nothing
        p3 = potentiometer(10e3)
    else
        p3 = potentiometer(10e3, level)
    end

    c1 = capacitor(47e-9)
    c2 = capacitor(18e-9)
    c4 = capacitor(47e-9)
    c5 = capacitor(18e-9)
    c6 = capacitor(27e-9)
    c7 = capacitor(10e-9)
    c8 = capacitor(1e-3)
    c9 = capacitor(47e-9)
    c10 = capacitor(1e-6)
    c11 = capacitor(100e-6)
    c12 = capacitor(47e-6)

    d1 = diode(is=4e-9, η=2) # 1N914
    d2 = diode(is=3e-9, η=2) # 1N914
    d3 = diode(is=5e-9, η=2) # 1N914
    d4 = diode(is=12e-9, η=2) # 1N4001
    q1 = bjt(:npn, is=80e-15, βf=500, βr=10)
    q2 = bjt(:npn, is=80e-15, βf=500, βr=10)

    ic1a = opamp()
    ic1b = opamp()

    j1 = voltagesource() # input
    j2 = voltageprobe() # output
    j3 = voltagesource(9) # 9V power supply

    circ = Circuit()

    # power supply
    connect!(circ, j3[:+], :vcc)
    connect!(circ, j3[:-], :gnd)
    connect!(circ, d4[:-], c11[1], :vcc)
    connect!(circ, d4[:+], c11[2], :gnd)
    connect!(circ, :vcc, r17[1])
    connect!(circ, r17[2], r18[1], c12[1], :vb)
    connect!(circ, :gnd, r18[2], c12[2])

    # input stage
    connect!(circ, j1[:-], r1[2], :gnd)
    connect!(circ, j1[:+], r1[1], c1[1])
    connect!(circ, c1[2], r2[1])
    connect!(circ, r2[2], q1[:base], r3[1])
    connect!(circ, r3[2], :vb)
    connect!(circ, q1[:collector], :vcc)
    connect!(circ, q1[:emitter], r4[1], c2[1])
    connect!(circ, r4[2], :gnd)
    connect!(circ, c2[2], r5[1])
    connect!(circ, r5[2], :vb)

    # distortion stage
    connect!(circ, c2[2], ic1a["in+"])
    connect!(circ, ic1a["out+"], d3[:+], d1[:-], p1[2], p1[3])
    connect!(circ, ic1a["in-"], d2[:-], d1[:+], r6[1], c4[1])
    connect!(circ, ic1a["out-"], :gnd)
    if sym
        connect!(circ, d3[:-], d3[:+])
    end
    connect!(circ, r6[2], p1[1])
    connect!(circ, d2[:+], d3[:-])
    connect!(circ, c4[2], r7[1])
    connect!(circ, r7[2], :vb)

    # tone control stage
    connect!(circ, r8[1], ic1a["out+"])
    connect!(circ, r8[2], c5[1], ic1b["in+"])
    connect!(circ, c5[2], :gnd)
    connect!(circ, ic1b["out+"], r10[1], c7[1])
    connect!(circ, ic1b["in-"], r10[2], c7[2], p2[3])
    connect!(circ, ic1b["out-"], :gnd)
    connect!(circ, ic1b["in+"], p2[1])
    connect!(circ, p2[2], c6[1])
    connect!(circ, c6[2], r11[1])
    connect!(circ, r11[2], :gnd)

    # output stage
    connect!(circ, ic1b["out+"], c8[1])
    connect!(circ, c8[2], r12[1])
    connect!(circ, r12[2], p3[3])
    connect!(circ, p3[2], r20[1])
    connect!(circ, p3[1], :vb)
    connect!(circ, r20[2], c9[1])
    connect!(circ, c9[2], r13[1], q2[:base])
    connect!(circ, r13[2], :vb)
    connect!(circ, q2[:collector], :vcc)
    connect!(circ, q2[:emitter], r14[1], r15[1])
    connect!(circ, r14[2], :gnd)
    connect!(circ, r15[2], c10[1])
    connect!(circ, c10[2], r16[1])
    connect!(circ, r16[2], :gnd)
    connect!(circ, c10[2], j2[:+])
    connect!(circ, :gnd, j2[:-])

    return circ
end

superover{T<:DiscreteModel}(::Type{T}=DiscreteModel; drive=nothing, tone=nothing, level=nothing, sym::Bool=false, fs=44100) =
    T(superover(Circuit, drive=drive, tone=tone, level=level, sym=sym), 1/fs)
