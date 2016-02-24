# Copyright 2016 Martin Holters
# See accompanying license file.

using ACME

# Model of "Der Birdie", see
# http://diy.musikding.de/wp-content/uploads/2013/06/birdieschalt.pdf
# for schematics.
# Note: Contrary to the schematics, in the real circuit as well as in this
# model, C5 and D1 are connected directly to the power supply; R6 is in series
# with LED1 only, which has been omitted from the model, though.

function birdie(::Type{Circuit}, vol)
    c1 = capacitor(2.2e-9)
    c3 = capacitor(2.2e-9)
    c5 = capacitor(100e-6)

    r1 = resistor(1e6)
    r2 = resistor(43e3)
    r3 = resistor(430e3)
    r4 = resistor(390)
    r5 = resistor(10e3)
    p1a = resistor(vol*100e3)
    p1b = resistor((1-vol)*100e3)

    d1 = diode(is=350e-12, η=1.6);
    t1 = bjt(:npn, isc=154.1e-15, ise=64.53e-15, ηc=1.10, ηe=1.06, βf=500, βr=12)

    j1 = voltagesource() # input
    j2 = voltageprobe() # output
    j3 = voltagesource(9) # 9V power supply

    circ = Circuit()
    connect!(circ, j3[:+], c5[2], d1[:-], :vcc)
    connect!(circ, j3[:-], c5[1], d1[:+], :gnd)
    connect!(circ, j1[:-], r1[2], :gnd)
    connect!(circ, j1[:+], r1[1], c1[1])
    connect!(circ, c1[2], r2[1], r3[1], t1[:base])
    connect!(circ, r2[2], :gnd)
    connect!(circ, r3[2], :vcc)
    connect!(circ, r4[1], t1[:emitter])
    connect!(circ, r4[2], :gnd)
    connect!(circ, r5[1], c3[1], t1[:collector])
    connect!(circ, r5[2], :vcc)
    connect!(circ, c3[2], p1b[1])
    connect!(circ, p1a[1], p1b[2])
    connect!(circ, p1b[2], j2[:+])
    connect!(circ, p1a[2], j2[:-], :gnd)

    return circ
end

birdie(::Type{DiscreteModel}, vol; fs=44100) =
    DiscreteModel(birdie(Circuit, vol), 1/fs)
birdie(vol; fs=44100) = birdie(DiscreteModel, vol, fs=fs)
