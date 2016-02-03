# Copyright 2016 Martin Holters
# See accompanying license file.

using ACME

function birdie(vol)
    v_supply = voltagesource(9)
    c5 = capacitor(100e-6)
    d1 = diode(is=350e-12, η=1.6);
    j_in = voltagesource()
    r1 = resistor(1e6)
    c1 = capacitor(2.2e-9)
    r2 = resistor(43e3)
    r3 = resistor(430e3)
    q1 = bjt(:npn, isc=154.1e-15, ise=64.53e-15, ηc=1.10, ηe=1.06, βf=500, βr=12)
    r4 = resistor(390)
    r5 = resistor(10e3)
    c3 = capacitor(2.2e-9)
    p1a = resistor(vol*100e3)
    p1b = resistor((1-vol)*100e3)
    j_out = voltageprobe()

    circ = Circuit()
    connect!(circ, v_supply[:+], c5[2], d1[:-], :vcc)
    connect!(circ, v_supply[:-], c5[1], d1[:+], :gnd)
    connect!(circ, j_in[:-], r1[2], :gnd)
    connect!(circ, j_in[:+], r1[1], c1[1])
    connect!(circ, c1[2], r2[1], r3[1], q1[:base])
    connect!(circ, r2[2], :gnd)
    connect!(circ, r3[2], :vcc)
    connect!(circ, r4[1], q1[:emitter])
    connect!(circ, r4[2], :gnd)
    connect!(circ, r5[1], c3[1], q1[:collector])
    connect!(circ, r5[2], :vcc)
    connect!(circ, c3[2], p1b[1])
    connect!(circ, p1a[1], p1b[2])
    connect!(circ, p1b[2], j_out[:+])
    connect!(circ, p1a[2], j_out[:-], :gnd)

    return DiscreteModel(circ, 1./44100)
end
