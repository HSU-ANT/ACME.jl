# Copyright 2016, 2017, 2018, 2021 Martin Holters
# See accompanying license file.

using ACME2

# Model of "Der Super Over", see
# http://diy.musikding.de/wp-content/uploads/2013/06/superoverschalt.pdf
# for schematics.
# Note: R19 and LED D5 have been omitted from the model.

function superover(::Type{Circuit}; drive=nothing, tone=nothing, level=nothing, sym::Bool=false)
    circ = @circuit begin
        # power supply
        j3 = voltagesource(9), [+] ⟷ vcc, [-] ⟷ gnd # 9V power supply
        d4 = diode(is=12e-9, η=2), [-] ⟷ vcc, [+] ⟷ gnd # 1N4001
        c11 = capacitor(100e-6), [1] ⟷ vcc, [2] ⟷ gnd
        r17 = resistor(33e3), [1] ⟷ vcc, [2] ⟷ vb
        r18 = resistor(33e3), [1] ⟷ vb, [2] ⟷ gnd
        c12 = capacitor(47e-6), [1] ⟷ vb, [2] ⟷ gnd

        # input stage
        j1 = voltagesource(), [-] ⟷ gnd # input
        r1 = resistor(2.2e6), [1] ⟷ j1[+], [2] ⟷ gnd
        c1 = capacitor(47e-9), [1] ⟷ j1[+]
        r2 = resistor(10e3), [1] ⟷ c1[2]
        r3 = resistor(470e3), [1] ⟷ r2[2], [2] ⟷ vb
        q1 = bjt(:npn, is=80e-15, βf=500, βr=10), [base] ⟷ r2[2], [collector] ⟷ vcc
        r4 = resistor(10e3), [1] ⟷ q1[emitter], [2] ⟷ gnd
        c2 = capacitor(18e-9), [1] ⟷ q1[emitter]
        r5 = resistor(100e3), [1] ⟷ c2[2], [2] ⟷ vb

        # distortion stage
        ic1a = opamp(), ["in+"] ⟷ c2[2], ["out-"] ⟷ gnd
        d1 = diode(is=4e-9, η=2), [-] ⟷ ic1a["out+"], [+] ⟷ ic1a["in-"] # 1N914
        d2 = diode(is=3e-9, η=2), [-] ⟷ ic1a["in-"] # 1N914
        d3 = diode(is=5e-9, η=2), [+] ⟷ ic1a["out+"], [-] ⟷ d2[+] # 1N914
        p1 = potentiometer(1e6, (drive == nothing ? () : (drive,))...), [2] ⟷ [3] ⟷ ic1a["out+"]
        r6 = resistor(33e3), [1] ⟷ ic1a["in-"], [2] ⟷ p1[1]
        c4 = capacitor(47e-9), [1] ⟷ ic1a["in-"]
        r7 = resistor(4.7e3), [1] ⟷ c4[2], [2] ⟷ vb

        # tone control stage
        r8 = resistor(10e3), [1] ⟷ ic1a["out+"]
        ic1b = opamp(), ["in+"] ⟷ r8[2], ["out-"] ⟷ gnd
        c5 = capacitor(18e-9), [1] ⟷ ic1b["in+"], [2] ⟷ gnd
        r10 = resistor(10e3), [1] ⟷ ic1b["out+"], [2] ⟷ ic1b["in-"]
        c7 = capacitor(10e-9), [1] ⟷ ic1b["out+"], [2] ⟷ ic1b["in-"]
        p2 = potentiometer(20e3, (tone == nothing ? () : (tone,))...), [1] ⟷ ic1b["in+"], [3] ⟷ ic1b["in-"]
        c6 = capacitor(27e-9), [1] ⟷ p2[2]
        r11 = resistor(470), [1] ⟷ c6[2], [2] ⟷ gnd

        # output stage
        c8 = capacitor(1e-3), [1] ⟷ ic1b["out+"]
        r12 = resistor(4.7e3), [1] ⟷ c8[2]
        p3 = potentiometer(10e3, (level == nothing ? () : (level,))...), [1] ⟷ vb, [3] ⟷ r12[2]
        r20 = resistor(22e3), [1] ⟷ p3[2]
        c9 = capacitor(47e-9), [1] ⟷ r20[2]
        r13 = resistor(1e6), [1] ⟷ c9[2], [2] ⟷ vb
        q2 = bjt(:npn, is=80e-15, βf=500, βr=10), [base] ⟷ c9[2], [collector] ⟷ vcc
        r14 = resistor(10e3), [1] ⟷ q2[emitter], [2] ⟷ gnd
        r15 = resistor(1e3), [1] ⟷ q2[emitter]
        c10 = capacitor(1e-6), [1] ⟷ r15[2]
        r16 = resistor(100e3), [1] ⟷ c10[2], [2] ⟷ gnd
        j2 = voltageprobe(), [+] ⟷ c10[2], [-] ⟷ gnd # output
    end

    if sym
        connect!(circ, (:d3, :-), (:d3, :+))
    end

    return circ
end

superover(::Type{DiscreteModel}=DiscreteModel;drive=nothing, tone=nothing, level=nothing, sym::Bool=false, fs=44100, solver=nothing) =
    solver === nothing ?
        DiscreteModel(superover(Circuit, drive=drive, tone=tone, level=level, sym=sym), 1//fs) :
        DiscreteModel(superover(Circuit, drive=drive, tone=tone, level=level, sym=sym), 1//fs, solver)
