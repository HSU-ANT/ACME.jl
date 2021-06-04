# Copyright 2016, 2017, 2018, 2021 Martin Holters
# See accompanying license file.

using ACME2

# Model of "Der Birdie", see
# http://diy.musikding.de/wp-content/uploads/2013/06/birdieschalt.pdf
# for schematics.
# Note: Contrary to the schematics, in the real circuit as well as in this
# model, C5 and D1 are connected directly to the power supply; R6 is in series
# with LED1 only, which has been omitted from the model, though.

function birdie(::Type{Circuit}; vol=nothing)
    @circuit begin
        j3 = voltagesource(9), [-] ⟷ gnd, [+] ⟷ vcc # 9V power supply
        c5 = capacitor(100e-6), [1] ⟷ gnd, [2] ⟷ vcc
        d1 = diode(is=350e-12, η=1.6), [-] ⟷ vcc, [+] ⟷ gnd
        j1 = voltagesource(), [-] ⟷ gnd # input
        r1 = resistor(1e6), [1] ⟷ j1[+], [2] ⟷ gnd
        c1 = capacitor(2.2e-9), [1] ⟷ j1[+]
        r2 = resistor(43e3), [1] ⟷ c1[2], [2] ⟷ gnd
        r3 = resistor(430e3), [1] ⟷ c1[2], [2] ⟷ vcc
        t1 = bjt(:npn, isc=154.1e-15, ise=64.53e-15, ηc=1.10, ηe=1.06, βf=500, βr=12),
             [base] ⟷ c1[2]
        r4 = resistor(390), [1] ⟷ t1[emitter], [2] ⟷ gnd
        r5 = resistor(10e3), [1] ⟷ t1[collector], [2] ⟷ vcc
        c3 = capacitor(2.2e-9), [1] ⟷ t1[collector]
        p1 = potentiometer(100e3, (vol == nothing ? () : (vol,))...), [1] ⟷ gnd, [3] ⟷ c3[2]
        j2 = voltageprobe(), [-] ⟷ gnd, [+] ⟷ p1[2]# output
    end
end

birdie(::Type{DiscreteModel}=DiscreteModel; vol=nothing, fs=44100, solver=nothing) =
    solver === nothing ?
        DiscreteModel(birdie(Circuit, vol=vol), 1//fs) :
        DiscreteModel(birdie(Circuit, vol=vol), 1//fs, solver)
