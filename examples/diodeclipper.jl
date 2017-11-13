# Copyright 2015, 2016, 2017 Martin Holters
# See accompanying license file.

using ACME

function diodeclipper(::Type{Circuit})
    @circuit begin
        j_in = voltagesource(), [-] ⟷ gnd
        r1 = resistor(1e3), [1] ⟷ j_in[+]
        c1 = capacitor(47e-9), [1] ⟷ r1[2], [2] ⟷ gnd
        d1 = diode(is=1e-15), [-] ⟷ gnd, [+] ⟷ r1[2]
        d2 = diode(is=1.8e-15), [-] ⟷ r1[2], [+] ⟷ gnd
        j_out = voltageprobe(), [-] ⟷ gnd, [+] ⟷ r1[2]
    end
end

diodeclipper{T<:DiscreteModel}(::Type{T}=DiscreteModel; fs=44100) =
    T(diodeclipper(Circuit), 1//fs)
