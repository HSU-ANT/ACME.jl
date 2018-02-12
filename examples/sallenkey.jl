# Copyright 2015, 2016, 2017, 2018 Martin Holters
# See accompanying license file.

using ACME

function sallenkey(::Type{Circuit})
    @circuit begin
        j_in = voltagesource(), [-] ⟷ gnd
        r1 = resistor(10e3), [1] ⟷ j_in[+]
        r2 = resistor(10e3), [1] ⟷ r1[2]
        c1 = capacitor(10e-9), [1] ⟷ r1[2]
        u1 = opamp(), ["in+"] ⟷ r2[2], ["in-"] ⟷ ["out+"] ⟷ c1[2], ["out-"] ⟷ gnd
        #u1 = opamp(Val{:macak}, 1000, -4, 4)
        c2 = capacitor(10e-9), [1] ⟷ u1["in+"], [2] ⟷ gnd
        j_out = voltageprobe(), [-] ⟷ gnd, [+] ⟷ u1["out+"]
    end
end

sallenkey{T<:DiscreteModel}(::Type{T}=DiscreteModel; fs=44100) =
    T(sallenkey(Circuit), 1//fs)
