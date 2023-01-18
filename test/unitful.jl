# Copyright 2023 Martin Holters
# See accompanying license file.

using Unitful: @u_str, DimensionError

@testset "element constructors supporting Unitful" begin
    @test resistor(10) == resistor(10u"Ω")
    @test resistor(1000) == resistor(1u"kΩ")
    @test_throws DimensionError resistor(1u"kV")

    @test potentiometer(10) == potentiometer(10u"Ω")
    @test potentiometer(1000) == potentiometer(1u"kΩ")
    @test_throws DimensionError potentiometer(1u"kV")
    @test potentiometer(10, 0.8) == potentiometer(10u"Ω", 0.8)
    @test potentiometer(1000, 0.8) == potentiometer(1u"kΩ", 0.8)
    @test_throws DimensionError potentiometer(1u"kV", 0.8)

    @test capacitor(10) == capacitor(10u"F")
    @test capacitor(1e-9) == capacitor(1.0u"nF")
    @test_throws DimensionError capacitor(1u"kV")

    @test inductor(10) == inductor(10u"H")
    @test inductor(1e-9) == inductor(1.0u"nH")
    @test_throws DimensionError inductor(1u"kV")

    @test inductor(Val{:JA}; D=10e-3) == inductor(Val{:JA}; D=10.0u"mm")
    @test inductor(Val{:JA}; D=10e-3, A=45e-6, n=200, a=14, α=4e-5, c=0.5, k=17, Ms=275e3) ==
        inductor(Val{:JA}; D=10.0u"mm", A=45e-6u"m^2", n=200, a=14u"A/m", α=4e-5, c=0.5, k=17u"A/m", Ms=275.0u"kA/m")
    @test_throws DimensionError inductor(Val{:JA}; D=10.0u"mm^2")

    @test transformer(10, 5) == transformer(10u"H", 5u"H")
    @test transformer(1e-9, 5e-9) == transformer(1.0u"nH", 5.0u"nH")
    @test transformer(4, 9; coupling_coefficient = 0.8) == transformer(4u"H", 9u"H"; coupling_coefficient = 0.8)
    @test transformer(1e-9, 5e-9; mutual_coupling = 2e-9) == transformer(1.0u"nH", 5.0u"nH"; mutual_coupling = 2.0u"nH")
    @test_throws DimensionError transformer(1u"kV", 2u"kV")
    @test_throws TypeError transformer(1u"nH", 2u"nH"; mutual_coupling = 2e-9)
    @test_throws DimensionError transformer(1u"nH", 2u"nH"; mutual_coupling = 2u"kV")

    @test transformer(Val{:JA}; D=10e-3) == transformer(Val{:JA}; D=10.0u"mm")
    @test transformer(Val{:JA}; D=10e-3, A=45e-6, ns=[20], a=14, α=4e-5, c=0.5, k=17, Ms=275e3) ==
        transformer(Val{:JA}; D=10.0u"mm", A=45e-6u"m^2", ns=[20], a=14u"A/m", α=4e-5, c=0.5, k=17u"A/m", Ms=275.0u"kA/m")
    @test_throws DimensionError transformer(Val{:JA}; D=10.0u"mm^2")

    @test voltagesource(10) == voltagesource(10u"V")
    @test voltagesource(1e-3) == voltagesource(1.0u"mV")
    @test voltagesource(5; rs=0.3) == voltagesource(5.0u"V"; rs=0.3u"Ω")
    @test voltagesource(; rs=0.3) == voltagesource(; rs=0.3u"Ω")
    @test_throws DimensionError voltagesource(1u"kA")
    @test_throws TypeError voltagesource(1u"V"; rs=1)
    @test_throws DimensionError voltagesource(1u"V"; rs=1u"V")
    @test_throws DimensionError voltagesource(; rs=1u"V")

    @test currentsource(10) == currentsource(10u"A")
    @test currentsource(1e-3) == currentsource(1.0u"mA")
    @test currentsource(5; gp=0.3) == currentsource(5.0u"A"; gp=0.3u"S")
    @test currentsource(; gp=0.3) == currentsource(; gp=0.3u"S")
    @test_throws DimensionError currentsource(1u"kV")
    @test_throws TypeError currentsource(1u"A"; gp=1)
    @test_throws DimensionError currentsource(1u"A"; gp=1u"A")
    @test_throws DimensionError currentsource(; gp=1u"A")

    @test voltageprobe(; gp=0.3) == voltageprobe(; gp=0.3u"S")
    @test_throws DimensionError voltageprobe(; gp=1u"A")

    @test currentprobe(; rs=0.3) == currentprobe(; rs=0.3u"Ω")
    @test_throws DimensionError currentprobe(; rs=1u"V")

    @test diode(; is=1e-15) == diode(; is=1.0u"fA")
    @test_throws DimensionError diode(; is=1u"V")
    @test_throws MethodError diode(; η=1u"V")

    @test bjt(:npn; is=1e-15) == bjt(:npn; is=1.0u"fA")
    @test bjt(:npn; isc=2e-12, ise=3e-12, ηc=1.1, ηe=1.2, βf=1000, βr=10, ile=4e-15,
              ilc=5e-15, ηcl=1.15, ηel=1.18, vaf=40, var=30, ikf=5e-12, ikr=6e-12, re=1,
              rc=2, rb=3) ==
        bjt(:npn; isc=2e-12u"A", ise=3e-12u"A", ηc=1.1, ηe=1.2, βf=1000, βr=10,
            ile=4e-15u"A", ilc=5e-15u"A", ηcl=1.15, ηel=1.18, vaf=40u"V", var=30u"V",
            ikf=5e-12u"A", ikr=6e-12u"A", re=1u"Ω", rc=2u"Ω", rb=3u"Ω")
    @test_throws DimensionError bjt(:npn; is=1u"V")
    @test_throws DimensionError bjt(:npn; βr=1u"A")

    @test mosfet(:n; vt=0.5) == mosfet(:n; vt=500.0u"mV")
    @test mosfet(:n; vt=0.5, α=25e-6, λ=0.1) ==
        mosfet(:n; vt=500.0u"mV", α=25.0e-6u"A/V^2", λ=0.1u"V^-1")
    @test mosfet(:n; vt=(-1.2454, -0.199, -0.0483), α=(0.0205, -0.0017), λ=0.1) ==
        mosfet(:n; vt=(-1.2454u"V", -0.199, -0.0483u"V^-1"), α=(0.0205u"A/V^2", -0.0017u"A/V^3"), λ=0.1u"V^-1")
    @test_throws DimensionError mosfet(:n; vt=500.0u"mA")
    @test_throws DimensionError mosfet(:n; vt=(500.0u"mV", 500.0u"mV"))

    @test opamp(Val{:macak}, 10, -3, 5) == opamp(Val{:macak}, 10, -3u"V", 5u"V")
    @test_throws DimensionError opamp(Val{:macak}, 10, -3u"V", 5u"A")
    @test_throws MethodError opamp(Val{:macak}, 10, -3u"V", 5)
    @test_throws MethodError opamp(Val{:macak}, 10u"V", -3u"V", 5u"V")
end
