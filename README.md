# ACME.jl - Analog Circuit Modeling and Emulation for Julia

ACME is at present a proof of concept implementation of the method described in
[M. Holters, U. Zölzer, "A Generalized Method for the Derivation of Non-Linear
State-Space Models from Circuit
Schematics"](http://www.eurasip.org/Proceedings/Eusipco/Eusipco2015/papers/1570103545.pdf).

To get started despite the complete lack of any other documentation, try this
source code as a basis, which models a simple diode clipper:

```Julia
using ACME

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

model = DiscreteModel(circ, 1./44100)

y = run(model, sin(2π*1000/44100*(0:44099)'))

```

Note that the solver used to solve the non-linear equation when running the
model saves solutions to use as starting points in the future. Model execution
will therefore become faster after an initial learning phase.  Nevertheless,
ACME is at present more geared towards computing all the model matrices than to
actually running the model.
