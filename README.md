# ACME.jl - Analog Circuit Modeling and Emulation for Julia

[![Join the chat at https://gitter.im/HSU-ANT/ACME.jl](https://badges.gitter.im/HSU-ANT/ACME.jl.svg)](https://gitter.im/HSU-ANT/ACME.jl?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://hsu-ant.github.io/ACME.jl/stable)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://hsu-ant.github.io/ACME.jl/dev/)
[![DOI](https://zenodo.org/badge/48224425.svg)](https://zenodo.org/badge/latestdoi/48224425)

[![Build Status](https://travis-ci.org/HSU-ANT/ACME.jl.svg?branch=develop)](https://travis-ci.org/HSU-ANT/ACME.jl)
[![codecov](https://codecov.io/gh/HSU-ANT/ACME.jl/branch/develop/graph/badge.svg)](https://codecov.io/gh/HSU-ANT/ACME.jl)
[![Coverage Status](https://coveralls.io/repos/github/HSU-ANT/ACME.jl/badge.svg?branch=develop)](https://coveralls.io/github/HSU-ANT/ACME.jl)

ACME is a [Julia](http://julialang.org/) package for the simulation of
electrical circuits, focusing on audio effect circuits. It allows to
programmatically describe a circuit in terms of elements and connections
between them and then automatically derive a model for the circuit. The model
can then be run on varying input data.

ACME is based on the method described in
[M. Holters, U. Zölzer, "A Generalized Method for the Derivation of Non-Linear
State-Space Models from Circuit
Schematics"](http://www.eurasip.org/Proceedings/Eusipco/Eusipco2015/papers/1570103545.pdf).

## Installation

If you have not done so already, [download and install
Julia](http://julialang.org/downloads/). (Any version starting with 1.0 should
be fine; earlier ACME versions support Julia starting with version 0.3.)

To install ACME, start Julia and run:

```Julia
Pkg.add("ACME")
```

This will download ACME and all of its dependencies.

## First steps

We will demonstrate ACME by modeling a simple diode clipper. The first step is
to load ACME:

```Julia
using ACME
```

Now we create the circuit description:

```Julia
circ = @circuit begin
    j_in = voltagesource()
    r1 = resistor(1e3)
    c1 = capacitor(47e-9)
    d1 = diode(is=1e-15)
    d2 = diode(is=1.8e-15)
    j_out = voltageprobe()
    j_in[+] ⟷ r1[1]
    j_in[-] ⟷ gnd
    r1[2] ⟷ c1[1] ⟷ d1[+] ⟷ d2[-] ⟷ j_out[+]
    gnd ⟷ c1[2] ⟷ d1[-] ⟷ d2[+] ⟷ j_out[-]
end
```

The first six lines inside the `begin`/`end` block instantiate circuit elements.
Specifying a `voltagesource()` sets up a voltage source as an input, i.e. the
voltage it sources will be specified when running the model. Alternatively, one
can instantiate a constant voltage source for say 9V with  `voltagesource(9)`.
The `resistor` and `capacitor` calls take the resistance in ohm and the
capacitance in farad, respectively, as arguments. For the `diode`, one may
specify the saturation current `is` as done here and/or the emission
coefficient `η`. Finally, desired outputs are denoted by adding probes to the
circuit; in this case a `voltageprobe()` will provide voltage as output.

The remaining four lines specify connections, either among element pins as in
`j_in[+] ⟷ r1[1]`, which connects the `+` pin of the input voltage to pin `1` of
the resistor, or among pins and named nets as in `j_in[-] ⟷ gnd`, which
connects the `-` pin of the input voltage source to a net named `gnd`. Note that
naming nets is only for the sake of readability; there is nothing special about
them and the names are arbitrary. As can be seen in the last two lines, multiple
pins can be connected at once.

It is also possible to specify connections following the element definition
(separated by commas), in which case the element name may be omitted. However,
one can only connect to elements defined before. Thus, above circuit could also
be entered as:

```Julia
circ = @circuit begin
    j_in = voltagesource(), [-] ⟷ gnd
    r1 = resistor(1e3), [1] ⟷ j_in[+]
    c1 = capacitor(47e-9), [1] ⟷ r1[2], [2] ⟷ gnd
    d1 = diode(is=1e-15), [+] ⟷ r1[2], [-] ⟷ gnd
    d2 = diode(is=1.8e-15), [+] ⟷ gnd, [-] ⟷ r1[2]
    j_out = voltageprobe(), [+] ⟷ r1[2], [-] ⟷ gnd
end
```

Now that the circuit has been set up, we need to turn it into a model. This
could hardly be any easier:

```Julia
model = DiscreteModel(circ, 1/44100)
```

The second argument specifies the sampling interval, the reciprocal of the
sampling rate, here assumed to be the typical 44100 Hz.

Now we can process some input data. It has to be provided as a matrix with one
row per input (just one in the example) and one column per sample. So for a
sinusoid at 1 kHz lasting one second, we do

```Julia
y = run!(model, sin(2π*1000/44100*(0:44099).'))
```

The output `y` now likewise is a matrix with one row for the one probe we have
added to the circuit and one column per sample.

In the likely event that you would like to process real audio data, take a look
at the [WAV](https://github.com/dancasimiro/WAV.jl) package for reading writing
WAV files.

Note that the solver used to solve the non-linear equation when running the
model saves solutions to use as starting points in the future. Model execution
will therefore become faster after an initial learning phase.  Nevertheless,
ACME is at present more geared towards computing all the model matrices than to
actually running the model. More complex circuits may run intolerably slow or
fail to run altogether.

## Moving on

There is some [documentation](https://hsu-ant.github.io/ACME.jl/dev/)
available for how
to use ACME. Additionally, you can take a look at the examples that can be found
in the `examples` directory below `Pkg.dir("ACME")`.

If you would like to extend and improve ACME, that's great! But unfortunately,
there is no developer documentation as of now, so you will to delve into the
source code to figure out how things work, or try to ask on
[gitter](https://gitter.im/HSU-ANT/ACME.jl).
