# Getting Started

## Installation

If you have not done so already, [download and install
Julia](http://julialang.org/downloads/). (Any version starting with 1.0 should
be fine; earlier ACME versions also support Julia 0.3 and later.)

To install ACME, start Julia and run:

```Julia
Pkg.add("ACME")
```

This will download ACME and all of its dependencies.

## First Steps

We will demonstrate ACME by modeling a simple diode clipper. The first step is
to load ACME:

```jldoctest firststeps; output = false
using ACME

# output

```

Now we create the circuit description:

```jldoctest firststeps; output = false, filter = r"(ACME\.)?Circuit\(.*"s
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

# output

Circuit(...)
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

```jldoctest firststeps; output = false, filter = r"(ACME\.)?Circuit\(.*"s
circ = @circuit begin
    j_in = voltagesource(), [-] ⟷ gnd
    r1 = resistor(1e3), [1] ⟷ j_in[+]
    c1 = capacitor(47e-9), [1] ⟷ r1[2], [2] ⟷ gnd
    d1 = diode(is=1e-15), [+] ⟷ r1[2], [-] ⟷ gnd
    d2 = diode(is=1.8e-15), [+] ⟷ gnd, [-] ⟷ r1[2]
    j_out = voltageprobe(), [+] ⟷ r1[2], [-] ⟷ gnd
end

# output

Circuit(...)
```

Now that the circuit has been set up, we need to turn it into a model. This
could hardly be any easier:

```jldoctest firststeps; output = false, filter = r"(ACME\.)?DiscreteModel{.*"s
model = DiscreteModel(circ, 1/44100)

# output

DiscreteModel{...}(...)
```

The second argument specifies the sampling interval, the reciprocal of the
sampling rate, here assumed to be the typical 44100 Hz.

Now we can process some input data. It has to be provided as a matrix with one
row per input (just one in the example) and one column per sample. So for a
sinusoid at 1 kHz lasting one second, we do:

```jldoctest firststeps; filter = r"\r?Running model:.*"
y = run!(model, sin.(2π*1000/44100*(0:44099)'))

# output

Running model: 100%|████████████████████████████████████| Time: 0:00:01
1×44100 Array{Float64,2}:
 0.0  0.0275964  0.0990996  0.195777  …  -0.537508  -0.462978  -0.36521
```

The output `y` now likewise is a matrix with one row for the one probe we have
added to the circuit and one column per sample.

More interesting circuits can be found in the examples located at
`Pkg.dir("ACME/examples")`.

In the likely event that you would like to process real audio data, take a look
at the [WAV](https://github.com/dancasimiro/WAV.jl) package for reading writing
WAV files.

Note that the solver used to solve the non-linear equation when running the
model saves solutions to use as starting points in the future. Model execution
will therefore become faster after an initial learning phase.  Nevertheless,
ACME is at present more geared towards computing all the model matrices than to
actually running the model. More complex circuits may run intolerably slow or
fail to run altogether.
