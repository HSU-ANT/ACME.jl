# ACME.jl - Analog Circuit Modeling and Emulation for Julia

ACME is a [Julia](http://julialang.org/) package for the simulation of
electrical circuits, focusing on audio effect circuits. It allows to
programmatically describe a circuit in terms of elements and connections
between them and then automatically derive a model for the circuit. The model
can then be run on varying input data.

ACME is based on the method described in [M. Holters, U. Zölzer, "A Generalized
Method for the Derivation of Non-Linear State-Space Models from Circuit
Schematics"](http://www.eurasip.org/Proceedings/Eusipco/Eusipco2015/papers/1570103545.pdf).

```@contents
Pages = [
      "gettingstarted.md",
      "ug.md",
      "elements.md",
  ]
```
