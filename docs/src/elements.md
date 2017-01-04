# Element Reference
## Passives

```@docs
resistor
capacitor
inductor(l)
inductor(::Type{Val{:JA}})
transformer(l1, l2; coupling_coefficient=1,mutual_coupling=coupling_coefficient*sqrt(l1*l2))
transformer(::Type{Val{:JA}})
```

## Independent Sources
```@docs
voltagesource
currentsource
```

## Probes
```@docs
voltageprobe
currentprobe
```

## Semiconductors
```@docs
diode
bjt
```

## Integrated Circuits

```@docs
opamp()
opamp(::Type{Val{:macak}}, gain, vomin, vomax)
```
