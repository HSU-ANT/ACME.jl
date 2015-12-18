# Copyright 2015 Martin Holters
# See accompanying license file.

export resistor, capacitor, inductor, transformer, voltagesource, currentsource,
       voltageprobe, diode, bjt, opamp_ideal, opamp_macak

resistor(r) = Element(mv=-1, mi=r)

capacitor(c) = Element(mv=[c,0], mi=[0, 1], mx=[-1,0], mxd=[0,-1])
inductor(l) = Element(mv=[1,0], mi=[0, l], mx=[0,-1], mxd=[-1,0])
transformer(l1, l2; coupling_coefficient=1,
            mutual_coupling=coupling_coefficient*sqrt(l1*l2)) =
    Element(mv=[1 0; 0 1; 0 0; 0 0],
            mi=[0 0; 0 0; l1 mutual_coupling; mutual_coupling l2],
            mx=[0 0; 0 0; -1 0; 0 -1], mxd=[-1 0; 0 -1; 0 0; 0 0],
            pins=[:primary1, :primary2, :secondary1, :secondary2])

voltagesource(v) = Element(mv=1, u0=v, pins=[:+, :-])
voltagesource() = Element(mv=1, mu=1, pins=[:+, :-])
currentsource(i) = Element(mi=-1, u0=i, pins=[:+, :-])
currentsource() = Element(mi=-1, mu=1, pins=[:+, :-])

voltageprobe() = Element(mi=1, pv=1, pins=[:+, :-])

diode(is::Number, η::Number = 1) =
  Element(mv=[1,0], mi=[0,1], mq=[-1 0; 0 -1], pins=[:+, :-], nonlinear_eq =
    quote
      let v = q[1], i = q[2], ex = exp(v*$(1./(25e-3 * η)))
        res[1] = $(is) * (ex - 1) - i
        J[1,1] = $(is/(25e-3 * η)) * ex
        J[1,2] = -1
      end
    end
  )

function bjt(typ; is=1e-12, η=1, isc=is, ise=is, ηc=η, ηe=η, βf=1000, βr=10)
    local polarity
    if typ == :npn
        polarity = 1
    elseif typ == :pnp
        polarity = -1
    else
        throw(ArgumentError(string("Unknown bjt type ", typ,
                                   ", must be :npn or :pnp")))
    end
    αf = -βf/(1+βf)
    αr = -βr/(1+βr)
    nonlinear_eq =
        quote
            let vE = q[1], vC = q[2], iE = q[3], iC = q[4],
                expE=exp(vE*$(1/(25e-3*ηe))), expC=exp(vC*$(1/(25e-3*ηc)))

                res[1] = $(ise) * (expE-1) + $(αr*isc) * (expC-1) - iE
                res[2] = $(αf*ise) * (expE-1) + $(isc) * (expC-1) - iC
                J[1,1] = $(ise/(25e-3*ηe)) * expE
                J[1,2] = $(αr*isc/(25e-3*ηc)) * expC
                J[1,3] = -1
                J[1,4] = 0
                J[2,1] = $(αf*ise/(25e-3*ηe)) * expE
                J[2,2] = $(isc/(25e-3*ηc)) * expC
                J[2,3] = 0
                J[2,4] = -1
            end
        end
    return Element(mv=[1 0; 0 1; 0 0; 0 0], mi = [0 0; 0 0; 1 0; 0 1],
                   mq = -polarity*speye(4), nonlinear_eq = nonlinear_eq,
                   pins = [:base, :emitter, :base, :collector])
end

opamp_ideal() = Element(mv=[0 0; 1 0], mi=[1 0; 0 0],
                        pins=[symbol("in+"), symbol(:"in-"),
                              symbol(:"out+"), symbol(:"out-")])
function opamp_macak(gain, vomin, vomax)
    offset = 0.5 * (vomin + vomax)
    scale = 0.5 * (vomax - vomin)
    nonlinear_eq =
        quote
            let vi = q[1], vo = q[2], vi_scaled = vi * $(gain/scale)

                res[1] = tanh(vi_scaled) * $(scale) - vo
                #res[1] = tanh(q[1]) * $(scale) - vo
                J[1,1] = $(gain) / cosh(vi_scaled)^2
                J[1,2] = -1
            end
        end
    return Element(mv=[0 0; 1 0; 0 1], mi=[1 0; 0 0; 0 0], mq=[0 0; -1 0; 0 -1],
                   u0=[0; 0; offset],
                   nonlinear_eq = nonlinear_eq,
                   pins=[symbol("in+"), symbol(:"in-"),
                         symbol(:"out+"), symbol(:"out-")])
end
