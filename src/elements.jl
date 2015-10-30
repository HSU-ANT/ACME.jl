export resistor, capacitor, voltagesource, voltageprobe, diode

resistor(r) = Element(mv=-1, mi=r)
capacitor(c) = Element(mv=[c,0], mi=[0, 1], mx=[-1,0], mxd=[0,-1])
voltagesource(v) = Element(mv=-1, u0=v, pins=[:+, :-])
voltagesource() = Element(mv=1, mu=1, pins=[:+, :-])
voltageprobe() = Element(mi=1, pv=1, pins=[:+, :-])
diode(is::Number, η::Number = 1) =
  Element(mv=[1,0], mi=[0,1], mq=[-1 0; 0 -1], pins=[:+, :-], nonlinear_eq =
    quote
      v::Float64 = q[1]
      i::Float64 = q[2]
      ex::Float64 = exp(v*$(1./(25e-3 * η)))
      res[1] = $(is) * (ex - 1) - i
      J[1,1] = $(is/(25e-3 * η)) * ex
      J[1,2] = -1
    end
  )
