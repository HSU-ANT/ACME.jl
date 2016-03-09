# Copyright 2015, 2016 Martin Holters
# See accompanying license file.

export resistor, capacitor, inductor, transformer, voltagesource, currentsource,
       voltageprobe, currentprobe, diode, bjt, opamp

resistor(r) = Element(mv=-1, mi=r)

capacitor(c) = Element(mv=[c;0], mi=[0; 1], mx=[-1;0], mxd=[0;-1])
inductor(l) = Element(mv=[1;0], mi=[0; l], mx=[0;-1], mxd=[-1;0])
transformer(l1, l2; coupling_coefficient=1,
            mutual_coupling=coupling_coefficient*sqrt(l1*l2)) =
    Element(mv=[1 0; 0 1; 0 0; 0 0],
            mi=[0 0; 0 0; l1 mutual_coupling; mutual_coupling l2],
            mx=[0 0; 0 0; -1 0; 0 -1], mxd=[-1 0; 0 -1; 0 0; 0 0],
            pins=[:primary1; :primary2; :secondary1; :secondary2])

function transformer(::Type{Val{:JA}}; D=2.4e-2, A=4.54e-5, ns=[],
                     a=14.1, α=5e-5, c=0.55/(1-0.55), k=17.8, Ms=2.75e5)
    const μ0 = 1.2566370614e-6
    nonlinear_eq =quote
        L_q1 = abs(q[1]) < 1e-4 ? q[1]/3 : coth(q[1])-1/q[1]
        Ld_q1 = abs(q[1]) < 1e-4 ? 1/3 : 1/q[1]^2-coth(q[1])^2+1
        Ld2_q1 = abs(q[1]) < 1e-3 ? -2/15*q[1] :
                                    2*coth(q[1])*(coth(q[1])^2 - 1) - 2/q[1]^3
        δ = q[3] > 0 ? 1.0 : -1.0

        Man = $(Ms)*L_q1
        Man_dq1 = $(Ms)*Ld_q1
        δM = sign(q[3]) == sign(Man - q[2]) ? 1.0 : 0.0

        den = δ*$(k)-$(α)*(Man-q[2])
        # at present, the error needs to be scaled to be comparable to those of
        # the other elements, hence the factor 1e-4/Ms
        res[1] = $(1e-4/Ms) * ($(1/(1+c)) * (δM*(Man-q[2])/den * q[3] +
                                             $(c*Ms/a)*(q[3]+$(α)*q[4])*Ld_q1)
                               - q[4])
        J[1,1] = $(1e-4/Ms) * $(1/(1+c)) * (δM*Man_dq1*δ*$(k)/den^2 * q[3] +
                                            $(c*Ms/a)*(q[3]+$(α)*q[4])*Ld2_q1)
        J[1,2] = $(1e-4/Ms) * -$(1/(1+c)) * δM*δ*$(k)/den^2 * q[3]
        J[1,3] = $(1e-4/Ms) * $(1/(1+c)) * (δM*(Man-q[2])/den + $(c*Ms/a)*Ld_q1)
        J[1,4] = $(1e-4/Ms) * ($(c/(1+c) * Ms/a * α)*Ld_q1 - 1)
    end
    Element(mv=[speye(length(ns)); spzeros(5, length(ns))],
            mi=[spzeros(length(ns), length(ns)); ns.'; spzeros(4, length(ns))],
            mx=[spzeros(length(ns), 2); -π*D 0; -1/a -α/a; 0 -1; 0 0; 0 0],
            mxd=[-μ0*A*ns -μ0*ns*A; 0 0; 0 0; 0 0; -1 0; 0 -1],
            mq=[zeros(length(ns)+1,4); eye(4)], nonlinear_eq = nonlinear_eq)
end

inductor(::Type{Val{:JA}}; n=230, args...) =
    transformer(Val{:JA}; ns=[n], args...)

voltagesource(v) = Element(mv=1, u0=v, pins=[:+; :-])
voltagesource() = Element(mv=1, mu=1, pins=[:+; :-])
currentsource(i) = Element(mi=-1, u0=i, pins=[:+; :-])
currentsource() = Element(mi=-1, mu=1, pins=[:+; :-])

voltageprobe() = Element(mi=1, pv=1, pins=[:+; :-])
currentprobe() = Element(mv=1, pi=1, pins=[:+; :-])

diode(;is::Number=1e-12, η::Number = 1) =
  Element(mv=[1;0], mi=[0;1], mq=[-1 0; 0 -1], pins=[:+; :-], nonlinear_eq =
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
                   pins = [:base; :emitter; :base; :collector])
end

opamp() = Element(mv=[0 0; 1 0], mi=[1 0; 0 0],
                  pins=[symbol("in+"); symbol(:"in-");
                        symbol(:"out+"); symbol(:"out-")])

function opamp(::Type{Val{:macak}}, gain, vomin, vomax)
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
                   pins=[symbol("in+"); symbol(:"in-");
                         symbol(:"out+"); symbol(:"out-")])
end

@Base.deprecate(opamp_ideal, opamp)
@Base.deprecate(opamp_macak(gain, vomin, vomax), opamp(Val{:macak}, gain, vomin, vomax))
