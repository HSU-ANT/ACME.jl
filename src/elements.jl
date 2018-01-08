# Copyright 2015, 2016, 2017 Martin Holters
# See accompanying license file.

export resistor, potentiometer, capacitor, inductor, transformer,
       voltagesource, currentsource,
       voltageprobe, currentprobe, diode, bjt, opamp


"""
    resistor(r)

Creates a resistor obeying Ohm’s law. The resistance `r` has to be given in Ohm.

Pins: `1`, `2`
"""
resistor(r) = Element(mv=-1, mi=r)

potentiometer(r, pos) = Element(mv=Matrix{Int}(-I, 2, 2), mi=[r*pos 0; 0 r*(1-pos)],
                                pins=[1, 2, 2, 3])
potentiometer(r) =
    Element(mv=[Matrix{Int}(I, 2, 2); zeros(3, 2)],
            mi=[zeros(2, 2); Matrix{Int}(I, 2, 2); zeros(1, 2)],
            mq=Matrix{Int}(-I, 5, 5), mu=[zeros(4, 1); -1],
            nonlinear_eq = quote
                let v1=q[1], v2=q[2], i1=q[3], i2=q[4], pos=q[5]
                    res[1] = v1 - $(r)*pos*i1
                    res[2] = v2 - $(r)*(1-pos)*i2
                    J[1,1] = 1
                    J[1,2] = 0
                    J[1,3] = $(-r)*pos
                    J[1,4] = 0
                    J[1,5] = $(-r)*i1
                    J[2,1] = 0
                    J[2,2] = 1
                    J[2,3] = 0
                    J[2,4] = $(-r)*(1-pos)
                    J[2,5] = $(-r)*i2
                end
            end,
            pins=[1, 2, 2, 3])

"""
    capacitor(c)

Creates a capacitor. The capacitance `c` has to be given in Farad.

Pins: `1`, `2`
"""
capacitor(c) = Element(mv=[c;0], mi=[0; 1], mx=[-1;0], mxd=[0;-1])

"""
    inductor(l)

Creates an inductor. The inductance `l` has to be given in Henri.

Pins: `1`, `2`
"""
inductor(l) = Element(mv=[1;0], mi=[0; l], mx=[0;-1], mxd=[-1;0])

"""
    transformer(l1, l2; coupling_coefficient=1, mutual_coupling=coupling_coefficient*sqrt(l1*l2))

Creates a transformer with two windings having inductances. The primary
self-inductance `l1` and the secondary self-inductance `l2` have to be given in
Henri. The coupling can either be specified using `coupling_coefficient` (0 is
not coupled, 1 is closely coupled) or by `mutual_coupling`, the mutual
inductance in Henri, where the latter takes precedence if both are given.

Pins: `1` and `2` for primary winding, `3` and `4` for secondary winding
"""
transformer(l1, l2; coupling_coefficient=1,
            mutual_coupling=coupling_coefficient*sqrt(l1*l2)) =
    Element(mv=[1 0; 0 1; 0 0; 0 0],
            mi=[0 0; 0 0; l1 mutual_coupling; mutual_coupling l2],
            mx=[0 0; 0 0; -1 0; 0 -1], mxd=[-1 0; 0 -1; 0 0; 0 0],
            pins=[:primary1; :primary2; :secondary1; :secondary2])

"""
    transformer(Val{:JA}; D, A, ns, a, α, c, k, Ms)

Creates a non-linear transformer based on the Jiles-Atherton model of
magnetization assuming a toroidal core thin compared to its diameter. The
parameters are set using named arguments:

| parameter | description |
|:--------- |:------------|
| `D`  | Torus diameter (in meters) |
| `A`  | Torus cross-sectional area (in square-meters) |
| `ns` | Windings' number of turns as a vector with one entry per winding |
| `a`  | Shape parameter of the anhysteretic magnetization curve (in Ampere-per-meter) |
| `α`  | Inter-domain coupling |
| `c`  | Ratio of the initial normal to the initial anhysteretic differential susceptibility |
| `k`  | amount of hysteresis (in Ampere-per-meter) |
| `Ms` | saturation magnetization (in Ampere-per-meter) |

A detailed discussion of the parameters can be found in D. C. Jiles and D. L.
Atherton, “Theory of ferromagnetic hysteresis,” J. Magn. Magn. Mater., vol.
61, no. 1–2, pp. 48–60, Sep. 1986 and J. H. B. Deane, “Modeling the dynamics
of nonlinear inductor circuits,” IEEE Trans. Magn., vol. 30, no. 5, pp.
2795–2801, 1994, where the definition of `c` is taken from the latter. The ACME
implementation is discussed in [M. Holters, U. Zölzer, "Circuit Simulation with
Inductors and Transformers Based on the Jiles-Atherton Model of
Magnetization"](http://ant-s4.unibw-hamburg.de/paper-archive/2016/dafxpapers/08-DAFx-16_paper_10-PN.pdf).

Pins: `1` and `2` for primary winding, `3` and `4` for secondary winding, and so
on
"""
function transformer(::Type{Val{:JA}}; D=2.4e-2, A=4.54e-5, ns=[],
                     a=14.1, α=5e-5, c=0.55, k=17.8, Ms=2.75e5)
    μ0 = 1.2566370614e-6
    nonlinear_eq = quote
        coth_q1 = coth(q[1])
        a_q1 = abs(q[1])
        L_q1 = a_q1 < 1e-4 ? q[1]/3 : coth_q1 - 1/q[1]
        Ld_q1 = a_q1 < 1e-4 ? 1/3 : 1/q[1]^2 - coth_q1^2 + 1
        Ld2_q1 = a_q1 < 1e-3 ? -2/15*q[1] : 2*coth_q1*(coth_q1^2 - 1) - 2/q[1]^3
        δ = q[3] > 0 ? 1.0 : -1.0

        Man = $(Ms)*L_q1
        δM = sign(q[3]) == sign(Man - q[2]) ? 1.0 : 0.0

        den = δ*$(k*(1-c))-$(α)*(Man-q[2])
        # at present, the error needs to be scaled to be comparable to those of
        # the other elements, hence the factor 1e-4/Ms
        res[1] = $(1e-4/Ms) * ($(1-c) * δM*(Man-q[2])/den * q[3]
                               + $(c*Ms/a)*(q[3]+$(α)*q[4])*Ld_q1 - q[4])
        J[1,1] = $(1e-4/Ms) * ($((1-c)^2*k*Ms) * δM*Ld_q1*δ/den^2 * q[3]
                               + $(c*Ms/a)*(q[3]+$(α)*q[4])*Ld2_q1)
        J[1,2] = $(1e-4/Ms) * -$((1-c)^2*k) * δM*δ/den^2 * q[3]
        J[1,3] = $(1e-4/Ms) * ($(1-c) * δM*(Man-q[2])/den + $(c*Ms/a)*Ld_q1)
        J[1,4] = $(1e-4/Ms) * ($(c * Ms/a * α)*Ld_q1 - 1)
    end
    Element(mv=[speye(length(ns)); spzeros(5, length(ns))],
            mi=[spzeros(length(ns), length(ns)); ns'; spzeros(4, length(ns))],
            mx=[spzeros(length(ns), 2); -π*D 0; -1/a -α/a; 0 -1; 0 0; 0 0],
            mxd=[-μ0*A*ns -μ0*ns*A; 0 0; 0 0; 0 0; -1 0; 0 -1],
            mq=[zeros(length(ns)+1,4); eye(4)], nonlinear_eq = nonlinear_eq)
end


"""
    inductor(Val{:JA}; D, A, n, a, α, c, k, Ms)

Creates a non-linear inductor based on the Jiles-Atherton model of
magnetization assuming a toroidal core thin compared to its diameter. The
parameters are set using named arguments:

| parameter | description |
|:--------- |:------------|
| `D`  | Torus diameter (in meters) |
| `A`  | Torus cross-sectional area (in square-meters) |
| `n`  | Winding's number of turns |
| `a`  | Shape parameter of the anhysteretic magnetization curve (in Ampere-per-meter) |
| `α`  | Inter-domain coupling |
| `c`  | Ratio of the initial normal to the initial anhysteretic differential susceptibility |
| `k`  | amount of hysteresis (in Ampere-per-meter) |
| `Ms` | saturation magnetization (in Ampere-per-meter) |

A detailed discussion of the paramters can be found in D. C. Jiles and D. L.
Atherton, “Theory of ferromagnetic hysteresis,” J. Magn. Magn. Mater., vol.
61, no. 1–2, pp. 48–60, Sep. 1986 and J. H. B. Deane, “Modeling the dynamics
of nonlinear inductor circuits,” IEEE Trans. Magn., vol. 30, no. 5, pp.
2795–2801, 1994, where the definition of `c` is taken from the latter. The ACME
implementation is discussed in [M. Holters, U. Zölzer, "Circuit Simulation with
Inductors and Transformers Based on the Jiles-Atherton Model of
Magnetization"](http://ant-s4.unibw-hamburg.de/paper-archive/2016/dafxpapers/08-DAFx-16_paper_10-PN.pdf).

Pins: `1`, `2`
"""
inductor(::Type{Val{:JA}}; n=230, args...) =
    transformer(Val{:JA}; ns=[n], args...)

"""
    voltagesource(; rs=0)
    voltagesource(v; rs=0)

Creates a voltage source. The source voltage `v` has to be given in Volt. If
omitted, the source voltage will be an input of the circuit. Optionally, an
internal series resistance `rs` (in Ohm) can be given which defaults to zero.

Pins: `+` and `-` with `v` being measured from `+` to `-`
"""
function voltagesource end
voltagesource(v; rs=0) = Element(mv=1, mi=-rs, u0=v, pins=[:+; :-])
voltagesource(; rs=0) = Element(mv=1, mi=-rs, mu=1, pins=[:+; :-])

"""
    currentsource(; gp=0)
    currentsource(i; gp=0)

Creates a current source. The source current `i` has to be given in Ampere. If
omitted, the source current will be an input of the circuit. Optionally, an
internal parallel conductance `gp` (in Ohm⁻¹) can be given which defaults to
zero.

Pins: `+` and `-` where `i` measures the current leaving source at the `+` pin
"""
function currentsource end
currentsource(i; gp=0) = Element(mv=gp, mi=-1, u0=i, pins=[:+; :-])
currentsource(; gp=0) = Element(mv=gp, mi=-1, mu=1, pins=[:+; :-])

"""
    voltageprobe()

Creates a voltage probe, providing the measured voltage as a circuit output.
Optionally, an internal parallel conductance `gp` (in Ohm⁻¹) can be given which
defaults to zero.

Pins: `+` and `-` with the output voltage being measured from `+` to `-`
"""
voltageprobe(;gp=0) = Element(mv=-gp, mi=1, pv=1, pins=[:+; :-])

"""
    currentprobe()

Creates a current probe, providing the measured current as a circuit output.
Optionally, an internal series resistance `rs` (in Ohm) can be given which
defaults to zero.

Pins: `+` and `-` with the output current being the current entering the probe
at `+`
"""
currentprobe(;rs=0) = Element(mv=1, mi=-rs, pi=1, pins=[:+; :-])

doc"""
    diode(;is=1e-12, η = 1)

Creates a diode obeying Shockley's law
$i=I_S\cdot(e^{v/(\eta v_T)}-1)$ where $v_T$ is fixed at 25 mV.
The reverse saturation current `is` has to be given in Ampere, the emission
coefficient `η` is unitless.

Pins: `+` (anode) and `-` (cathode)
"""
diode(;is::Real=1e-12, η::Real = 1) =
  Element(mv=[1;0], mi=[0;1], mq=[-1 0; 0 -1], pins=[:+; :-], nonlinear_eq =
    quote
      let v = q[1], i = q[2], ex = exp(v*$(1 / (25e-3 * η)))
        res[1] = $(is) * (ex - 1) - i
        J[1,1] = $(is/(25e-3 * η)) * ex
        J[1,2] = -1
      end
    end
  )

doc"""
    bjt(typ; is=1e-12, η=1, isc=is, ise=is, ηc=η, ηe=η, βf=1000, βr=10,
        ile=0, ilc=0, ηcl=ηc, ηel=ηe, vaf=Inf, var=Inf, ikf=Inf, ikr=Inf)

Creates a bipolar junction transistor obeying the Gummel-Poon model

$i_f = \frac{\beta_f}{1+\beta_f} I_{S,E} \cdot (e^{v_E/(\eta_E v_T)}-1)$

$i_r = \frac{\beta_r}{1+\beta_r} I_{S,C} \cdot (e^{v_C/(\eta_C v_T)}-1)$

$i_{cc} = \frac{2(1-\frac{V_E}{V_{ar}}-\frac{V_C}{V_{af}})}
               {1+\sqrt{1+4(\frac{i_f}{I_{KF}}+\frac{i_r}{I_{KR}})}}
          (i_f - i_r)$

$i_{BE} = \frac{1}{\beta_f} i_f + I_{L,E} \cdot (e^{v_E/(\eta_{EL} v_T)}-1)$

$i_{BC} = \frac{1}{\beta_r} i_r + I_{L,C} \cdot (e^{v_C/(\eta_{CL} v_T)}-1)$

$i_E = i_{cc} + i_{BE} \qquad i_C=-i_{cc} + i_{BC}$

where $v_T$ is fixed at 25 mV. For

$I_{L,E}=I_{L,C}=0,\quad V_{ar}=V_{af}=I_{KF}=I_{KR}=∞,$

this reduces to the Ebers-Moll equation

$i_E = I_{S,E} \cdot (e^{v_E/(\eta_E v_T)}-1)
           - \frac{\beta_r}{1+\beta_r} I_{S,C} \cdot (e^{v_C/(\eta_C v_T)}-1)$

$i_C = -\frac{\beta_f}{1+\beta_f} I_{S,E} \cdot (e^{v_E/(\eta_E v_T)}-1)
           + I_{S,C} \cdot (e^{v_C/(\eta_C v_T)}-1).$

Additionally, terminal series resistances are supported.

The parameters are set using named arguments:

| parameter | description |
|:--------- |:------------|
| `typ` | Either `:npn` or `:pnp`, depending on desired transistor type
| `is`  | Reverse saturation current in Ampere
| `η`   | Emission coefficient
| `isc` | Collector reverse saturation current in Ampere (overriding `is`)
| `ise` | Emitter reverse saturation current in Ampere (overriding `is`)
| `ηc`  | Collector emission coefficient (overriding `η`)
| `ηe`  | Emitter emission coefficient (overriding `η`)
| `βf`  | Forward current gain
| `βr`  | Reverse current gain
| `ilc` | Base-collector junction leakage current in Ampere
| `ile` | Base-emitter junction leakage current in Ampere
| `ηcl` | Base-collector junction leakage emission coefficient (overriding `η`)
| `ηel` | Base-emitter junction leakage emission coefficient (overriding `η`)
| `vaf` | Forward Early voltage in Volt
| `var` | Reverse Early voltage in Volt
| `ikf` | Forward knee current (gain roll-off) in Ampere
| `ikr` | Reverse knee current (gain roll-off) in Ampere
| `re`  | Emitter terminal resistance
| `rc`  | Collector terminal resistance
| `rb`  | Base terminal resistance

Pins: `base`, `emitter`, `collector`
"""
function bjt(typ; is=1e-12, η=1, isc=is, ise=is, ηc=η, ηe=η, βf=1000, βr=10,
             ile=0, ilc=0, ηcl=ηc, ηel=ηe, vaf=Inf, var=Inf, ikf=Inf, ikr=Inf,
             re=0, rc=0, rb=0)
    local polarity
    if typ == :npn
        polarity = 1
    elseif typ == :pnp
        polarity = -1
    else
        throw(ArgumentError(string("Unknown bjt type ", typ,
                                   ", must be :npn or :pnp")))
    end
    kernel = quote
        i_f = $(βf/(1+βf)*ise) * (expE - 1)
        i_r = $(βr/(1+βr)*isc) * (expC - 1)
        di_f1 = $(βf/(1+βf)*ise/(25e-3*ηe)) * expE
        di_r2 = $(βr/(1+βr)*isc/(25e-3*ηc)) * expC
    end
    if var == Inf && vaf == Inf && ikf == Inf && ikr == Inf
        append!(kernel.args, (quote
            i_cc = i_f-i_r
            di_cc1 = di_f1
            di_cc2 = -di_r2
        end).args)
    elseif (var ≠ Inf || vaf ≠ Inf) && ikf == Inf && ikr == Inf
        append!(kernel.args, (quote
            # inverse Early voltage factor
            q₁⁻¹ = 1 - vE*$(1/var) - vC*$(1/vaf)
            i_cc = q₁⁻¹ * (i_f-i_r)
            # partial derivatives without high level injection effect
            dq₁⁻¹1 = $(-1/var)
            dq₁⁻¹2 = $(-1/vaf)
            di_cc1 = dq₁⁻¹1*(i_f-i_r) + q₁⁻¹*di_f1
            di_cc2 = dq₁⁻¹2*(i_f-i_r) - q₁⁻¹*di_r2
        end).args)
    elseif var == Inf && vaf == Inf && (ikf ≠ Inf || ikr ≠ Inf)
        append!(kernel.args, (quote
            # high level injection effect
            q₂ = i_f*$(1/ikf) + i_r*$(1/ikr)
            qden = 1+sqrt(1+4q₂)
            qfact = 2/qden
            i_cc = qfact * (i_f-i_r)
            # partial derivatives without Early effect
            dq₂1 = di_f1*$(1/ikf)
            dq₂2 = di_r2*$(1/ikr)
            dqfact1 = -4dq₂1/(qden-1) / (qden^2)
            dqfact2 = -4dq₂2/(qden-1) / (qden^2)
            di_cc1 = dqfact1*(i_f-i_r) + qfact*di_f1
            di_cc2 = dqfact2*(i_f-i_r) - qfact*di_r2
        end).args)
    else
        append!(kernel.args, (quote
            # inverse Early voltage factor
            q₁⁻¹ = 1 - vE*$(1/var) - vC*$(1/vaf)
            # high level injection effect
            q₂ = i_f*$(1/ikf) + i_r*$(1/ikr)
            qden = 1+sqrt(1+4q₂)
            qfact = 2q₁⁻¹/qden
            i_cc = qfact * (i_f-i_r)
            # partial derivatives with high level injection effect and Early effect
            dq₁⁻¹1 = $(-1/var)
            dq₁⁻¹2 = $(-1/vaf)
            dq₂1 = di_f1*$(1/ikf)
            dq₂2 = di_r2*$(1/ikr)
            dqfact1 = (2dq₁⁻¹1*qden - q₁⁻¹*4dq₂1/(qden-1)) / (qden^2)
            dqfact2 = (2dq₁⁻¹2*qden - q₁⁻¹*4dq₂2/(qden-1)) / (qden^2)
            di_cc1 = dqfact1*(i_f-i_r) + qfact*di_f1
            di_cc2 = dqfact2*(i_f-i_r) - qfact*di_r2
        end).args)
    end

    if ile ≠ 0
        if ηel ≠ ηe
            append!(kernel.args, (quote
                expEl = exp(vE*$(1/(25e-3*ηel)))
                iBE = $(1/βf)*i_f + $ile*(expEl - 1)
                diBE1 = $(1/βf)*di_f1 + $(ile/(25e-3*ηe))*expEl
            end).args)
        else
            append!(kernel.args, (quote
                iBE = $(1/βf)*i_f + $ile*(expE - 1)
                diBE1 = $(1/βf)*di_f1 + $(ile/(25e-3*ηe))*expE
            end).args)
        end
    else
        append!(kernel.args, (quote
            iBE = $(1/βf)*i_f
            diBE1 = $(1/βf)*di_f1
        end).args)
    end
    if ilc ≠ 0
        if ηcl ≠ ηc
            append!(kernel.args, (quote
                expCl = exp(vC*$(1/(25e-3*ηcl)))
                iBC = $(1/βr)*i_r + $ilc*(expCl - 1)
                diBC2 = $(1/βr)*di_r2 + $(ilc/(25e-3*ηc))*expCl
            end).args)
        else
            append!(kernel.args, (quote
                iBC = $(1/βr)*i_r + $ilc*(expC - 1)
                diBC2 = $(1/βr)*di_r2 + $(ilc/(25e-3*ηc))*expC
            end).args)
        end
    else
        append!(kernel.args, (quote
            iBC = $(1/βr)*i_r
            diBC2 = $(1/βr)*di_r2
        end).args)
    end

    nonlinear_eq =
        quote
            let vE = q[1], vC = q[2], iE = q[3], iC = q[4],
                expE=exp(vE*$(1/(25e-3*ηe))), expC=exp(vC*$(1/(25e-3*ηc)))

                $kernel

                res[1] = i_cc + iBE - iE
                res[2] = -i_cc + iBC - iC
                J[1,1] = di_cc1 + diBE1
                J[1,2] = di_cc2
                J[1,3] = -1.0
                J[1,4] = 0.0
                J[2,1] = -di_cc1
                J[2,2] = -di_cc2 + diBC2
                J[2,3] = 0.0
                J[2,4] = -1.0
            end
        end
    return Element(mv=[1 0; 0 1; 0 0; 0 0],
                   mi = [-(re+rb) -rb; -rb -(rc+rb); 1 0; 0 1],
                   mq = Matrix{Int}(-polarity*I, 4, 4), nonlinear_eq = nonlinear_eq,
                   pins = [:base; :emitter; :base; :collector])
end

"""
    opamp()

Creates an ideal operational amplifier. It enforces the voltage between the
input pins to be zero without sourcing any current while sourcing arbitrary
current on the output pins wihtout restricting their voltage.

Note that the opamp has two output pins, one of which will typically be
connected to a ground node and has to provide the current sourced on the other
output pin.

Pins: `in+` and `in-` for input, `out+` and `out-` for output
"""
opamp() = Element(mv=[0 0; 1 0], mi=[1 0; 0 0],
                  pins=["in+", "in-", "out+", "out-"])

doc"""
    opamp(Val{:macak}, gain, vomin, vomax)

Creates a clipping operational amplifier where input and output voltage are
related by

$v_\text{out} = \frac{1}{2}\cdot(v_\text{max}+v_\text{min})
                   +\frac{1}{2}\cdot(v_\text{max}-v_\text{min})\cdot
                    \tanh\left(\frac{g}{\frac{1}{2}\cdot(v_\text{max}-v_\text{min})}\cdot  v_\text{in}\right).$

The input current is zero, the output current is arbitrary.

Note that the opamp has two output pins, one of which will typically be
connected to a ground node and has to provide the current sourced on the other
output pin.

Pins: `in+` and `in-` for input, `out+` and `out-` for output
"""
function opamp(::Type{Val{:macak}}, gain, vomin, vomax)
    offset = 0.5 * (vomin + vomax)
    scale = 0.5 * (vomax - vomin)
    nonlinear_eq =
        quote
            let vi = q[1], vo = q[2], vi_scaled = vi * $(gain/scale)
                res[1] = tanh(vi_scaled) * $(scale) - vo
                J[1,1] = $(gain) / cosh(vi_scaled)^2
                J[1,2] = -1
            end
        end
    return Element(mv=[0 0; 1 0; 0 1], mi=[1 0; 0 0; 0 0], mq=[0 0; -1 0; 0 -1],
                   u0=[0; 0; offset],
                   nonlinear_eq = nonlinear_eq,
                   pins=["in+", "in-", "out+", "out-"])
end

@Base.deprecate(opamp_ideal, opamp)
@Base.deprecate(opamp_macak(gain, vomin, vomax), opamp(Val{:macak}, gain, vomin, vomax))
