# Copyright 2015, 2016, 2017, 2018, 2019, 2020 Martin Holters
# See accompanying license file.

export resistor, potentiometer, capacitor, inductor, transformer,
       voltagesource, currentsource,
       voltageprobe, currentprobe, diode, bjt, mosfet, opamp


"""
    resistor(r)

Creates a resistor obeying Ohm’s law. The resistance `r` has to be given in Ohm.

Pins: `1`, `2`
"""
resistor(r) = Element(mv=-1, mi=r)

potentiometer(r, pos) = Element(mv=Matrix{Int}(-I, 2, 2), mi=[r*pos 0; 0 r*(1-pos)],
                                ports=[1 => 2, 2 => 3])
potentiometer(r) =
    Element(mv=[Matrix{Int}(I, 2, 2); zeros(3, 2)],
            mi=[zeros(2, 2); Matrix{Int}(I, 2, 2); zeros(1, 2)],
            mq=Matrix{Int}(-I, 5, 5), mu=[zeros(4, 1); -1],
            nonlinear_eq =
                (@inline function (q)
                    v1, v2, i1, i2, pos=q
                    res = @SVector [v1 - r*pos*i1, v2 - r*(1-pos)*i2]
                    J = @SMatrix [1 0 -r*pos 0 -r*i1; 0 1 0 -r*(1-pos) -r*i2]
                    return (res, J)
                end),
            ports=[1 => 2, 2 => 3])

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
            ports=[:primary1 => :primary2, :secondary1 => :secondary2])

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
    nonlinear_eq = @inline function (q)
        coth_q1 = coth(q[1])
        a_q1 = abs(q[1])
        L_q1 = a_q1 < 1e-4 ? q[1]/3 : coth_q1 - 1/q[1]
        Ld_q1 = a_q1 < 1e-4 ? 1/3 : 1/q[1]^2 - coth_q1^2 + 1
        Ld2_q1 = a_q1 < 1e-3 ? -2/15*q[1] : 2*coth_q1*(coth_q1^2 - 1) - 2/q[1]^3
        δ = q[3] > 0 ? 1.0 : -1.0

        Man = Ms*L_q1
        δM = sign(q[3]) == sign(Man - q[2]) ? 1.0 : 0.0

        den = δ*(k*(1-c))-α*(Man-q[2])
        # at present, the error needs to be scaled to be comparable to those of
        # the other elements, hence the factor 1e-4/Ms
        res = @SVector [(1e-4/Ms) * ((1-c) * δM*(Man-q[2])/den * q[3]
                                     + (c*Ms/a)*(q[3]+α*q[4])*Ld_q1 - q[4])]
        J_1_1 = (1e-4/Ms) * (((1-c)^2*k*Ms) * δM*Ld_q1*δ/den^2 * q[3]
                               + (c*Ms/a)*(q[3]+α*q[4])*Ld2_q1)
        J_1_2 = (1e-4/Ms) * -(1-c)^2*k * δM*δ/den^2 * q[3]
        J_1_3 = (1e-4/Ms) * ((1-c) * δM*(Man-q[2])/den + (c*Ms/a)*Ld_q1)
        J_1_4 = (1e-4/Ms) * ((c * Ms/a * α)*Ld_q1 - 1)
        return (res, @SMatrix [J_1_1 J_1_2; J_2_1 J_2_2])
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
voltagesource(v; rs=0) = Element(mv=1, mi=-rs, u0=v, ports=[:+ => :-])
voltagesource(; rs=0) = Element(mv=1, mi=-rs, mu=1, ports=[:+ => :-])

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
currentsource(i; gp=0) = Element(mv=gp, mi=-1, u0=i, ports=[:+ => :-])
currentsource(; gp=0) = Element(mv=gp, mi=-1, mu=1, ports=[:+ => :-])

"""
    voltageprobe()

Creates a voltage probe, providing the measured voltage as a circuit output.
Optionally, an internal parallel conductance `gp` (in Ohm⁻¹) can be given which
defaults to zero.

Pins: `+` and `-` with the output voltage being measured from `+` to `-`
"""
voltageprobe(;gp=0) = Element(mv=-gp, mi=1, pv=1, ports=[:+ => :-])

"""
    currentprobe()

Creates a current probe, providing the measured current as a circuit output.
Optionally, an internal series resistance `rs` (in Ohm) can be given which
defaults to zero.

Pins: `+` and `-` with the output current being the current entering the probe
at `+`
"""
currentprobe(;rs=0) = Element(mv=1, mi=-rs, pi=1, ports=[:+ => :-])

@doc raw"""
    diode(;is=1e-12, η = 1)

Creates a diode obeying Shockley's law
$i=I_S\cdot(e^{v/(\eta v_T)}-1)$ where $v_T$ is fixed at 25 mV.
The reverse saturation current `is` has to be given in Ampere, the emission
coefficient `η` is unitless.

Pins: `+` (anode) and `-` (cathode)
""" diode(;is::Real=1e-12, η::Real = 1) =
  Element(mv=[1;0], mi=[0;1], mq=[-1 0; 0 -1], ports=[:+ => :-], nonlinear_eq =
        @inline function(q)
            v, i = q
            ex = exp(v*(1 / (25e-3 * η)))
            res = @SVector [is * (ex - 1) - i]
            J = @SMatrix [is/(25e-3 * η) * ex -1]
            return (res, J)
        end
  )

@doc raw"""
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
""" function bjt(typ; is=1e-12, η=1, isc=is, ise=is, ηc=η, ηe=η, βf=1000, βr=10,
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

    nonlinear_eq =
        @inline function (q)
            vE, vC, iE, iC = q
            expE=exp(vE*(1/(25e-3*ηe)))
            expC=exp(vC*(1/(25e-3*ηc)))
            i_f = (βf/(1+βf)*ise) * (expE - 1)
            i_r = (βr/(1+βr)*isc) * (expC - 1)
            di_f1 = (βf/(1+βf)*ise/(25e-3*ηe)) * expE
            di_r2 = (βr/(1+βr)*isc/(25e-3*ηc)) * expC
            if var == Inf && vaf == Inf && ikf == Inf && ikr == Inf
                i_cc = i_f-i_r
                di_cc1 = di_f1
                di_cc2 = -di_r2
            elseif (var ≠ Inf || vaf ≠ Inf) && ikf == Inf && ikr == Inf
                # inverse Early voltage factor
                q₁⁻¹ = 1 - vE*(1/var) - vC*(1/vaf)
                i_cc = q₁⁻¹ * (i_f-i_r)
                # partial derivatives without high level injection effect
                dq₁⁻¹1 = (-1/var)
                dq₁⁻¹2 = (-1/vaf)
                di_cc1 = dq₁⁻¹1*(i_f-i_r) + q₁⁻¹*di_f1
                di_cc2 = dq₁⁻¹2*(i_f-i_r) - q₁⁻¹*di_r2
            elseif var == Inf && vaf == Inf && (ikf ≠ Inf || ikr ≠ Inf)
                # high level injection effect
                q₂ = i_f*(1/ikf) + i_r*(1/ikr)
                qden = 1+sqrt(1+4q₂)
                qfact = 2/qden
                i_cc = qfact * (i_f-i_r)
                # partial derivatives without Early effect
                dq₂1 = di_f1*(1/ikf)
                dq₂2 = di_r2*(1/ikr)
                dqfact1 = -4dq₂1/(qden-1) / (qden^2)
                dqfact2 = -4dq₂2/(qden-1) / (qden^2)
                di_cc1 = dqfact1*(i_f-i_r) + qfact*di_f1
                di_cc2 = dqfact2*(i_f-i_r) - qfact*di_r2
            else
                # inverse Early voltage factor
                q₁⁻¹ = 1 - vE*(1/var) - vC*(1/vaf)
                # high level injection effect
                q₂ = i_f*(1/ikf) + i_r*(1/ikr)
                qden = 1+sqrt(1+4q₂)
                qfact = 2q₁⁻¹/qden
                i_cc = qfact * (i_f-i_r)
                # partial derivatives with high level injection effect and Early effect
                dq₁⁻¹1 = -1/var
                dq₁⁻¹2 = -1/vaf
                dq₂1 = di_f1*(1/ikf)
                dq₂2 = di_r2*(1/ikr)
                dqfact1 = (2dq₁⁻¹1*qden - q₁⁻¹*4dq₂1/(qden-1)) / (qden^2)
                dqfact2 = (2dq₁⁻¹2*qden - q₁⁻¹*4dq₂2/(qden-1)) / (qden^2)
                di_cc1 = dqfact1*(i_f-i_r) + qfact*di_f1
                di_cc2 = dqfact2*(i_f-i_r) - qfact*di_r2
            end
            iBE = (1/βf)*i_f
            diBE1 = (1/βf)*di_f1
            if ile ≠ 0
                if ηel ≠ ηe
                    expEl = exp(vE*(1/(25e-3*ηel)))
                else
                    expEl = expE
                end
                iBE += ile*(expEl - 1)
                diBE1 += (ile/(25e-3*ηe))*expEl
            end
            iBC = (1/βr)*i_r
            diBC2 = (1/βr)*di_r2
            if ilc ≠ 0
                if ηcl ≠ ηc
                    expCl = exp(vC*(1/(25e-3*ηcl)))
                else
                    expCl = expC
                end
                iBC += ilc*(expCl - 1)
                diBC2 += (ilc/(25e-3*ηc))*expCl
            end
            res = @SVector [i_cc + iBE - iE, -i_cc + iBC - iC]
            J = @SMatrix [di_cc1 + diBE1 di_cc2          -1.0 0.0;
                          -di_cc1        -di_cc2 + diBC2 0.0  -1.0]
            return (res, J)
        end
    return Element(mv=[1 0; 0 1; 0 0; 0 0],
                   mi = [-(re+rb) -rb; -rb -(rc+rb); 1 0; 0 1],
                   mq = Matrix{Int}(-polarity*I, 4, 4), nonlinear_eq = nonlinear_eq,
                   ports = [:base => :emitter, :base => :collector])
end

@doc raw"""
    mosfet(typ; vt=0.7, α=2e-5, λ=0)

Creates a MOSFET transistor with the simple model

$i_D=\begin{cases}
  0 & \text{if } v_{GS} \le v_T \\
  \alpha \cdot (v_{GS} - v_T - \tfrac{1}{2}v_{DS})\cdot v_{DS}
  \cdot (1 + \lambda v_{DS})
  & \text{if } v_{DS} \le v_{GS} - v_T \cap v_{GS} > v_T \\
  \frac{\alpha}{2} \cdot (v_{GS} - v_T)^2 \cdot (1 + \lambda v_{DS})
  & \text{otherwise.}
\end{cases}$

The `typ` parameter chooses between NMOS (`:n`) and PMOS (`:p`). The threshold
voltage `vt` is given in Volt, `α` (in A/V²) is a constant depending on the
physics and dimensions of the device, and `λ` (in V⁻¹) controls the channel
length modulation.

Optionally, it is possible to specify tuples of coefficients for `vt` and `α`.
These will be used as polynomials in $v_{GS}$ to determine $v_T$ and $\alpha$,
respectively. E.g. with `vt=(0.7, 0.1, 0.02)`, the $v_{GS}$-dpendent threshold
voltage $v_T = 0.7 + 0.1\cdot v_{GS} + 0.02\cdot v_{GS}^2$ will be used.

Pins: `gate`, `source`, `drain`
""" function mosfet(typ; vt=0.7, α=2e-5, λ=0)
    if typ == :n
        polarity = 1
    elseif typ == :p
        polarity = -1
    else
        throw(ArgumentError("Unknown mosfet type $(typ), must be :n or :p"))
    end
    vt = (vt...,)
    α = (α...,)
    dvt = vt[2:end] .* (1:length(vt)-1...,)
    dα = α[2:end] .* (1:length(α)-1...,)
    let polarity = polarity, α = α, vt = vt
        return Element(mv=[-1 0; 0 -1; 0 0; 0 0],
            mi=[0 0; 0 0; 0 -1; 1 0],
            mq=polarity*[1 0 0; 0 1 0; 0 0 1; 0 0 0],
            ports=[:gate => :source, :drain => :source],
            nonlinear_eq = @inline function (q)
                vgs, vds, id = q
                α´ = evalpoly(polarity*vgs, α)
                if !isempty(dα)
                    dα´_dvgs = evalpoly(polarity*vgs, dα)
                else
                    dα´_dvgs = 0
                end
                vt´ = evalpoly(polarity*vgs, vt)
                if !isempty(dvt)
                    dvt´_dvgs = evalpoly(polarity*vgs, dvt)
                else
                    dvt´_dvgs = 0
                end
                λ´ = vds ≥ 0 ? λ : zero(λ)
                if vgs <= vt´
                    res = @SVector [-id]
                    J = @SMatrix [0.0 0.0 -1.0]
                elseif vds <= vgs - vt´ # && vgs > vt´
                    res = @SVector [α´ * (vgs-vt´-0.5*vds)*vds*(1+λ´*vds) - id]
                    J = @SMatrix [α´*(1-dvt´_dvgs)*vds*(1+λ´*vds) + dα´_dvgs * (vgs-vt´-0.5*vds)*vds*(1+λ´*vds)  α´*(vgs-vt´+vds*(2*λ´*(vgs-vt´-0.75*vds)-1))  -1.0]
                else # 0 < vgs - vt´ < vds
                    res = @SVector [(α´/2) * (vgs-vt´)^2*(1+λ´*vds) - id]
                    J = @SMatrix [α´*(vgs-vt´)*(1-dvt´_dvgs)*(1+λ´*vds) + dα´_dvgs/2 * (vgs-vt´)^2*(1+λ´*vds) λ´*α´/2*(vgs-vt´)^2 -1.0]
                end
                return (res, J)
            end)
    end
end

@doc raw"""
    opamp(;maxgain=Inf, gain_bw_prod=Inf)

Creates a linear operational amplifier as a voltage-controlled voltage source.
The input current is zero while the input voltage is mapped to the output
voltage according to the transfer function

$H(f) = \frac{A_\text{max}}{\sqrt{A_\text{max}^2-1} i \frac{f}{f_\text{UG}} + 1}$

where $f$ is the signal frequency, $A_\text{max}$ (`maxgain`) is the maximum
open loop gain and $f_\text{UG}$ (`gain_bw_prod`) is the gain/bandwidth
product (unity gain bandwidth). For `gain_bw_prod=Inf` (the default), this
corresponds to a frequency-independent gain of `maxgain`. For `maxgain=Inf`
(the default), the amplifier behaves as a perfect integrator.

For both `maxgain=Inf` and `gain_bw_prod=Inf`, i.e. just `opamp()`, an ideal
operational amplifier is obtained that enforces the voltage between the input
pins to be zero while sourcing arbitrary current on the output pins without
restricting their voltage.

Note that the opamp has two output pins, where the negative one will typically
be connected to a ground node and has to provide the current sourced on the
positive one.

Pins: `in+` and `in-` for input, `out+` and `out-` for output
""" function opamp(;maxgain=Inf, gain_bw_prod=Inf)
    if gain_bw_prod==Inf # special case to avoid unnecessary state
        Element(mv=[0 0; 1 -1/maxgain], mi=[1 0; 0 0],
                ports=["in+" => "in-", "out+" => "out-"])
    else
        Element(mv=[0 0; -1/sqrt(1-1/maxgain^2) 0; 0 -1], mi=[1 0; 0 0; 0 0],
                mx=[0; 1/sqrt(maxgain^2-1); 1], mxd=[0; 1/(2π*gain_bw_prod); 0],
                ports=["in+" => "in-", "out+" => "out-"])
    end
end

@doc raw"""
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
""" function opamp(::Type{Val{:macak}}, gain, vomin, vomax)
    offset = 0.5 * (vomin + vomax)
    scale = 0.5 * (vomax - vomin)
    nonlinear_eq =
        @inline function (q)
            vi, vo = q
            vi_scaled = vi * (gain/scale)
            res = @SVector [tanh(vi_scaled) * scale - vo]
            J = @SMatrix [gain / cosh(vi_scaled)^2 -1.0]
            return (res, J)
        end
    return Element(mv=[0 0; 1 0; 0 1], mi=[1 0; 0 0; 0 0], mq=[0 0; -1 0; 0 -1],
                   u0=[0; 0; offset],
                   nonlinear_eq = nonlinear_eq,
                   ports=["in+" => "in-", "out+" => "out-"])
end
