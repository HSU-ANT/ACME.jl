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
mosfet
```

## Integrated Circuits

```@docs
opamp()
opamp(::Type{Val{:macak}}, gain, vomin, vomax)
```
# Adding custom elements

One advantage of ACME is that it is relatively easy to define one's own circuit
elements. However, the required formalism is a bit non-obvious. So this section
is meant to be an introduction into steps required to implement a custom circuit
element. As a running example, we shall consider an operational transconductance
amplifier (OTA). As the name suggests, it is similar to an operation amplifier
in that it has two input pins, acting as a differential input, and a single
output pin. For the OTA, the output current is (approximately) proportional to
the differential input voltage. The factor of proportionality (transconductance)
can be controlled by the current drawn at a further pin. Additionally, the OTA
has power supply pins, of course.

We shall start with a simple linear model with fixed transconductance. We
idealize the input current to be zero adn the output voltage as arbitrary, i.e.
the OTA acts as an ideal, voltage-controlled current source. The first step is
to define the ports (pairs of pins) of the element, such that the behavior of
the element can be defined in terms of the voltages across and currents through
these ports. For the OTA, the first port can obviously be chosen as the input
pins. For the second port, the output pin has to be paired with a pin that
provides the current to be sourced at the output to maintain a net current sum
of zero. As for the opamp models, we therefore introduce a negative output pin.
For the mathematical formulation, the port voltages and currents are collected
in vectors $\bm{v}$ and $\bm{i}$, respectively. Following the convention that
the input port is the first entry, the OTA is therefore described with
$i_2=-g\cdot v_1$ (where $g$ is the transconductance) and $i_1=0$. The negative
sign in the first equation is required because ACME considers the currents on
the element's inside, but here, we want a positive input voltage to yield a
positive output current on the outside.

ACME uses a matrix notation that in full generality looks like
```math
\bm{M}_\text{v}\bm{v} + \bm{M}_\text{i}\bm{i} + \bm{M}_\text{x}\bm{x}
+ \bm{M}_{\dot{\text{x}}}\dot{\bm{x}} + \bm{M}_\text{q}\bm{q}
= \bm{M}_\text{u}\bm{u} + \bm{u}_0.
```

The matrices $\bm{M}_{\cdot}$ describe the element and will have to be
determined in the following. The vector $\bm{x}$ holds the internal states of
the element and $\dot{\bm{x}}$ their derivatives. We model a stateless OTA, so
these parts can be ignored by choosing $\bm{M}_\text{x}$ and
$\bm{M}_{\dot{\text{x}}}$ as zero-column matrices. The vector $\bm{q}$ becomes
relevant for nonlinear elements, but can likewise be ignored for the moment
choosing $\bm{M}_\text{q}$ as a zero-column matrix. The vector $\bm{u}$ contains
externally controlled inputs, which we don't have for the OTA, so
$\bm{M}_\text{u}$ is also chosen to have zero columns. Finally, $\bm{u}_0$ is
again used for describing the element. For the linear OTA, we can write the two
equations as
```math
\begin{pmatrix} g & 0 \\ 0 & 0 \end{pmatrix}\bm{v}
+ \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}\bm{i}
= \begin{pmatrix} 0 \\ 0 \end{pmatrix}
```
where the first row translates to $i_2=-g\cdot v_1$ and the second row to $i_1=0$.
Thus, we obtain
```math
\bm{M}_\text{v}=\begin{pmatrix} g & 0 \\ 0 & 0 \end{pmatrix}
\quad \bm{M}_\text{i} = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}
\quad \bm{u}_0 = \begin{pmatrix} 0 \\ 0 \end{pmatrix}.
```

Translation to ACME is straight-forward, as the `Element` constructor takes
self-explanatorily named keyword arguments:
```julia
g=10e-3
ACME.Element(mv=[g 0; 0 0], mi=[0 1; 1 0],  ports=["in+" => "in-", "out+" => "out-"])
```
Omitted matrices are taken to be empty/all-zero matrices of suitable size.
The ports are specified as pairs of pins. The names can be chosen arbitrarily.
If the same pin name is used in multiple ports, these ports share a pin.
For example, the [`bjt`](@ref) has a base-collector and a base-emitter port,
sharing base as a common pin.

It is usually convenient to wrap the element creation in a suitable function
like it is done for the elements defined as part of ACME:
```jldoctest ota_linear; output=false, setup=:(using ACME)
ota(g) = ACME.Element(mv=[g 0; 0 0], mi=[0 1; 1 0],  ports=["in+" => "in-", "out+" => "out-"])

# output

ota (generic function with 1 method)
```
A simple test circuit for this OTA could then be run with
```jldoctest ota_linear
circ = @circuit begin
    Jin = voltagesource()
    U = ota(10e-3), ["in+"] ⟷ Jin[+], ["in-"] ⟷ Jin[-]
    Jout = currentprobe(), [+] ⟷ U["out+"], [-] ⟷ U["out-"]
end
model = DiscreteModel(circ, 1/44100)
run!(model, [-2.0 -1.0 0.0 1.0 2.0])

# output

1×5 Matrix{Float64}:
 -0.02  -0.01  0.0  0.01  0.02
```

Real OTAs exhibit nonlinear behavior. A full treatment of the all the sources of
nonlinearities is beyond this manual, but a first refinement is given by the
input/output relationship $i_2=-i_\text{bias}\cdot\tanh(v_1/(2v_\text{T}))$,
where $v_\text{T}$ is the thermal voltage (approximately 25 mV at room
temperature) and the transconductance is controlled with the bias current
$i_\text{bias}$ by $g=i_\text{bias}/(2v_\text{T})$ for small input voltages. For
higher input voltages, saturation occurs and the output current cannot exceed
the bias current.

To encode nonlinear behavior with ACME, all quantities that participate in the
nonlinearity have to be collected in the $\bm{q}$ vector, i.e. we now need a
non-empty matrix $\bm{M}_\text{q}$, and the nonlinear equation has be provided
via a function $\bm{f}$ such that $\bm{f}(\bm{q})=\bm{0}$ is the desired
condition. For the present example, we choose $q_1=v_1$ and $q_2=i_2$ with
```math
\begin{pmatrix} 1 & 0 \\ 0 & 0 \\ 0 & 0 \end{pmatrix}\bm{v}
+ \begin{pmatrix} 0 & 0 \\ 0 & 1 \\ 1 & 0 \end{pmatrix}\bm{i}
+ \begin{pmatrix} -1 & 0 \\ 0 & -1 \\ 0 & 0 \end{pmatrix}\bm{q}
= \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}
```
i.e.
```math
\bm{M}_\text{v}=\begin{pmatrix} 1 & 0 \\ 0 & 0 \\ 0 & 0 \end{pmatrix}
\quad \bm{M}_\text{i} = \begin{pmatrix} 0 & 0 \\ 0 & 1 \\ 1 & 0 \end{pmatrix}
\quad \bm{M}_\text{q}=\begin{pmatrix} -1 & 0 \\ 0 & -1 \\ 0 & 0 \end{pmatrix}
\quad \bm{u}_0 = \begin{pmatrix} 0 \\ 0 \\ 0 \end{pmatrix}
```
and
```math
\bm{f}(\bm{q}) = \begin{pmatrix} i_\text{bias}\cdot\tanh(q_1/(2v_\text{T})) + q_2 \end{pmatrix}.
```
In order to do nonlinear equation solving, ACME requires not only the function
$\bm{f}$, but also its Jacobian, in this case
```math
\bm{J}(\bm{q}) = \begin{pmatrix} \frac{df_1}{dq_1} & \frac{df_1}{dq_2} \end{pmatrix}
= \begin{pmatrix} \frac{i_\text{bias}}{2v_\text{T}\cdot\cosh^2(q_1/(2v_\text{T}))} & 1 \end{pmatrix}.
```
The final bit of information is that the function value (residual) and the value
of the Jacobian have to returned as an `@SVector` and `@SMatrix` (from the
StaticArrays package), respectively, combined in a tuple. Thus, the implementation becomes:

```@example ota_nonlinear
using ACME, StaticArrays # hide
ota(i_bias) = ACME.Element(
    mv=[1 0; 0 0; 0 0], mi=[0 0; 0 1; 1 0], mq=[-1 0; 0 -1; 0 0],
    ports=["in+" => "in-", "out+" => "out-"],
    nonlinear_eq = function (q)
        res = @SVector [i_bias * tanh(q[1]/50e-3) + q[2]]
        J = @SMatrix [i_bias / (50e-3*coth(q[1]/50e-3)^2) 1]
        return (res, J)
    end
)
nothing # hide
```
With the same test circuit as before, we can create a plot of the input/output relationship:
```@example ota_nonlinear
circ = @circuit begin
    Jin = voltagesource()
    U = ota(10e-3), ["in+"] ⟷ Jin[+], ["in-"] ⟷ Jin[-]
    Jout = currentprobe(), [+] ⟷ U["out+"], [-] ⟷ U["out-"]
end
model = DiscreteModel(circ, 1/44100)
u = range(-0.2, 0.2, length=200)
y = run!(model, u')'
using Plots
plot(u, y; xlabel="input voltage", ylabel="output current", legend=false)
```
For the small input voltages, the theoretical transconductance of 10 mA / 50 mV
= 0.2 Ω⁻¹ is obtained, while for the larger voltages, the saturation effect
becomes clearly visible.