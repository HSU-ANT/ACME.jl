.. _element-reference:

*******************
 Element Reference
*******************

Passives
========

.. function:: resistor(r)

    Creates a resistor obeying Ohm's law.

    :param r: Resistance in Ohm

    Pins: ``1``, ``2``

.. function:: capacitor(c)

    Creates a capacitor.

    :param c: Capacitance in Farad

    Pins: ``1``, ``2``

.. function:: inductor(l)

    Creates an inductor.

    :param l: Inductance in Henri

    Pins: ``1``, ``2``

.. function:: inductor(Val{:JA}; D, A, n, a, α, c, k, Ms)

    Creates a non-linear inductor based on the Jiles-Atherton model of
    magnetization assuming a toroidal core thin compared to its diameter.

    :param D: Torus diameter (in meters)
    :param A: Torus cross-sectional area (in square-meters)
    :param n: Winding's number of turns
    :param a: Shape parameter of the anhysteretic magnetization curve (in
              Ampere-per-meter)
    :param α: Inter-domain coupling
    :param c: Ratio of the initial normal to the initial anhysteretic
              differential susceptibility
    :param k: amount of hysteresis (in Ampere-per-meter)
    :param Ms: saturation magnetization (in Ampere-per-meter)

    A detailed discussion of the paramters can be found in D. C. Jiles and D. L.
    Atherton, “Theory of ferromagnetic hysteresis,” J. Magn. Magn. Mater., vol.
    61, no. 1–2, pp. 48–60, Sep. 1986 and J. H. B. Deane, “Modeling the dynamics
    of nonlinear inductor circuits,” IEEE Trans. Magn., vol. 30, no. 5, pp.
    2795–2801, 1994.

    Pins: ``1``, ``2``

.. function:: transformer(l1, l2; [coupling_coefficient=1,] [mutual_coupling=coupling_coefficient*sqrt(l1*l2)])

    Creates a transformer with two windings having inductances.

    :param l1: Primary self-inductance in Henri
    :param l2: Secondary self-inductance in Henri
    :param coupling_coefficient: Coupling coefficient (0 is not coupled, 1 is
                                 closely coupled)
    :param mutual_coupling: Mutual inductance in Henri; overrides
                            ``coupling_coefficient`` if both are given

    Pins: ``1`` and ``2`` for primary winding, ``3`` and ``4`` for secondary
    winding

.. function:: transformer(Val{:JA}; D, A, ns, a, α, c, k, Ms)

    Creates a non-linear transformer based on the Jiles-Atherton model of
    magnetization assuming a toroidal core thin compared to its diameter.

    :param D: Torus diameter (in meters)
    :param A: Torus cross-sectional area (in square-meters)
    :param ns: Windings' number of turns as a vector with one entry per winding
    :param a: Shape parameter of the anhysteretic magnetization curve (in
              Ampere-per-meter)
    :param α: Inter-domain coupling
    :param c: Ratio of the initial normal to the initial anhysteretic
              differential susceptibility
    :param k: amount of hysteresis (in Ampere-per-meter)
    :param Ms: saturation magnetization (in Ampere-per-meter)

    A detailed discussion of the paramters can be found in D. C. Jiles and D. L.
    Atherton, “Theory of ferromagnetic hysteresis,” J. Magn. Magn. Mater., vol.
    61, no. 1–2, pp. 48–60, Sep. 1986 and J. H. B. Deane, “Modeling the dynamics
    of nonlinear inductor circuits,” IEEE Trans. Magn., vol. 30, no. 5, pp.
    2795–2801, 1994.

    Pins: ``1`` and ``2`` for primary winding, ``3`` and ``4`` for secondary
    winding, and so on

Independent Sources
===================

.. function:: voltagesource([v])

    Creates a voltage source.

    :param v: Source voltage in Volt. If omitted, the source voltage will be an
              input of the circuit.

    Pins: ``+`` and ``-`` with ``v`` being measured from ``+`` to ``-``

.. function:: currentsource([i])

    Creates a current source.

    :param i: Source current in Ampere. If omitted, the source current will be an
              input of the circuit.

    Pins: ``+`` and ``-`` where ``i`` measures the current leaving source at the
    ``+`` pin

Probes
======

.. function:: voltageprobe()

    Creates a voltage probe, provding the measured voltage as a circuit output.

    Pins: ``+`` and ``-`` with the output voltage being measured from ``+`` to
    ``-``

.. function:: currentprobe()

    Creates a current probe, provding the measured current as a circuit output.

    Pins: ``+`` and ``-`` with the output current being the current entering the
    probe at ``+``

Semiconductors
==============

.. function:: diode(;[is=1e-12,] [η = 1])

    Creates a diode obeying Shockley's law
    :math:`i=I_S\cdot(e^{v/(\eta v_T)}-1)` where :math:`v_T` is fixed at 25 mV.

    :param is: Reverse saturation current in Ampere
    :param η: Emission coefficient

.. function:: bjt(typ; is=1e-12, η=1, isc=is, ise=is, ηc=η, ηe=η, βf=1000, βr=10)

    Creates a bipolar junction transistor obeying the Ebers-Moll equation

    .. math::
        i_E &= I_{S,E} \cdot (e^{v_E/(\eta_E v_T)}-1)
               - \frac{\beta_r}{1+\beta_r} I_{S,C} \cdot (e^{v_C/(\eta_C v_T)}-1)
        \\
        i_C &= -\frac{\beta_f}{1+\beta_f} I_{S,E} \cdot (e^{v_E/(\eta_E v_T)}-1)
               + I_{S,C} \cdot (e^{v_C/(\eta_C v_T)}-1)

    where :math:`v_T` is fixed at 25 mV.

    :param typ: Either ``:npn`` or ``:pnp``, depending on desired transistor type
    :param is: Reverse saturation current in Ampere
    :param η: Emission coefficient
    :param isc: Collector reverse saturation current in Ampere (overriding ``is``)
    :param ise: Emitter reverse saturation current in Ampere (overriding ``is``)
    :param ηc: Collector emission coefficient (overriding ``η``)
    :param ηe: Emitter emission coefficient (overriding ``η``)
    :param βf: Forward current gain
    :param βr: Reverse current gain

Integrated Circuits
===================

.. function:: opamp()

    Creates an ideal operational amplifier. It enforces the voltage between the
    input pins to be zero without sourcing any current while sourcing arbitrary
    current on the output pins wihtout restricting their voltage.

    Note that the opamp has two output pins, one of which will typically be
    connected to a ground node and has to provide the current sourced on the
    other output pin.

    Pins: ``in+`` and ``in-`` for input, ``out+`` and ``out-`` for output

.. function:: opamp(Val{:macak}, gain, vomin, vomax)

    Creates a clipping operational amplifier where input and output voltage are
    related by

    .. math::
        v_\text{out} = \frac{1}{2}\cdot(v_\text{max}+v_\text{min})
                       +\frac{1}{2}\cdot(v_\text{max}-v_\text{min})\cdot
                        \tanh\left(\frac{g}{\frac{1}{2}\cdot(v_\text{max}-v_\text{min})}\cdot  v_\text{in}\right).

    The input current is zero, the output current is arbitrary.

    Note that the opamp has two output pins, one of which will typically be
    connected to a ground node and has to provide the current sourced on the
    other output pin.

    Pins: ``in+`` and ``in-`` for input, ``out+`` and ``out-`` for output
