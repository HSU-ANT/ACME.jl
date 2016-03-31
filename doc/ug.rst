************
 User Guide
************

Element Creation
================

All circuit elements are created by calling corresponding functions; see the
:ref:`element-reference` for details.

Circuit Description
===================

Circuits are described using ``Circuit`` instances, created with ``Circuit()``.
Once a ``Circuit`` and elements have been created, the elements can be added to
the circuit using the ``add!`` method::

    r = resistor(1e3)
    c = capacitor(22e-9)
    circ = Circuit()
    add!(circ, r)
    add!(circ, c)

Multiple elements can be added also be at once; the last two lines could have
been replaced with ``add!(circ, r, c)``.

In many cases, however, explicitly calling ``add!`` is not necessary. All that
is needed is ``connect!``, which connects two (or more) element pins. The
elements to which these pins belong are automatically added to the circuit if
needed. The only reason to explicitly call ``add!`` is to control the insertion
order of sources or sinks, which determines the order in which inputs have to be
provided and outputs are obtained.

Pins are obtained from elements using ``[]``-style indexing, i.e. ``r[1]`` gives
the first pin of the resistor defined above. So this connects the first pin of
the resistor with the first pin of the capacitor::

    connect!(circ, r[1], c[1])

Further connections involving the same pins are possible and will *not* replace
existing ones. So this will effectively shorten the resistor, because now both
of its pins are connected to ``c[1]``::

    connect!(circ, r[2], c[1])

Note that not all elements have numbered pins. For elements with polarity, they
may be called ``+`` and ``-``, while a bipolar transistor has pins ``base``,
``collector``, and ``emitter``. The pins provided by each type of element are
described in the :ref:`element-reference`. Internally, the pin designators are
``Symbol``\s. However, not all symbols are conveniently entered in Julia:
``:base`` is nice, ``symbol("1")`` less so. Therefore, the ``[]`` operation on
elements also accepts integers and strings and converts them to the respective
``Symbol``\s. So ``r[symbol("1")]`` is equivalent to ``r[1]`` and (assuming
``d`` to be a diode) ``d[:+]`` is equivalent to ``d["+"]`` (but ``d[+]`` does
not work).

In addition to pins, ``connect!`` also accepts ``Symbol``\s as input. This
creates named nets which may improve readability for nets with many conneted
pins::

    connect!(c[2], :gnd)
    connect!(r[2], :gnd)

Again, this only adds connections, keeping existing ones, so together with the
above snippets, now all pins are connected to each other and to net named
``gnd``. It is even possible to connect multple named nets to each other, though
this will only rarely be useful.

Model Creation and Use
======================

A ``Circuit`` only stores elements and information about their connections. To
simulate a circuit, a model has to be derived from it. This can be as simple
as::

    model = DiscreteModel(circ, 1/44100)

Here, ``1/44100`` denotes the sampling interval, i.e. the reciprocal of the
sampling rate at which the model should run. Optionally, one can specify the
solver to use for solving the model's non-linear equation as a type parameter::

    model = DiscreteModel{HomotopySolver{SimpleSolver}}(circ, 1/44100)

See Solvers_ for more information about the available solvers.

Once a model is created, it can be run::

    y = run!(model, u)

The input ``u`` is matrix with one row for each of the circuit's inputs and one
column for each time step to simulate. Likewise, the output ``y`` will be a
matrix with one row for each of the circuit's outputs and one column for each
simulated time step. The order of the rows will correspond to the order in which
the respective input and output elements were added to the ``Circuit``. To
simulate a circuit without inputs, a matrix with zero rows may be passed::

    y = run!(model, zeros(0, 100))

The internal state of the model (e.g. capacitor charges) is preserved accross
calls to ``run!``. Initially, all states are zeroed. It is also possible to set
the states to a steady state (if one can be found) with::

    steadystate!(model)

This is often desirable for circuits were bias voltages are only slowly obtained
after turning them on.

Solvers
=======

``SimpleSolver``
----------------

The ``SimpleSolver`` is the simplest available solver. It uses Newton iteration
which features fast local convergence, but makes no guarantees about global
convergence. The initial solution of the iteration is obtained by extrapolating
the last solution found (or another solution provided externally) using the
available Jacobians. Due to the missing global convergence, the ``SimpleSolver``
is rarely useful as such.

``HomotopySolver``
------------------

The ``HomotopySolver`` extends an existing solver (provided as a type parameter)
by applying homotopy to (at least theoretically) ensure global convergence. It
can be combined with the ``SimpleSolver`` as ``HomotopySolver{SimpleSolver}`` to
obtain a useful Newton homtopy solver with generally good convergence
properties.

``CachingSolver``
-----------------

The ``CachingSolver`` extends an existing solver (provided as a type parameter)
by storing found solutions in a k-d tree to use as initial solutions in the
future. Whenever the underlying solver needs more than a preset number of
iterations (defaults to five), the solution will be stored. Storing new
solutions is a relatively expensive operation, so until the stored solutions
suffice to ensure convergence in few iterations throughout, use of a
``CachingSolver`` may actually slow things down.

The default solver used is a ``HomotopySolver{CachingSolver{SimpleSolver}}``.
