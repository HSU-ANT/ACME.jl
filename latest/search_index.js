var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#ACME.jl-Analog-Circuit-Modeling-and-Emulation-for-Julia-1",
    "page": "Home",
    "title": "ACME.jl - Analog Circuit Modeling and Emulation for Julia",
    "category": "section",
    "text": "ACME is a Julia package for the simulation of electrical circuits, focusing on audio effect circuits. It allows to programmatically describe a circuit in terms of elements and connections between them and then automatically derive a model for the circuit. The model can then be run on varying input data.ACME is based on the method described in M. Holters, U. Zölzer, \"A Generalized Method for the Derivation of Non-Linear State-Space Models from Circuit Schematics\"."
},

{
    "location": "gettingstarted.html#",
    "page": "Getting Started",
    "title": "Getting Started",
    "category": "page",
    "text": ""
},

{
    "location": "gettingstarted.html#Getting-Started-1",
    "page": "Getting Started",
    "title": "Getting Started",
    "category": "section",
    "text": ""
},

{
    "location": "gettingstarted.html#Installation-1",
    "page": "Getting Started",
    "title": "Installation",
    "category": "section",
    "text": "If you have not done so already, download and install Julia. (Any version starting with 0.4 should be fine; earlier ACME versions also support Julia 0.3.)To install ACME, start Julia and run:Pkg.add(\"ACME\")This will download ACME and all of its dependencies."
},

{
    "location": "gettingstarted.html#First-Steps-1",
    "page": "Getting Started",
    "title": "First Steps",
    "category": "section",
    "text": "We will demonstrate ACME by modeling a simple diode clipper. The first step is to load ACME:using ACMENow we create all the necessary circuit elements:j_in = voltagesource()\nr1 = resistor(1e3)\nc1 = capacitor(47e-9)\nd1 = diode(is=1e-15)\nd2 = diode(is=1.8e-15)\nj_out = voltageprobe()Specifying a voltagesource() sets up a voltage source as an input, i.e. the voltage it sources will be specified when running the model. Alternatively, one can instantiate a constant voltage source for say 9V with  voltagesource(9). The resistor and capacitor calls take the resistance in ohm and the capacitance in farad, respectively, as arguments. For the diode, one may specify the saturation current is as done here and/or the emission coefficient η. Finally, desired outputs are denoted by adding probes to the circuit; in this case a voltageprobe() will provide voltage as output.Next we need a Circuit instance to keep track of how the elements connect to each other:circ = Circuit()Connections can be specified by naming element pins that are connected:connect!(circ, j_in[\"+\"], r1[1])This connects the positive output of the input voltage source with pin 1 of the resistor. Alternatively, one can introduce named nets to which element pins connect. This may increase readability for nets with many connected elements, like supply voltages. Here, we use it for the ground net where we connect the negative side of the input voltage:connect!(circ, j_in[\"-\"], :gnd)One can also connect multiple pins at once:connect!(circ, r1[2], c1[1], d1[\"+\"], d2[\"-\"], j_out[\"+\"])\nconnect!(circ, :gnd, c1[2], d1[\"-\"], d2[\"+\"], j_out[\"-\"])Now that all connections have been set up, we need to turn the circuit description into a model. This could hardly be any easier:model = DiscreteModel(circ, 1./44100)The second argument specifies the sampling interval, the reciprocal of the sampling rate, here assumed to be the typical 44100 Hz.Now we can process some input data. It has to be provided as a matrix with one row per input (just one in the example) and one column per sample. So for a sinusoid at 1 kHz lasting one second, we do::y = run!(model, sin(2π*1000/44100*(0:44099).'))The output y now likewise is a matrix with one row for the one probe we have added to the circuit and one column per sample.More interesting circuits can be found in the examples located at Pkg.dir(\"ACME/examples\").In the likely event that you would like to process real audio data, take a look at the WAV package for reading writing WAV files.Note that the solver used to solve the non-linear equation when running the model saves solutions to use as starting points in the future. Model execution will therefore become faster after an initial learning phase.  Nevertheless, ACME is at present more geared towards computing all the model matrices than to actually running the model. More complex circuits may run intolerably slow or fail to run altogether."
},

{
    "location": "ug.html#",
    "page": "User Guide",
    "title": "User Guide",
    "category": "page",
    "text": ""
},

{
    "location": "ug.html#User-Guide-1",
    "page": "User Guide",
    "title": "User Guide",
    "category": "section",
    "text": ""
},

{
    "location": "ug.html#Element-Creation-1",
    "page": "User Guide",
    "title": "Element Creation",
    "category": "section",
    "text": "All circuit elements are created by calling corresponding functions; see the Element Reference for details."
},

{
    "location": "ug.html#Circuit-Description-1",
    "page": "User Guide",
    "title": "Circuit Description",
    "category": "section",
    "text": "Circuits are described using Circuit instances, created with Circuit(). Once a Circuit and elements have been created, the elements can be added to the circuit using the add! method:r = resistor(1e3);\nc = capacitor(22e-9);\ncirc = Circuit();\nadd!(circ, r)\nadd!(circ, c)Multiple elements can be added also be at once; the last two lines could have been replaced with add!(circ, r, c).In many cases, however, explicitly calling add! is not necessary. All that is needed is connect!, which connects two (or more) element pins. The elements to which these pins belong are automatically added to the circuit if needed. The only reason to explicitly call add! is to control the insertion order of sources or sinks, which determines the order in which inputs have to be provided and outputs are obtained.Pins are obtained from elements using []-style indexing, i.e. r[1] gives the first pin of the resistor defined above. So this connects the first pin of the resistor with the first pin of the capacitor:connect!(circ, r[1], c[1])Further connections involving the same pins are possible and will not replace existing ones. So this will effectively shorten the resistor, because now both of its pins are connected to c[1]:connect!(circ, r[2], c[1])Note that not all elements have numbered pins. For elements with polarity, they may be called + and -, while a bipolar transistor has pins base, collector, and emitter. The pins provided by each type of element are described in the Element Reference. Internally, the pin designators are Symbols. However, not all symbols are conveniently entered in Julia: :base is nice, symbol(\"1\") less so. Therefore, the [] operation on elements also accepts integers and strings and converts them to the respective Symbols. So r[symbol(\"1\")] is equivalent to r[1] and (assuming d to be a diode) d[:+] is equivalent to d[\"+\"] (but d[+] does not work).In addition to pins, connect! also accepts Symbols as input. This creates named nets which may improve readability for nets with many conneted pins:connect!(c[2], :gnd)\nconnect!(r[2], :gnd)Again, this only adds connections, keeping existing ones, so together with the above snippets, now all pins are connected to each other and to net a named gnd. It is even possible to connect multiple named nets to each other, though this will only rarely be useful."
},

{
    "location": "ug.html#Model-Creation-and-Use-1",
    "page": "User Guide",
    "title": "Model Creation and Use",
    "category": "section",
    "text": "A Circuit only stores elements and information about their connections. To simulate a circuit, a model has to be derived from it. This can be as simple as:model = DiscreteModel(circ, 1/44100)Here, 1/44100 denotes the sampling interval, i.e. the reciprocal of the sampling rate at which the model should run. Optionally, one can specify the solver to use for solving the model's non-linear equation:model = DiscreteModel(circ, 1/44100, HomotopySolver{SimpleSolver})See Solvers for more information about the available solvers.Once a model is created, it can be run:y = run!(model, u)The input u is matrix with one row for each of the circuit's inputs and one column for each time step to simulate. Likewise, the output y will be a matrix with one row for each of the circuit's outputs and one column for each simulated time step. The order of the rows will correspond to the order in which the respective input and output elements were added to the Circuit. To simulate a circuit without inputs, a matrix with zero rows may be passed:y = run!(model, zeros(0, 100))The internal state of the model (e.g. capacitor charges) is preserved accross calls to run!.Each invocation of run! in this way has to allocate some memory as temporary storage. To avoid this overhead when running the same model for many small input blocks, a ModelRunner instance can be created explicitly:runner = ModelRunner(model, false)\nrun!(runner, y, u)By using a pre-allocated output y as in the example, allocations in run! are reduced to a minimum.Upon creation of a DiscreteModel, its internal states (e.g. capacitor charges) are set to zero. It is also possible to set the states to a steady state (if one can be found) with:steadystate!(model)This is often desirable for circuits where bias voltages are only slowly obtained after turning them on."
},

{
    "location": "ug.html#ACME.SimpleSolver",
    "page": "User Guide",
    "title": "ACME.SimpleSolver",
    "category": "Constant",
    "text": "SimpleSolver\n\nThe SimpleSolver is the simplest available solver. It uses Newton iteration which features fast local convergence, but makes no guarantees about global convergence. The initial solution of the iteration is obtained by extrapolating the last solution found (or another solution provided externally) using the available Jacobians. Due to the missing global convergence, the SimpleSolver is rarely useful as such.\n\n\n\n"
},

{
    "location": "ug.html#ACME.HomotopySolver",
    "page": "User Guide",
    "title": "ACME.HomotopySolver",
    "category": "Constant",
    "text": "HomotopySolver{BaseSolver}\n\nThe HomotopySolver extends an existing solver (provided as the type parameter) by applying homotopy to (at least theoretically) ensure global convergence. It can be combined with the SimpleSolver as HomotopySolver{SimpleSolver} to obtain a useful Newton homtopy solver with generally good convergence properties.\n\n\n\n"
},

{
    "location": "ug.html#ACME.CachingSolver",
    "page": "User Guide",
    "title": "ACME.CachingSolver",
    "category": "Constant",
    "text": "CachingSolver{BaseSolver}\n\nThe CachingSolver extends an existing solver (provided as the type parameter) by storing found solutions in a k-d tree to use as initial solutions in the future. Whenever the underlying solver needs more than a preset number of iterations (defaults to five), the solution will be stored. Storing new solutions is a relatively expensive operation, so until the stored solutions suffice to ensure convergence in few iterations throughout, use of a CachingSolver may actually slow things down.\n\nSee M. Holters, U. Zölzer, \"A k-d Tree Based Solution Cache for the Non-linear Equation of Circuit Simulations\" for a more detailed discussion.\n\n\n\n"
},

{
    "location": "ug.html#Solvers-1",
    "page": "User Guide",
    "title": "Solvers",
    "category": "section",
    "text": "SimpleSolver\nHomotopySolver\nCachingSolverThe default solver used is a HomotopySolver{CachingSolver{SimpleSolver}}."
},

{
    "location": "elements.html#",
    "page": "Element Reference",
    "title": "Element Reference",
    "category": "page",
    "text": ""
},

{
    "location": "elements.html#Element-Reference-1",
    "page": "Element Reference",
    "title": "Element Reference",
    "category": "section",
    "text": ""
},

{
    "location": "elements.html#ACME.resistor",
    "page": "Element Reference",
    "title": "ACME.resistor",
    "category": "Function",
    "text": "resistor(r)\n\nCreates a resistor obeying Ohm’s law. The resistance r has to be given in Ohm.\n\nPins: 1, 2\n\n\n\n"
},

{
    "location": "elements.html#ACME.capacitor",
    "page": "Element Reference",
    "title": "ACME.capacitor",
    "category": "Function",
    "text": "capacitor(c)\n\nCreates a capacitor. The capacitance c has to be given in Farad.\n\nPins: 1, 2\n\n\n\n"
},

{
    "location": "elements.html#ACME.inductor-Tuple{Any}",
    "page": "Element Reference",
    "title": "ACME.inductor",
    "category": "Method",
    "text": "inductor(l)\n\nCreates an inductor. The inductance l has to be given in Henri.\n\nPins: 1, 2\n\n\n\n"
},

{
    "location": "elements.html#ACME.inductor-Tuple{Type{Val{:JA}}}",
    "page": "Element Reference",
    "title": "ACME.inductor",
    "category": "Method",
    "text": "inductor(Val{:JA}; D, A, n, a, α, c, k, Ms)\n\nCreates a non-linear inductor based on the Jiles-Atherton model of magnetization assuming a toroidal core thin compared to its diameter. The parameters are set using named arguments:\n\nparameter description\nD Torus diameter (in meters)\nA Torus cross-sectional area (in square-meters)\nn Winding's number of turns\na Shape parameter of the anhysteretic magnetization curve (in Ampere-per-meter)\nα Inter-domain coupling\nc Ratio of the initial normal to the initial anhysteretic differential susceptibility\nk amount of hysteresis (in Ampere-per-meter)\nMs saturation magnetization (in Ampere-per-meter)\n\nA detailed discussion of the paramters can be found in D. C. Jiles and D. L. Atherton, “Theory of ferromagnetic hysteresis,” J. Magn. Magn. Mater., vol. 61, no. 1–2, pp. 48–60, Sep. 1986 and J. H. B. Deane, “Modeling the dynamics of nonlinear inductor circuits,” IEEE Trans. Magn., vol. 30, no. 5, pp. 2795–2801, 1994, where the definition of c is taken from the latter. The ACME implementation is discussed in M. Holters, U. Zölzer, \"Circuit Simulation with Inductors and Transformers Based on the Jiles-Atherton Model of Magnetization\".\n\nPins: 1, 2\n\n\n\n"
},

{
    "location": "elements.html#ACME.transformer-Tuple{Any,Any}",
    "page": "Element Reference",
    "title": "ACME.transformer",
    "category": "Method",
    "text": "transformer(l1, l2; coupling_coefficient=1, mutual_coupling=coupling_coefficient*sqrt(l1*l2))\n\nCreates a transformer with two windings having inductances. The primary self-inductance l1 and the secondary self-inductance l2 have to be given in Henri. The coupling can either be specified using coupling_coefficient (0 is not coupled, 1 is closely coupled) or by mutual_coupling, the mutual inductance in Henri, where the latter takes precedence if both are given.\n\nPins: 1 and 2 for primary winding, 3 and 4 for secondary winding\n\n\n\n"
},

{
    "location": "elements.html#ACME.transformer-Tuple{Type{Val{:JA}}}",
    "page": "Element Reference",
    "title": "ACME.transformer",
    "category": "Method",
    "text": "transformer(Val{:JA}; D, A, ns, a, α, c, k, Ms)\n\nCreates a non-linear transformer based on the Jiles-Atherton model of magnetization assuming a toroidal core thin compared to its diameter. The parameters are set using named arguments:\n\nparameter description\nD Torus diameter (in meters)\nA Torus cross-sectional area (in square-meters)\nns Windings' number of turns as a vector with one entry per winding\na Shape parameter of the anhysteretic magnetization curve (in Ampere-per-meter)\nα Inter-domain coupling\nc Ratio of the initial normal to the initial anhysteretic differential susceptibility\nk amount of hysteresis (in Ampere-per-meter)\nMs saturation magnetization (in Ampere-per-meter)\n\nA detailed discussion of the parameters can be found in D. C. Jiles and D. L. Atherton, “Theory of ferromagnetic hysteresis,” J. Magn. Magn. Mater., vol. 61, no. 1–2, pp. 48–60, Sep. 1986 and J. H. B. Deane, “Modeling the dynamics of nonlinear inductor circuits,” IEEE Trans. Magn., vol. 30, no. 5, pp. 2795–2801, 1994, where the definition of c is taken from the latter. The ACME implementation is discussed in M. Holters, U. Zölzer, \"Circuit Simulation with Inductors and Transformers Based on the Jiles-Atherton Model of Magnetization\".\n\nPins: 1 and 2 for primary winding, 3 and 4 for secondary winding, and so on\n\n\n\n"
},

{
    "location": "elements.html#Passives-1",
    "page": "Element Reference",
    "title": "Passives",
    "category": "section",
    "text": "resistor\ncapacitor\ninductor(l)\ninductor(::Type{Val{:JA}})\ntransformer(l1, l2; coupling_coefficient=1,mutual_coupling=coupling_coefficient*sqrt(l1*l2))\ntransformer(::Type{Val{:JA}})"
},

{
    "location": "elements.html#ACME.voltagesource",
    "page": "Element Reference",
    "title": "ACME.voltagesource",
    "category": "Function",
    "text": "voltagesource()\nvoltagesource(v)\n\nCreates a voltage source. The source voltage v has to be given in Volt. If omitted, the source voltage will be an input of the circuit.\n\nPins: + and - with v being measured from + to -\n\n\n\n"
},

{
    "location": "elements.html#ACME.currentsource",
    "page": "Element Reference",
    "title": "ACME.currentsource",
    "category": "Function",
    "text": "currentsource()\ncurrentsource(i)\n\nCreates a current source. The source current i has to be given in Ampere. If omitted, the source current will be an input of the circuit.\n\nPins: + and - where i measures the current leaving source at the + pin\n\n\n\n"
},

{
    "location": "elements.html#Independent-Sources-1",
    "page": "Element Reference",
    "title": "Independent Sources",
    "category": "section",
    "text": "voltagesource\ncurrentsource"
},

{
    "location": "elements.html#ACME.voltageprobe",
    "page": "Element Reference",
    "title": "ACME.voltageprobe",
    "category": "Function",
    "text": "voltageprobe()\n\nCreates a voltage probe, provding the measured voltage as a circuit output.\n\nPins: + and - with the output voltage being measured from + to -\n\n\n\n"
},

{
    "location": "elements.html#ACME.currentprobe",
    "page": "Element Reference",
    "title": "ACME.currentprobe",
    "category": "Function",
    "text": "currentprobe()\n\nCreates a current probe, provding the measured current as a circuit output.\n\nPins: + and - with the output current being the current entering the probe at +\n\n\n\n"
},

{
    "location": "elements.html#Probes-1",
    "page": "Element Reference",
    "title": "Probes",
    "category": "section",
    "text": "voltageprobe\ncurrentprobe"
},

{
    "location": "elements.html#ACME.diode",
    "page": "Element Reference",
    "title": "ACME.diode",
    "category": "Function",
    "text": "diode(;is=1e-12, η = 1)\n\nCreates a diode obeying Shockley's law i=I_Scdot(e^v(eta v_T)-1) where v_T is fixed at 25 mV. The reverse saturation current is has to be given in Ampere, the emission coefficient η is unitless.\n\nPins: + (anode) and - (cathode)\n\n\n\n"
},

{
    "location": "elements.html#ACME.bjt",
    "page": "Element Reference",
    "title": "ACME.bjt",
    "category": "Function",
    "text": "bjt(typ; is=1e-12, η=1, isc=is, ise=is, ηc=η, ηe=η, βf=1000, βr=10)\n\nCreates a bipolar junction transistor obeying the Ebers-Moll equation\n\ni_E = I_SE cdot (e^v_E(eta_E v_T)-1)\n           - fracbeta_r1+beta_r I_SC cdot (e^v_C(eta_C v_T)-1)\n\ni_C = -fracbeta_f1+beta_f I_SE cdot (e^v_E(eta_E v_T)-1)\n           + I_SC cdot (e^v_C(eta_C v_T)-1)\n\nwhere v_T is fixed at 25 mV. The parameters are set using named arguments:\n\nparameter description\ntyp Either :npn or :pnp, depending on desired transistor type\nis Reverse saturation current in Ampere\nη Emission coefficient\nisc Collector reverse saturation current in Ampere (overriding is)\nise Emitter reverse saturation current in Ampere (overriding is)\nηc Collector emission coefficient (overriding η)\nηe Emitter emission coefficient (overriding η)\nβf Forward current gain\nβr Reverse current gain\n\nPins: base, emitter, collector\n\n\n\n"
},

{
    "location": "elements.html#Semiconductors-1",
    "page": "Element Reference",
    "title": "Semiconductors",
    "category": "section",
    "text": "diode\nbjt"
},

{
    "location": "elements.html#ACME.opamp-Tuple{}",
    "page": "Element Reference",
    "title": "ACME.opamp",
    "category": "Method",
    "text": "opamp()\n\nCreates an ideal operational amplifier. It enforces the voltage between the input pins to be zero without sourcing any current while sourcing arbitrary current on the output pins wihtout restricting their voltage.\n\nNote that the opamp has two output pins, one of which will typically be connected to a ground node and has to provide the current sourced on the other output pin.\n\nPins: in+ and in- for input, out+ and out- for output\n\n\n\n"
},

{
    "location": "elements.html#ACME.opamp-Tuple{Type{Val{:macak}},Any,Any,Any}",
    "page": "Element Reference",
    "title": "ACME.opamp",
    "category": "Method",
    "text": "opamp(Val{:macak}, gain, vomin, vomax)\n\nCreates a clipping operational amplifier where input and output voltage are related by\n\nv_textout = frac12cdot(v_textmax+v_textmin)\n                   +frac12cdot(v_textmax-v_textmin)cdot\n                    tanhleft(fracgfrac12cdot(v_textmax-v_textmin)cdot  v_textinright)\n\nThe input current is zero, the output current is arbitrary.\n\nNote that the opamp has two output pins, one of which will typically be connected to a ground node and has to provide the current sourced on the other output pin.\n\nPins: in+ and in- for input, out+ and out- for output\n\n\n\n"
},

{
    "location": "elements.html#Integrated-Circuits-1",
    "page": "Element Reference",
    "title": "Integrated Circuits",
    "category": "section",
    "text": "opamp()\nopamp(::Type{Val{:macak}}, gain, vomin, vomax)"
},

]}
