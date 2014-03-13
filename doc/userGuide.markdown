User Guide – the level-1 interface of C++QED {#userguide}
============================================

\par Abstract
This Guide allows for understanding and using the framework on the highest level, that is, for building simulations (writing \ref glossaryscript "“scripts”") out of the already existing [modules](\coreMainPage) and \ref glossaryelement "elementary systems" as building blocks. These modules constitute the level-1 interface to C++QED, while the \ref structurebundleguide "level-2 interface" allows for the implementation of new elementary physical systems.

\tableofcontents

In the following we will cover how to use the framework on the highest level. The highest level is a C++ program, which is called a \ref glossaryscript. Such scripts can be written either in C++ or in Python, using the Python frontend of the framework.

\see [CpypyQED documentation](\cpypyqedMainPage)

A script creates an executable which defines and simulates a system of a particular layout. All information pertaining to the layout of the system is processed at compile time. Our compile-time algorithms can be regarded as C++ programs generating C++ programs in such a way that in the resulting executable all the compile-time information is encoded to yield a maximally efficient executable for the given layout.

A \ref glossaryscript is composed of a part where the system is specified, and another, where we do something with the system, most notably simulate its time evolution in one of the ways described in the \ref userguidebasicspecification.

In this section, the examples are taken from the field of (moving-particle) cavity QED.

\see \ref genericelements and \cite vukics07a for more information about the appearing physical elements

An elementary example {#userguideelementary}
=====================

The simplest case is that of a \ref glossaryelementfree system alone. Assume that this free system is a lossy resonator mode, which may be driven (“pumped” with a laser) as well. Its evolution can be described by the following Master equation:

\f[\dot\rho=\frac1{i\hbar}\lp\HnH\,\rho-\rho\HnH^\dag\rp+\kappa\,a\rho a^\dag,\quad\HnH=-\hbar(\delta+i\kappa)\adagger a+i\hbar\lp\eta\adagger-\hermConj\rp.\f]

Proof of principle {#userguideelementaryproofofprinciple}
------------------

\note The code fragments in this subsection are meant as proofs of principle. Cf. \ref userguideelementaryparameters for how a \ref glossaryscript will look like in practice.

We begin with defining the system, which is trivial in this case because the required free element is already given:

    PumpedLossyMode mode(delta,kappa,eta,cutoff);

where `cutoff` is the cutoff of the mode’s Fock space.

\see Mode.h, PumpedLossyMode

The system defined, we can already use it for several things, e.g. to calculate the action of the system Hamiltonian on a given state vector, or export the matrix of the Hamiltonian for exact diagonalization. Suppose we want to do something more: run a single MCWF trajectory. The system is started from a pure initial state, specified as

    quantumdata::StateVector<1> psi(mode::coherent(alpha,cutoff));

that is, the mode is in a coherent state with amplitude `alpha`. `StateVector<1>` means that it is a state vector of a system featuring a single quantum number (unary state vector).

\see quantumdata::StateVector and mode::coherent

Next, we define our Monte Carlo wave-function trajectory driver:

    quantumtrajectory::MCWF_Trajectory<1> trajectory(psi,mode,...a lot more parameters...);

The first two parameters are clear, the only thing to note is that `psi` is taken as a reference here, so the initial `psi` will be actually evolved as we evolve the trajectory. A lot more parameters are needed, pertaining to the ordinary differential equation stepper, the random number generation, etc., for which default values are also provided (see \ref userguideelementaryparameters "below").

\see quantumtrajectory::MCWF_Trajectory

### Running the trajectory ### {#userguideelementaryproofofprinciplerun}

The trajectory can be run with

    run(trajectory,time,dt);

This will evolve `trajectory` for time `time`, and \link trajectory::Trajectory::display display\endlink information about the state of the system after every time interval `dt`. The set of information (usually a \link structure::Averaged set of quantum averages\endlink) displayed is defined by the system. There is another version, which can be invoked as

    run(trajectory,time,dc);

Here, `dc` is an integer standing for the number of (adaptive) timesteps between two displays.

\note We believe that this version is actually more suited to the physics of the problem, since the timesteps will be small when many things are happening, and this is when we want more output, too.

\see trajectory::run


A full-fledged script {#userguideelementaryparameters}
---------------------

In the above, the necessary parameters must be previously defined somewhere. Parameters can come from several sources, but the most useful alternative is to have sensible defaults for all parameters, and to be able to override each of them separately in the command line when actually executing the \ref glossaryscript with a given set of parameters. This allows for a fine-grained control over what we want to accept as default and what we want to override, and the command line never gets crowded by those parameters for which the default is fine.

This possibility is indeed supported by the framework. Consider the following program: \include tutorialMode.cc \see The \link parameters parameter-bundle\endlink

This is a full-fledged script, so if you copy this into directory `CustomScriptsExample` as `temp.cc`, it will compile. Don't forget to re-run `cmake` in order to pick up the new script and call `make temp` in the build directory. (In fact, the script is already there in `CPPQEDscripts` as `tutorialMode.cc`.)

Let us analyse the script line by line: \dontinclude tutorialMode.cc \until ParameterTable

parameters::ParameterTable is a utility module storing all the parameters of the problem, and enabling them to be manipulated from the command line. It can store and handle any type for which i/o streaming operations are defined.

\skip evolution
\until Lossy
We instantiate the actual parameters for the time-evolution driver(s) and the mode, respectively. All the modules in the framework provide corresponding `Pars...` classes. \see evolution::Pars, mode::Pars, mode::ParsPumped, mode::ParsPumpedLossy

\skip evolution
\until follow
We specify the default values, e.g. that the default \link evolution::Method evolution method\endlink should be \link quantumtrajectory::Master Master equation\endlink.

\note Ultimate defaults are anyway given by the framework at the points where the `Pars...` class is defined, but since at those points there is no knowledge about the details of the problem, these cannot always qualify as “sensible”.

\skipline update

The command line is parsed by the parameters::update() function.

###Runtime implementation-selection for elements### {#userguideelementaryruntimeimplementationselection}

\skipline mode

We instantiate our mode, but now instead of the concrete PumpedLossyMode class, we are using a \link mode::make maker (dispatcher) function\endlink, which selects the best mode class corresponding to the parameters. The possibilities are numerous: a mode can be pumped or unpumped, lossy or conservative, finite or zero temperature, and it can be considered in several quantum mechanical pictures. Each of these possibilities is represented by a class, and these classes all have a common base ModeBase.

\note As demonstrated in Mode_.h, to implement such sets of classes in these situations, we use a combination of class composition, template,
and \refBoost{preprocessor metaprogramming,preprocessor}.

The mode::Ptr is a smart pointer to this type, so that it can “store” either of these classes. For example, if \f$\kappa=0\f$ and \f$\eta=0\f$, then we will have a #Mode; if \f$\eta\f$ is nonzero, a PumpedMode; and if both are nonzero, a PumpedLossyMode.

\note The significance of this is that e.g. if the mode is not lossy, then it is not that the *probability* of a quantum jump will be calculated and found zero, but rather that the *possibility* of a quantum jump will not even be considered during time evolution, which speeds the evolution up.

\see The full set of classes: Mode, ModeSch, PumpedMode, PumpedModeSch, LossyMode, LossyModeUIP, LossyModeSch, PumpedLossyMode, PumpedLossyModeUIP, PumpedLossyModeSch, and \ref genericelementsfreesmode "this summary".

`Sch`, `UIP` and no suffix mean Schrödinger picture, unitary interaction picture, and “full” interaction picture, respectively. If the system is not lossy, then the latter two coincide. Schrödinger picture is provided mostly for testing purposes, the performance is usually not optimal in this case.

\see \ref masterequationlimitations "This description" and structure::Exact::applicableInMaster

What we are telling the maker function in the same line is that the picture should be unitary interaction picture. Alternatively, we could add this as a parameter as well, which can be achieved by replacing `update` with `QM_Picture& qmp=updateWithPicture(p,argc,argv);`

###Initial condition, evolution### {#userguideelementaryinitalconditionevolution}

\skipline init

mode::StateVector is just another name for a unary quantumdata::StateVector, and mode::init is just another dispatcher, this time for the initial condition of the mode.

\skipline evolve

The function #evolve is a dispatcher for different evolution methods and the different versions of trajectory::run. So with this, the evolution mode can be changed from the command line, e.g. depending on the dimension of the problem – while for low dimension we may use Master equation, for larger dimension we may need to switch to an ensemble of quantum trajectories. The possibilities are:

- `--evol master` \link quantumtrajectory::Master Full Master equation\endlink

- `--evol single` Single \link quantumtrajectory::MCWF_Trajectory Monte Carlo wave-function trajectory\endlink

- `--evol ensemble` \link quantumtrajectory::EnsembleMCWF Ensemble of Monte Carlo wave-function trajectories\endlink

\note The function #evolve is overloaded in such a way that a \link quantumdata::DensityOperator density operator\endlink initial condition can also be given. In this case the evolved trajectory will always be a quantumtrajectory::Master

###More on command-line parameters### {#userguideelementarygenericparameters}

If we specify the `--help` option in the command line, the program will display all the available parameters together with their types (the less readable the more complex the type is), a short description, and the default value.

An example command line then looks like:

    PumpedLossyModeScript --eps 1e-12 --dc 100 --deltaC -10 --cutoff 20 --eta "(2,-1)" ...

\note There are some parameters that are “stronger” than others. E.g. if `--dc` is nonzero, then always the dc-mode trajectory::run function will be selected by the #evolve function above, regardless of the value of `--Dt`. The latter will only be considered if `--dc` is zero because in this case the deltaT-mode trajectory::run function will be selected. There is a similar relationship between `--minitFock` and `--minit`: the former will always override the latter.

Example: a binary system {#userguidebinary}
========================

Imagine we would like to define a more complex system, where a two-level atom (\ref genericelementsfreesqbit "qbit") interacts with a single cavity mode with a Jaynes-Cummings type interaction. Both the qbit and the mode may be pumped, and they may also be lossy. The full non-Hermitian Hamiltonian of the system in the frame rotating with the pump frequency reads:
\f[\HnH=-\hbar\lp\delta_\text{Q}+i\gamma\rp\,\sigma^\dagger\sigma+i\lp\eta_\text{Q}\,\sigma^\dagger-\hermConj\rp-\hbar\lp\delta_\text{M}+i\kappa\rp a^\dagger a+i\hbar\lp\eta_\text{M}\, a^\dagger-\hermConj\rp+i\hbar\lp g^*\sigma a^\dagger-\hermConj\rp\f]

Proof of principle {#userguidebinaryproofofprinciple}
------------------

\note The code fragments in this subsection are meant as proofs of principle. Cf. \ref userguidebinaryfullfledged for how a \ref glossaryscript will look like in practice.

First, we need to define the \ref glossaryelementfree elements of the system:

    PumpedLossyQbit qbit(deltaA,gamma,etaA);
    PumpedLossyMode mode(deltaC,kappa,etaC,cutoff);

Next, we define the interaction between them:

    JaynesCummings act(qbit,mode,g);

Now we have to bind together the two free subsystems with the interaction, which is simply accomplished by:

    binary::Ptr system(binary::make(act));

In the case of a BinarySystem the complete layout of the system can be trivially figured out from the single interaction element. Similarly to mode::Ptr, binary::Ptr is an entity which can store all the versions of BinarySystem (for example, conservative or lossy).

\note BinarySystem is an extremely powerful module, whose design reflects the basic idea behind the framework. It internally handles all the \link blitzplusplus::basi loops and slicing\endlink that are necessary to calculate e.g. the effect of the Hamiltonian of the qbit component if it is part of a binary system. It acts and feels like any other system, like e.g., Qbit itself, the difference being that the latter has only one quantum number, while BinarySystem has two.

Our next task is to define the initial condition:

    StateVector<2> psi(qbit::state0()*mode::coherent(alpha,cutoff));

This is to say that the qbit is in its \f$\ket0\f$ state, and the mode is in a coherent state \f$\ket\alpha\f$. Both states are of type quantumdata::StateVector `<1>`, meaning that they are state vectors featuring a single quantum number, and `*` means direct product of state vectors, so the result here is a quantumdata::StateVector `<2>`. Direct product is *not commutative*, and we have to comply with the order defined above for the free systems, which is Qbit followed by Mode.

Alternatively, we could have said

    StateVector<2> psi(init(pq)*init(pm));

From this point on, usage is the same as we have seen above for the mode example. Since in this case the system is a BinarySystem, it will reach into its constituents for what system characteristics to display, supplying either with the corresponding \link quantumdata::LazyDensityOperator reduced density operator\endlink (cf. \ref slicinganldo), which contains all information on the state of the subsystem.

A full-fledged script for a binary system {#userguidebinaryfullfledged}
-----------------------------------------

If the system is not to be used for anything else, just for being \link evolve evolved\endlink, we can avoid having to invent all these redundant names like `qbit`, `mode`, `act`, `system`, `trajectory`, and create everything in place. In this case a full-fledged script can be as terse as (cf. `CPPQEDscripts/tutorialBinary.cc`):

\include tutorialBinary.cc

Here qbit::make will dispatch exactly \ref genericelementsfreesqbit "the same possibilities" that we have seen for the mode above.


More complex examples {#userguidemorecomplex}
=====================

If there are more than two \ref glossaryelementfree subsystems, the system can be much more complex. The number of possible interactions rises combinatorically with the number of frees. This is the situation when the full potential of C++QED is displayed.

Ring cavity {#userguidemorecomplexring}
-----------

Assume we want to define a system where a particle is moving along the axis of a ring cavity, and is interacting with two counterpropagating running-wave modes of the cavity. Both of the modes are lossy, and one of them is pumped; the particle is not pumped \cite niedenzu10. This system consists of three subsystems, a \ref genericelementsfreesparticle "particle", and the two \ref genericelementsfreesmode "modes". There are three interactions:

- (1-2) The particle can absorb from either of the modes and emit it in a stimulated way into the same mode. This yields dipole force for the particle and a corresponding light shift for the mode. It is implemented by the interaction element ParticleAlongCavity.

- (3) The particle can emit into the other mode. This yields a ternary interaction between all the frees, implemented by ParticleTwoModes.

We can lay out the system as the following simple network:

\image html NetworkRing.png "Layout of composite system in the Ring-cavity example"

The Composite module of the framework is designed to represent such a network. Assume the following definitions are in effect:

    // Instantiate Frees
    Particle        part (...);
    LossyMode       plus (...);
    PumpedLossyMode minus(...);
    // Instantiate Interactions
    ParticleAlongCavity actP(plus ,part,...,MFT_PLUS );
    ParticleAlongCavity actM(minus,part,...,MFT_MINUS);
    ParticleTwoModes    act3(plus,minus,part,...);

Here `MFT_` means the type of the mode function and can be `PLUS`, `MINUS`, `COS`, and `SIN` \cite vukics07a.

Then the system can be created by invoking the maker function for Composite with the special syntax:

    composite::make(_<1,0>  (actP),
                    _<2,0>  (actM),
                    _<1,2,0>(act3));

What we are expressing here e.g. with the specification `_<1,2,0>(act3)` is that the 0th “leg” of the interaction element ParticleTwoModes, which is the mode `plus`, is the 1st in our row of frees in the network above.

\note Following C/C++ convention, all ordinals begin with 0 in the framework.

The 1st leg, the mode `minus` is the 2nd in the row; and the 2nd leg, the particle is the 0th in the row of frees. The legs of interaction elements cannot be interchanged, and we also have to be consistent with our preconceived order of frees throughout. The three `_` objects above contain all the information needed by the framework to figure out the full layout of the system.

Any inconsistency in the layout will result in a compile-time or runtime error (cf. the description @ Composite). The user is encouraged to play around creating layout errors deliberately, and see what effect they have.

The actual C++ type of a Composite object returned by such an invocation of composite::make is quite complex, but a maker metafunction facilitates referencing it

    composite::result_of::Make<_<1,0>,_<2,0>,_<1,2,0> >::type // type of the above Composite

\note We can use the `auto` keyword to let the compiler figure out the type returned by composite::make, if we need to store the object in a variable.

A full-fledged script in the terse way may read as:

\dontinclude tutorialComposite.cc
\until pe);
\skip }
\skip ptm
\skip }
\until }

A notable additional feature as compared to previous examples is that since now we have two modes in the system, we somehow have to differentiate between their parameters in the command line. This is achieved by the `"P"` and `"M"` modifiers added to the constructors of `Pars...` objects, so that e.g. instead of `--cutoff` we now have the separate options `--cutoffP` and `--cutoffM`.

Although now all the frees have the general types contained by the mode::Ptr classes, their possible types are still restricted by the `Pars...` classes, such that e.g. `plus` can never become pumped.

\note Most of the interactions in C++QED will be binary, here we have seen an example for a ternary interaction. So far, we have encountered only a single example for a quaternary interaction.

Multi-particle example {#userguidemorecomplexmultiparticle}
----------------------

Finally, we are reviewing one more example, which displays another feature reflecting a basic principle of quantum physics: if two systems are identical, they are indistinguishable. In the language of C++QED this means that a single object is enough to represent them.

Consider two identical pumped \ref genericelementsfreesparticle "particles" moving in a direction orthogonal to the axis of a cavity sustaining a single lossy mode \cite vukics07b . The layout of the system is:

\image html NetworkSelforg.png "Layout of composite system in the two-particle example"

The kernel of the corresponding script reads:

\dontinclude tutorialComposite.cc
\skip LossyMode<>
\until pe);

The element IdenticalParticles has no dynamical role, it can be used to calculate and display characteristics of the quantum state of the two-particle subsystem.


Input/output of scripts {#userguideio}
=======================

Each trajectory output has a header part. First, the full command line is echoed, then C++QED and other linked library versions are displayed.

This is followed by a summary of the parameters of the given simulation. The header displays system characteristics, and lists the decay channels of the system, each decay channel endowed with an ordinal number. This list facilitates the interpretation of \ref userguideiologging "log output".

Following the header part, the time-dependent simulated data is displayed, organized into columns. The first two columns are time and timestep, then separated by tab characters, the data stemming from the different subsystems follows. A key to the data with a short description of each column is provided at the end of the header part.

\anchor userguideioprecision

The output of real numbers has a precision of three digits by default, this can be overridden by the `--precision` option. This characteristic propagates throughout the framework starting from trajectory headers:

- Trajectory header: real-number parameters of overall precision, but at least 6
- Time and timestep: of overall precision, but at least 6
- Quantum averages: of overall precision

Logging {#userguideiologging}
-------

Log output during or at the end of the trajectory is governed by the `--logLevel` option. When nonzero, the trajectory drivers will record certain data during the run.

### Single tarjectory ### {#userguideiologgingsingle}

For a \link quantumtrajectory::MCWF_Trajectory single MCWF trajectory\endlink, there are several logging levels:

- `logLevel<=1` No log output during the trajectory.
- `logLevel>0` At the trajectory’s end, a summary and a list of the jump time instants and jump kinds (ordinal number of the jump which occured) – this list is usually considered the actual quantum trajectory.
- `logLevel>1` Reporting jumps also during the trajectory evolution.
- `logLevel>2` Reporting dpLimit overshoots and the resulting stepsize decrease also during trajectory evolution.
- `logLevel>3` Reporting the number of failed ODE steps in the given step of the trajectory evolution.

\see The discussion at quantumtrajectory::MCWF_Trajectory and the page \ref mcwftrajectory. Furthermore, the paragraph "An important note concerning sampling" at quantumtrajectory::TimeAveragingMCWF_Trajectory.

\note If the timestep-control is ODE-dominated, we can get several failed ODE steps per timestep, on the other hand, if it’s dpLimit-dominated, we will get dpLimit overshoots in the majority of timesteps.

### Ensemble of tarjectories ### {#userguideiologgingensemble}

For an \link quantumtrajectory::EnsembleMCWF ensemble of MCWF trajectories\endlink, if `logLevel>0`, we get accumulated values of the logged characteristics of the element \link quantumtrajectory::MCWF_Trajectory MCWF trajectories\endlink.

This is followed by a time histogram of the occured quantum jumps in all decay channels. The preparation of the histogram is governed by two parameters: `--nBin`, the number of bins, and `--nJumpsPerBin`. If the former is zero, then on the basis of latter a heuristic determination of the number of bins occurs in such a way that the average number of jumps equal `--nJumpsPerBin`. If the total number of jumps is smaller than `2*nJumpsPerBin`, no histogram is created.

Time Averaging {#userguideiotimeaveraging}
--------------

For a single \link quantumtrajectory::MCWF_Trajectory MCWF trajectory\endlink, time averaging of the displayed system characteristics can be turned on with the `--timeAverage` switch. The time averages will be displayed at the end of the trajectory organized into columns in the same way as the displayed system characteristics are during the run.

A relaxation time which is left out from averaging at the beginning of the trajectory can be specified with the `--relaxationTime` option.

Output into file {#userguideiointofile}
----------------

Trajectory output goes to the standard output by default.

Alternatively, an output file can be specified with the `--o` option. In this case, \ref userguideiologging "log output" during the trajectory will still go to the standard output, which can be piped into another file.

Input/output of trajectory state {#userguideiotrajectorystate}
--------------------------------

The feature described in this subsection relies on \refBoost{Boost.Serialization,serialization} installed on the system. If this library is not available, these features get silently disabled.

For all the trajectory types (quantumtrajectory::Master, quantumtrajectory::MCWF_Trajectory, quantumtrajectory::EnsembleMCWF) it is possible to save the complete state of the trajectory. E.g., for a single \link quantumtrajectory::MCWF_Trajectory MCWF trajectory\endlink, this involves the random-number generator internal state, the ODE-stepper internal state, and the full state vector of the system (plus log state and eventually also the time-averaging data).

If the trajectory is started with a \ref userguideiointofile "specified output file" as `--o <filename>`, then at the end of the trajectory the final state will be saved into a file named `<filename>.state`. This allows for eventually resuming the trajectory to continue for a longer time, which results in exactly the same run as if the trajectory has not been halted.

Any state file can also be used as initial condition for a new run through the `--initFileName` option.

\note At trajectory-resumption or when starting a trajectory with an initial state resulting from a previous run, parameters that do not affect the dimensionality of the system can even be changed. One can simulate sudden parameter changes (e.g. the switching-off of a pump laser) in this way.

It is also possible to damp out trajectory states during the run. The frequency of state output is defined by the `--sdf` option. This allows for postanalysing the state of the system with machine precision.

\note Trajectory state i/o can be handled by the [Python interface](\cpypyqedMainPage) as well.

Assessing entanglement {#userguideentanglement}
======================

The feature described in this subsection relies on the [FLENS library](http://flens.sf.net), without which it gets silently disabled.

In a composite system we may want to assess the entanglement between two parts of the system. This can be done using the negativity of the density operator’s partial transpose \cite vidal02. Since this quantity depends on the density operator in a non-linear way, this makes sense only in the case of \link quantumtrajectory::Master Master-equation evolution\endlink or an \link quantumtrajectory::EnsembleMCWF ensemble of quantum trajectories\endlink.

If we want to assess the entanglement for a single trajectory, this can be achieved by choosing \link evolution::ENSEMBLE ensemble evolution \endlink and setting the number of trajectories to 1.

Since the negativity of the partial transpose can be used for a bipartite system only, we have to partition our composite system into two. The subsystem to be considered as one part of the two has to be specified in an additional compile-time integer list to the evolve() function.

Assume e.g. that in the \ref userguidemorecomplexmultiparticle "above example" we are looking for the entanglement between the two particles together as one part, and the mode as the other part. Then the invocation of #evolve is modified as

    evolve<1,2>(psi,system,pe);

We simply have to list in the compile-time argument list the ordinals of those \ref glossaryelementfree "frees" that consist one part. Of course, in this case this is equivalent to

    evolve<0>(psi,system,pe);

The negativity will appear as the last column in the system characteristics display during the run, separated by a tab character from the other columns.

\see quantumdata::negPT

Some common command-line options of scripts {#userguidescriptoptions}
===========================================

\par `--help`
display the list of options

\par `--version`
git or release version information about the framework and underlying libraries

Trajectory parameters
---------------------

\par `--T <double>`
simulated time

\par `--dc <int>, --Dt <double>`
output frequency in number of ODE steps or in time, respectively – when nonzero, the former is used, cf. \ref userguideelementaryproofofprinciplerun

\par `--o <filename>`
output file, cf. \ref userguideiointofile

\par `--precision <int>`
number of digits in the output of doubles throughout the framework, cf. \ref userguideioprecision "here"

\par `--displayInfo, --no_displayInfo`
switching header output on (default) and off, cf. \ref userguideio

\par `--sdf <int>`
state display frequency, cf. \ref userguideiotrajectorystate

\par `--eps <double>, --epsAbs <double>`
ODE stepper relative and absolute precision, cf. evolved::TimeStepBookkeeper

\see trajectory::Trajectory and trajectory::Adaptive

Stochastic trajectory parameters
--------------------------------

\par `--seed <long>`
random number generator seed

\par `--noise, --no_noise`
switching noise on (default) and off

\par `--nTraj <long>`
number of trajectories for ensemble averaging

\see trajectory::StochasticTrajectory

MCWF trajectory parameters
--------------------------

\par `--dpLimit <double> --overshootTolerance <double>`
total jump probability limit and jump probability overshoot tolerance factor, cf. the discussion at quantumtrajectory::MCWF_Trajectory and the page \ref mcwftrajectory

\par `--logLevel <int>`
governs the level of log output, cf. \ref userguideiologgingsingle

\par `--nBins <int> --nJumpsPerBin <int>`
governs the properties of jump histograms, cf. \ref userguideiologgingensemble

Highest-level evolution parameters
----------------------------------

\par `--evol <EvolutionMethod>`
evolution method, cf. \ref evolution::Method

\par `--negativity, --no_negativity`
switching entanglement calculation on and off (default), cf. \ref userguideentanglement

\par `--timeAverage, --no_timeAverage`
switching time averaging for a single MCWF trajectory on and off (default), cf. \ref userguideiotimeaveraging

\par `--relaxationTime <double>`
relaxation time for time averaging
