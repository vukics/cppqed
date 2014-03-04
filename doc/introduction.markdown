Introduction {#userguideintroduction}
============

\par Abstract
This guide gives a bird’s eye view of the background, history, and physical, computational, and design concepts of C++QED.

\tableofcontents

C++QED is a framework for simulating open quantum dynamics in general. It addresses the following problem: somebody who writes simulation code for a single quantum particle, mode, or spin today; tomorrow will want to simulate two or ten such systems interacting with each other. They want to reuse the code written for the single-system case, but this turns out to be rather difficult for naively written code, due to the algebraic structure of quantum mechanics. C++QED facilitates this task, providing a framework for representing elementary physical systems (particle motional degrees of freedom, harmonic oscillator modes, spins, etc.) in such a way that they can immediately be used as constituents of composite systems as well. Dynamical simulations of such systems can be performed with a set of tools provided by the framework.

Historically, C++QED originated from the simulation of systems in cavity quantum electrodynamics (CQED – hence the framework’s name). Its approach subsequently proved particularly useful in the wider context of quantum optics \cite vukics07b \cite vukics09 \cite niedenzu10, as this field typically addresses systems composed of several diverse “small” subsystems interacting with each other; but also in atomic physics or quantum many-body physics \cite vukics07b \cite maschler08 \cite nagy09.

The framework is capable of simulating fully quantum dynamics (Schrödinger equation) of any system, provided that its \link structure::Hamiltonian Hamiltonian\endlink is expressed in a finite discrete basis, and open (\link structure::Liouvillean Liouvillean\endlink) dynamics, if the Lindblad- or quantum-jump operators are expressed in the same basis. Apart from this, the only limitation on the simulated system is its dimensionality: with present-day hardware, a typical limit is
* a few millions dimensions for state-vector (e.g. \link quantumtrajectory::MCWF_Trajectory MCWF simulation\endlink), and
* a few thousands for density-operator manipulations (e.g. \link quantumtrajectory::Master Master equation evolution\endlink).

Since at present C++QED does not offer any special tools for addressing many-body problems, only “brute-force” calculations, this limitation means that such problems cannot be pursued very far.

C++QED saw a first, prototype version developed between 2006–2008, a purely object-oriented design, which was partially documented in \cite vukics07a. The second version, defined by the multi-array concept and compile-time algorithms, has been developed since 2008. By today it is quite mature, robust, and trustworthy, having passed thousands of tests also in real-life situations. It is an open-source project, hosted by [SourceForge.net](http://cppqed.sf.net), and it builds only on open-source libraries, such as
* [the Boost library collection](http://boost.org)
* [the GNU Scientific Library](http://www.gnu.org/software/gsl/manual) (GSL)
* [Blitz++](http://blitz.sf.net) and
* [the Flexible Library for Efficient Numerical Solutions](http://flens.sf.net) (FLENS).


Basic specification {#userguidebasicspecification}
===================

The framework provides building blocks (refered to as \ref mainpageelements "elements") for \link Composite composite quantum systems\endlink, which come in two types:
1. \link structure::Free elementary free subsystems\endlink, where “free” means that the system has only one quantum number;
2. \link structure::Interaction interactions\endlink between such free systems.

The representation of composites was designed in such a way that recursivity is possible in their definition: a composite system can act as a building block for an even more composite system. (As yet, no working example of this possibility has been implemented; but e.g. we may think about systems composed of several atoms with both external and internal degrees of freedom.)

Time-evolution drivers are provided for systems such composed. At the moment three kinds of time evolution are implemented:
1. \link quantumtrajectory::Master Full Master equation\endlink \cite carmichael solution for the density operator of the system. This is a simple \link evolved::Evolved adaptive-stepsize ordinary differential equation\endlink (ODE) evolution, whose dimension is the square of the dimension of the system Hilbert space.
2. \link quantumtrajectory::MCWF_Trajectory Single Monte Carlo wave-function (MCWF) trajectory\endlink for a stochastic state vector of the system \cite diosi86 \cite carmichael87 \cite dum92 \cite dalibard92 \cite gisin92 \cite molmer93 \cite soklakov02 \cite plenio98. \note This reproduces the normal Schrödinger-equation evolution in the case when the system is closed. \note We use a modified version of the original MCWF method in which the coherent (no-jump) evolution is performed with adaptive timestep, while the original method reduces to the completely inappropriate first-order Euler method for the coherent evolution \cite vukics07a. In the adaptive-timestep case, the timestep is limited not only by the \link evolved::Evolved ODE solver\endlink, but also by that the \link quantumtrajectory::mcwf::Pars::dpLimit total probability of quantum jumps\endlink must not be too high, in order to avoid two jumps in one timestep. Cf. \ref mcwftrajectory.
3. \link quantumtrajectory::EnsembleMCWF Ensemble of quantum (at present, MCWF) trajectories\endlink. Here, the trajectories are all evolved to certain time instants, at which points the approximate density operator of the system is assembled from the ensemble of stochastic state vectors. This density operator is then used to gain information about the state of the system, in exactly the same way as if it was obtained from a full Master-equation solution.

These modules together provide a high-level C++ application-programming interface. Normally, they are assembled in high-level C++ programs (which we will refer to as <em>script</em>s throughout), which specify the system and what to do with it. Example scripts are given in Sec. \ref userguidescripts. Hence, in normal usage to each physical system there corresponds a program which uses the framework as a library.

The framework strives to facilitate the implementation of new building blocks. New free systems and interactions are best implemented by inheritance from a class hierarchy (the classes making up the #structure namespace in the framework) which provides a lot of services. There are classes representing quantum operators with their appropriate algebra (the #quantumoperator namespace).

\see \ref structurebundleguide

A very important principle during the design and implementation of the framework was to process all information which is available at compile time, at compile time. This leads to optimal separation between the two types of information, the information available at compile time and the information becoming available only at runtime, and allows for maximal exploitation of compile time. Normally, a \ref userguidescripts "script" compiled once will be used many times for data collection.

We strive to recycle code as much as possible, the two most powerful ideas supporting this are
1. the reusal of the same \link #structure quantum-system interfaces\endlink for \ref masterequation "Master equation evolution" and \link quantumtrajectory::MCWF_Trajectory MCWF evolution\endlink, and
2. the common interface for \link quantumdata::StateVector state vectors\endlink and \link quantumdata::DensityOperator density operators\endlink for calculating quantum averages, cf. quantumdata::LazyDensityOperator.

The idea of composing systems out of elementary building blocks would fail if too many building blocks were required for realistic systems. Experience shows, however, that for a given problem domain, only a few such blocks are required for building arbitrarily complex systems in the domain. Furthermore, every such block has a clear physical meaning. The example of polarizable particles moving in optical (resonator or free) fields was presented in Sec. 5 of \cite vukics07a.


Compile-time algorithms {#userguidecompiletimealgorithms}
=======================
In compiled computer languages, source code goes through two stages to produce computational data:
1. the compilation resulting in an executable code and
2. the actual execution producing the data.

(In contrast, in interpreted languages these two stages are fused.) Usually, calculations are performed during stage 2 only. C++ template metaprogramming (TMP) provides a Turing-complete language \cite abrahams04 for stage 1, whose basic operands are the C++ *types*. Hence, stage 1 is useful not only for performing certain optimizations, but through TMP it becomes possible to shift calculations from stage 2 to here. This possibility was demonstrated in Appendix A in \cite Vukics2012CQEDv2.

Let us see how TMP comes about in the definition of composite quantum systems. Such systems have several quantum numbers, and their state vector is most conveniently represented by an entity having as many indices (this is referred to as the `rank` or `arity` of the system or state vector: un<em>ary,</em> bin<em>ary,</em> tern<em>ary,</em> quatern<em>ary,</em> etc.). Very many constructs in the framework have `RANK` as a template parameter (cf. the #quantumdata, #structure, #quantumtrajectory, and #quantumoperator namespaces).

Consider the following two possibilities:
1. The rank of the system is an information that becomes available only at runtime. This is how this information was treated in C++QEDv1.
2. Or, it is an information available already at compile time, the way it is treated in C++QEDv2.

If we have two state vectors of different arity \f$\Psi^{(\text{rank1})}\f$ and \f$\Psi^{(\text{rank2})}\f$, then in the first case they have to be represented by entities of the *same type,* while in the second case they *can be different types*. Therefore, a nonsensical expression like
\f[\Psi^{(\text{rank1})}+\Psi^{(\text{rank2})}\f]
causes a probably fatal error in the first case, which can be detected only at runtime, possibly after a lot of calculations. In the second case, however, such an expression can be a *grammar error* prohibiting compilation. Furthermore, for indexing such a state vector, in the first case a runtime loop is needed, which is not necessary in the second case, where this loop can be unravelled at compile time.

In C++QEDv2, since a script corresponds to a given physical system, the layout of the system is known at the time we compile the script, so that its arity is naturally also known. This then implies a lot of further compile-time calculations. Furthermore, many errors related to inconsistencies in the system layout can also be detected at compile time.

\see \ref metaprogrammingexample for a typical usecase of template metaprogramming in the framework.

Very roughly, we can think about C++QEDv2 scripts as C++ programs which exploit the template mechanism of C++, to generate (lower-level) C++ programs in such a way that in the resulting executable all compile-time information is encoded to yield a maximally effective and safe application for the given physical system.

In C++QEDv2, compile-time resource requirement scales with the complexity of the simulated system, while runtime resource requirement scales with the total dimensionality. Hence, compile time can be best exploited when the system is composed of a lot of subsystems, all with low dimensionality. As an example, we might think about several spin one-halves with complicated interactions.

\see Appendix A in \cite Vukics2012CQEDv2

Comparison with other software in the field {#userguideintroductioncomparison}
===========================================

In terms of simulation software available for open quantum systems, a very remarkable project is [QuTip](http://qutip.org/) \cite Johansson20121760 \cite Johansson20131234 , which appears to be a Python reimplementation of the originally Matlab-based Quantum Optics Toolbox \cite tan99. There, the interface is built on quantum operators, and general usage involves a considerable amount of boilerplate.

Since C++QED operates with elementary physical systems, the interface is of a much higher level than in the QuTip.

C++QEDv2 structure {#userguideintroductionstructure}
==================

\image html ProjectLayout.png "The different tracts of the framework and their relationships. (The set of displayed classes is rudimentary, with no pretension to completeness.)"

This structure is reflected on the build system, dependencies being very strictly defined and observed throughout.

The “very core” is composed of the following namespaces:
* #quantumdata: In this namespace the fundamental data structures are defined. It is completely autonomous, not depending on the rest of the framework, but all the rest directly or indirectly depends on it.
* #structure: In this namespace such interfaces are defined as quantum systems must or may present towards time-evolution drivers (or other clients). E.g. every system *must* be derived from the abstract interface class structure::QuantumSystem to be usable with the \link #quantumtrajectory Time-evolution drivers\endlink. A system *may* derive e.g. from structure::Hamiltonian, if its dynamics has a Hamiltonian part.
* #quantumtrajectory: Here, the time-evolution drivers are defined.

Within \ref mainpagecore "“core”", a separate tract is constituted by the #quantumoperator namespace, which defines classes for representing special operator structures to facilitate the implementation of new elements. At present, \link quantumoperator::Tridiagonal tridiagonal\endlink and \link quantumoperator::Sigma sparse\endlink matrices are implemented: this covers most of the problems we have come across in quantum optics so far.

Composites belong to “core” as they represent the most fundamental design concept of the framework. They are at the moment BinarySystem and Composite.

Elements come in two brands: in this tract, free elements are independent, while interactions depend on frees. Most free and interaction elements are implemented with the help of special operator structures, so that their implementation depends on the #quantumoperator namespace. Both brands of elements derive from at least one of the classes in namespace #structure.

Scripts use both the core and the elements of the framework for their implementation. Separate shared libraries are created out of the foregoing tracts, which may also link against the third-party libraries like GSL, FLENS, etc. Scripts in turn link against these C++QED libraries, so that the whole framework is stored in memory in a fully shared way.
