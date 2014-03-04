User Guide – the level-1 interface to C++QEDv2 {#userguide}
==============================================

* User guide
  * trajectory-state serialization
  * trajectory logging (especially EnsembleMCWF)
  * time averaging
  * version information
  * (heavier use of smart pointers)

\par Abstract
This Guide allows for understanding and using the framework on the highest level, that is, for building simulations out of the already existing modules and elementary systems as building blocks. These modules constitute the level-1 interface to C++QED, while the \ref structurebundleguide "level-2 interface" allows for the implementation of new elementary physical systems.

\tableofcontents


The basic idea is to allow users to build arbitrarily complex interacting quantum systems from elementary free subsystems and interactions (below these are commonly referred to as elements), and simulate their time evolution with a number of available time-evolution drivers. Operating with elementary physical systems, the interface is of a much higher level than that in the popular Quantum Optics Toolbox for Matlab \cite tan99, or the much less known, but venerable QSD (Quantum State Diffusion) C++ library \cite schack97. In both, the interface is built on quantum operators, and general usage involves a considerable amount of boilerplate.

C++QEDv2 specifies a small grammar to describe composite quantum systems. Apart from providing a number of elements out of the box, there are several tools which facilitate the implementation of new elements. These are being added continuously, as the need arises.

The framework saw a first release (v1), which was also partially documented in a journal article [EPJD]\_. As to physics, this paper still applies, and the reader should consult it if something is not clear from the physics point of view in the following. The differences between C++QEDv1 and v2 are so numerous that it is easier to describe the similarities, which essentially lie in the basic idea of the interface as described above. Because of this, and despite all the differences, experience has shown that it is not difficult migrate elements and |script| s written for v1 to v2. C++QEDv1 relied on a pure object-oriented design, while in v2 the emphasis is more on generic programming and template metaprogramming. In particular, the principal concept of the design of C++QEDv2 is that *all information available at compile time should be processed as such*, with the help of template metaprogramming. All in all, v1 can be regarded as a prototype of v2.

At the moment the following possibilities for time evolution are provided in the framework (cf. namespace #quantumtrajectory):

\link quantumtrajectory::Master Full Master equation\endlink
  Simple adaptive-stepsize evolution.

\link quantumtrajectory::MCWF_Trajectory Single Monte Carlo wave-function trajectory\endlink
  We use a modification of the original method with higher order adaptive-stepsize time evolution.

\link quantumtrajectory::EnsembleMCWF Ensemble of quantum (at present, Monte Carlo wave-function) trajectories\endlink
  These are evolved serially at the moment, parallelization should be implemented here.

A number of other methods, e.g. the quantum state diffusion, can be easily incorporated into the framework.

### Performance issues

The framework is very sensitive to performance both in terms of computer resources and coding/design. In the latter aspect the goal, as always in software architecture, is to create maximally reusable code. Perhaps the most spectacular example is that the very same code, if it is written generally enough, can be used to calculate Monte Carlo wave-function and full Master-equation evolution.

In the former aspect, there are physical and computational methods to increase performance. Among the physical ones, the most notable is the maximal use of interaction picture, which may help to get rid of very separate timescales in the problem. Among the computational ones we can mention

-   Maximal exploitation of special operator structures, i.e., sparse and tridiagonal matrices.

-   The use of adaptive-stepsize methods for evolving ordinary differential equations.

-   Judicious use of memory. The guideline is: If it is necessary to copy something, ask first whether it is *really* necessary.

In addition, we have to note that simulation of moving particles is inherently hard, since the Schrödinger equation is a partial differential equation, and we inevitably have to deal with both position and momentum representations, which are linked by Fourier transformation. In our problems, however, the particles are mostly moving in potentials created by electromagnetic fields, mainly standing and running waves. In this case we can stay in momentum space during the whole time evolution, and no |FFT| is necessary (cf. [EPJD]\_ and :ref:\`generalElements\_Particle\`).

A strange consequence is that in numerical physics the harmonic oscillator seems to be hard, while the cosine potential is easy.

Installation
------------


Installation {#userguideinstallation}
============

Writing and executing scripts {#userguidescripts}
=============================

An elementary example {#userguidescriptselementary}
---------------------

### Handling command-line parameters ### {#userguidescriptselementaryparameters}

Example: a binary system {#userguidescriptsbinary}
------------------------

Input/output of scripts {#userguidescriptsio}
-----------------------

Logging, also in EnsembleMCWF

Trajectory state i/o

More complex examples {#userguidescriptsmorecomplex}
---------------------

### Ring cavity ### {#userguidescriptsmorecomplexring}

### Multi-particle example ### {#userguidescriptsmorecomplexmultiparticle}
