***********************
Code map
***********************

=======================
Naming scheme
=======================


Each file in the framework revolve around a single concept, and is named after this concept. Cf. sec. :ref:`codeOrganization`.

``FilenameFwd.h``

  forward-declaring header

``Filename.h``

  normal header file

``impl/Filename.tcc``

  implementation header file containing template definitions, situated in subdirectory ``impl`` under the directory of the corresponding header

``Filename.cc``

  implementation file, containing non-template definitions

Headers under directories ``details`` are mostly very technical, often used for code generation.


==========================
``elements``
==========================

-------------------------
``composites``
-------------------------

``Act``

  the :class:`Act` helper class for defining :class:`Composite`\ s

``BinarySystem``

  class :class:`BinarySystem`

``Composite``

  the :class:`Composite` class and the maker function :func:`makeComposite`

-------------------
``frees``
-------------------

This directory contains the definitions of elementary free building blocks. 

.. note::

  Since most of the files are extremely special here, and the directory's content is continually extending and is subject to changes according to the framework's usage, we are listing only the most typical files here.

``Mode``

  the classes :class:`(Pumped)(Lossy)Mode <Mode>`, together with numerous helpers to handle harmonic-oscillator mode operators, states, and dynamics

``MultiLevel``

  the class :class:`MultiLevel` 

``TimeIndependentMatrixHamiltonian``

  defines class :class:`TimeIndependentMatrixHamiltonian`

--------------------------
``interactions``
--------------------------

This directory contains the definitions of elementary interaction building blocks. 


==========================
``quantumdata``
==========================

``ArrayBase``, ``DensityOperator``, ``DimensionsBookkeeper``, ``StateVector``, ``Types``

  define classes :class:`quantumdata::ArrayBase`, :class:`quantumdata::DensityOperator`, :class:`DimensionsBookkeeper`, :class:`quantumdata::StateVector`, and :class:`quantumdata::Types`, respectively

``LazyDensityOperator`` & ``LazyDensityOperatorSliceIterator``

  define classes :class:`quantumdata::LazyDensityOperator` and :class:`quantumdata::ldo::DiagonalIterator` together with the function :func:`quantumdata::partialTrace`

``NegPT``

  the function :func:`quantumdata::negPT`

``NonOrthogonalDensityOperator``, ``NonOrthogonalStateVector``, ``Transformation``

  define the infrastructure for :ref:`representing quantumdata in non-orthogonal bases <quantumdataNonOrthogonal>` (covariant-contravariant formalism) **INCOMPLETE**


==============================
``quantumoperator``
==============================

``Sigma``, ``Tridiagonal``

  define classes :class:`quantumoperator::Sigma` and :class:`quantumoperator::Tridiagonal`, respectively


==============================
``quantumtrajectory``
==============================

``DO_Display``, ``EnsembleMCWF``, ``Master``, ``MCWF_Trajectory``, ``ProjectingMCWF_Trajectory``, & ``Pars...``

  define classes :class:`quantumtrajectory::DO_Display`, :class:`quantumtrajectory::EnsembleMCWF`, :class:`quantumtrajectory::Master`, :class:`quantumtrajectory::MCWF_Trajectory`, :class:`quantumtrajectory::ProjectingMCWF_Trajectory`, respectively & the corresponding parameter classes


``Evolution``

  the function :func:`evolve`

``EvolutionHigh``

  collecting header, to be included by scripts


==============================
``scripts``
==============================

This directory contains the clients of the framework, the actual programs corresponding to various physical systems. The same note applies as for the ``frees`` above, we are describing only some typical files.

``1particle1mode``

  a very flexible script describing the motion of a single polarizable particle in a cavity mode with different geometries

``QbitMode_C++QED``

  a qbit interacting with a harmonic-oscillator mode through Jaynes-Cummings interaction

``Ring``

  a single polarizable particle moving along two modes (ring cavity)


==============================
``structure``
==============================

``Averaged``, ``DynamicsBase``, ``ElementLiouvillean``, ``Exact``, ``FreeExact``, ``Free``, ``Hamiltonian``, ``Interaction``, ``Liouvillean``, ``QuantumSystem``, ``TridiagonalHamiltonian``

  classes :class:`structure::Averaged`, :class:`structure::DynamicsBase`, :class:`structure::ElementLiouvillean`, :class:`structure::Exact`, :class:`structure::FreeExact`, :class:`structure::Free`, :class:`structure::Hamiltonian`, :class:`structure::Interaction`, :class:`structure::Liouvillean`, :class:`structure::QuantumSystem`, :class:`structure::TridiagonalHamiltonian`, respectively

``ElementAveraged``

  defines class :class:`structure::ElementAveraged` and a generic implementation :class:`structure::averaged::DiagonalDO`

``MatrixOfHamiltonian``

  defines function :func:`calculateMatrix`

``Structure``

  collecting header, to be included by trajectories and composites

``SubSystem``

  helper classes for composites


======================
``utils``
======================

.. toctree::
  :maxdepth: 1

  Utility modules <codeMapUtils>

