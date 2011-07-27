***********************
Code map
***********************

=======================
General naming scheme
=======================


Each file in the framework revolve around a single concept, and is named after this concept. Cf. sec. :ref:`codeOrganization`.

``FilenameFwd.h``

  forward-declaring header

``Filename.h``

  normal header file

``impl/Filename.tcc``

  implementation header file containing template definitions, situated in subdirectory ``impl`` under the directory of the corresponding header

``Filename.cc``

  implementation header, containing non-template definitions

Headers under directories ``details`` are mostly very technical, often used for code generation.

======================
Directory ``utils``
======================

Here, such modules are defined as are more general than C++QED proper, so that the collection of these modules may in time become a library on their own. Headers are under directory ``include``, implementation files under ``src``.

``Algorithm``

  some additional STL-style algorithms, e.g. binary for_each

``ArrayTraitsFwd.h``

  forward-declaring the :class:`ArrayMemoryTraits` and :class:`ArrayTraversalTraits` classes necessary for :class:`Evolved`

``Blitz2FLENS``

  utilities for converting ``blitz::Array``\ s to FLENS vectors and matrices

``BlitzArrayExtensions``

  creating vector and matrix views of ``blitz::Array``\ s

``BlitzArray``

  template typedefs for real and complex arrays + a germinal Array class for wrapping ``blitz::Array``\ s

``BlitzArraySliceIterator``

  definition of :class:`blitzplusplus::basi::Iterator` together with its helpers

``BlitzArrayTraits``

  specializing the traits classes :class:`ArrayMemoryTraits` and :class:`ArrayTraversalTraits` for various ``blitz::Array``\ s

``BlitzTinyExtensions``

  helpers for ``blitz::TinyVector``\ s

``BlitzTiny``

  template typedefs for ``blitz::TinyVector``\ s used for characterising the size of multi-arrays and indexing them

``BlitzTinyOfArrays``

  class for wrapping ``blitz::TinyVector<blitz::Array<T,RANK>,LENGTH>``, solving the problem of such default-constructed classes

``BooleanNegatedProxy``

  class :class:`cpputils::BooleanNegatedProxy` which stores a reference to a boolean and when (implicitly) converted to a boolean, it returns the negated value of the original boolean

``CMatrix``

  defining the typedef :type:`linalg::CMatrix` and some helpers

``Combinatorics``

  utilities for combinatorics, at the moment class :class:`cpputils::CWR_Dir` which stores combinations with repetitions

``ComplexArrayExtensions``

  helpers for complex ``blitz::Array``\ s, e.g. Hermitian conjugation of multi-matrices

``ComplexExtensions``

  additional helpers for ``std::complex``

``Conversions``

  utilities related to Boost.NumericConversion (mostly unused at the moment)

``CVector``

  defining the typedef :type:`linalg::CVector`

``Evolved``

  defining :class:`~evolved::Evolved` and helpers

``EvolvedGSL`` & ``EvolvedNR``

  implementations of the abstract interface :class:`~evolved::Evolved` relying on GSL and the Numerical Recipes recipe, respectively

``Exception``

  base classes for exception classes

``FFT``

  generic FFT function in similar vein as the generic :class:`~evolved::Evolved`

``FormDouble``

  collecting to one place all the issues of the formatting of doubles, necessary for output in various places (cf. Stroustrup: The C++ Programming Language (special edition) 21.4.6.3.)

``Functional``

  some additional STL-style functionals (which cannot be expressed with Boost.Lambda)

``FuzzyDouble``

  a "fuzzy" double class, whose comparison accounts for an eventual error interval

``Hermite`` & ``HermiteCoefficients``

  calculating Hermite polynomials with the coefficients pre-defined up to a certain limit

``Integration``

  wrappers for numerical integration (not nicely done at the moment)

``MathExtensions``

  mathematical utilities partly implemented via GSL

``MultiIndexIterator``

  class :class:`cpputils::MultiIndexIterator`

``Operators``

  additional operator groups in the style of (and based on) Boost.Operators

``Pars``

  classes :class:`~parameters::Parameter` and :class:`~parameters::ParameterTable`

``Profiling``

  progress monitoring based on Boost.Timer

``Randomized``

  the general (?) frontend :class:`~randomized::Randomized` for random number generation

``Range``

  additional algorithms to Boost.Range (like binary for_each)

``SharedMemoryAllocator``

  **EXOTIC** an STL-style allocator for shared memory

``Simulated``

  class :class:`trajectory::Simulated`, an easy-to-use :class:`~trajectory::Trajectory`

``SimulatedHigh``

  collecting header for the convenience of those who want to "simulate"

``StochasticTrajectory`` & ``ParsStochasticTrajectory``

  classes :class:`trajectory::StochasticTrajectory` and :class:`trajectory::EnsembleTrajectories` together with related utilities

``TMP_Tools``

  template metaprogramming tools, extending (and based on) Boost.MPL

``Trajectory`` & ``ParsTrajectory``

  classes :class:`~trajectory::TrajectoryBase` and :class:`~trajectory::Trajectory` together with related utilities

``VectorFromMatrixSliceIterator``

  special case of :class:`blitzplusplus::basi::Iterator` for iterating over rows or columns of (multi-)matrices

``VectorTraits``

  **EXOTIC** specializing the traits classes :class:`ArrayMemoryTraits` and :class:`ArrayTraversalTraits` for ``std::vector``\ s with shared memory **BROKEN**

----------------------
Directory ``range_ex``
----------------------

This contains the future Boost.RangeEx library (algorithms for Boost.Range). To be removed once Boost.RangeEx gets incorporated into Boost.



==========================
Directory ``elements``
==========================

-------------------------
directory ``composites``
-------------------------

-------------------
directory ``frees``
-------------------

--------------------------
directory ``interactions``
--------------------------



==========================
Directory ``quantumdata``
==========================

``ArrayBase``, ``DensityOperator``, ``DimensionsBookkeeper``, ``StateVector``, ``Types``

  define classes :class:`quantumdata::ArrayBase`, :class:`quantumdata::DensityOperator`, :class:`DimensionsBookkeeper`, :class:`quantumdata::StateVector`, and :class:`quantumdata::Types`, respectively

``LazyDensityOperator``, ``LazyDensityOperatorSliceIterator``

  define classes :class:`quantumdata::LazyDensityOperator` and :class:`quantumdata::ldo::DiagonalIterator` together with the function :func:`quantumdata::partialTrace`

``NegPT``

  defines the function :func:`quantumdata::negPT`

``NonOrthogonalDensityOperator``, ``NonOrthogonalStateVector``, ``Transformation``

  defines the infrastructure for :ref:`representing quantumdata in non-orthogonal bases <quantumdataNonOrthogonal>` (covariant-contravariant formalism) **INCOMPLETE**


==============================
Directory ``quantumoperator``
==============================

``Sigma``, ``Tridiagonal``

  define classes :class:`quantumoperator::Sigma` and :class:`quantumoperator::Tridiagonal`, respectively


==============================
Directory ``scripts``
==============================


==============================
Directory ``structure``
==============================
