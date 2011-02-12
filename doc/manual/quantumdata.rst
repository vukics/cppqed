.. _quantumdata:

=====================================
The ``quantumdata`` namespace
=====================================

::

  namespace quantumdata

The ``quantumdata`` namespace comprises classes which represent the state of composite quantum systems and provide various interfaces to manipulate this data. Some of its most important classes fit into a single class hierarchy, which may be sketched as follows:

.. image:: figures/quantumdata3.png
   :width: 608

.. _quantumdataTemplates:

.. note::

  Template argument definitions:
  
  ``int RANK``, ``int RANK1``, ``int RANK2``
    Positive integer standing for the number of elementary Hilbert spaces (arity of the state vector)

  ``bool IS_CONST``
    Governs the constness of the corresponding class.

  ``typename B``
    An optional base class for base-class chaining. An empty class (``quantumdata::details::Empty``) by default.

  ``typename F``
    A functor.
    
  ``typename V``
    A compile-time vector with the same functionality and characteristics as :ref:`above <basiTemplates>`.

  ``typename TRAFO``
    Type representing a metrical transformation

------------------------
Utilities
------------------------

.. toctree::

  quantumdataUtils

----------------------------------------------------
Common interface for calculating quantum averages
----------------------------------------------------

In a quantum-simulation framework, users should be given the possibility to write code able to calculate quantum expectation values from quantumdata, independently of whether this data is represented by state vectors or density operators, in orthogonal or non-orthogonal bases. One obvious solution is relying on the formula

.. math::

  \avr{A}=\Tr{A\rho}

(:math:`A` being an observable and :math:`\rho` the density operator of the system), to write code only for the density-operator case, and fall back to this in the state-vector case as well, by calculating a dyad from the state vector. This is, however, extremely wasteful, since usually not all the matrix elements of :math:`\rho` are needed for calculating the average, furthermore, for large dimensionality this solution may become outright unaffordable in terms of memory.

The solution adopted to this problem in C++QED is represented by the class :class:`~quantumdata::LazyDensityOperator`, which provides a common interface for all the four cases :class:`~quantumdata::StateVector`, :class:`~quantumdata::DensityOperator`, and their non-orthogonal versions, to calculate quantum averages from their data. The "laziness" is constituted by that in the case of state vectors, only those elements of the density operator are calculated that are asked for.

.. py:module:: LazyDensityOperator.h
   :synopsis: Defines LazyDensityOperator abstract interface

.. class:: quantumdata::LazyDensityOperator

  :ref:`template parameters <quantumdataTemplates>`: RANK; inherits publicly from :class:`DimensionsBookkeeper`

  This class is totally inmutable, all its member functions being constant.

  .. type:: Dimensions

    Inherited from :class:`DimensionsBookkeeper`

  .. type:: Idx

    The type used for indexing the "rows" and the "columns", this is either a normal, or a multi-index, in the latter case represented by a :c:macro:`TTD_IDXTINY`::

      typedef typename mpl::if_c<(RANK==1),int,TTD_IDXTINY(RANK)>::type Idx;
      // Idx is just an int if RANK==1, otherwise a tiny of ints

  .. function:: LazyDensityOperator(const Dimensions& dims)

    Constructor (protected)

  .. function:: const dcomp operator()(const Idx& i, const Idx& j) const

    The indexing function, purely virtual.

  .. function:: double operator()(const Idx& i) const

    An inline function for conveniently addressing the diagonal elements::

      double operator()(const Idx& i) const {return real((*this)(i,i));}

  .. function:: begin(V v) const

  .. function:: end(V v) const

    :ref:`template parameters <quantumdataTemplates>`: V

    These functions return the :class:`quantumdata::ldo::DiagonalIterator`\ s corresponding to the start and end of the sequence of slices defined by ``V`` (cf. the :ref:`section on slicing a LazyDensityOperator <ldoSlicing>`). 


**Semantics:**

  *Unary system*
    Imagine having a mode represented in Fock basis with ladder-operator :math:`a`. To calculate the quantum expectation value 

    .. math::

      \avr{a^2}=\Tr{a^\dag b\rho}=\sum_i\sqrt{i(i-1)}\,\rho_{i;i-2},

    one can write the following function::

      include "LazyDensityOperator.h"

      const dcomp calculateASqr(const LazyDensityOperator<1>& matrix)
      {
        dcomp res;
        for (int i=2; i<matrix.getTotalDimensions(); ++i) res+=sqrt(i*(i-2))*matrix(i,i-2);
        return res;
      }

.. _calculateADaggerB:

  *Binary system*
    Imagine having two modes represented in Fock bases with ladder-operators :math:`a` and :math:`b`, respectively. To calculate the quantum expectation value 

    .. math::

      \avr{a^\dag b}=\Tr{a^\dag b\rho}=\sum_{i,j}\sqrt{(i+1)j}\,\rho_{(i,j);(i+1,j-1)},

    one can write the following function::

      include "LazyDensityOperator.h"

      const dcomp calculateADaggerB(const LazyDensityOperator<2>& matrix)
      {
        typedef LazyDensityOperator<2>::Idx Idx;
        const LazyDensityOperator<2>::Dimensions dim(matrix.getDimensions());

        dcomp res;
        for (int i=0; i<dim[0]-1; ++i) for (int j=1; m<dim[1]; ++j) res+=sqrt((i+1)*j)*matrix(Idx(i,j),Idx(i+1,j-1));
        return res;
      }


---------------------------
High-level data structures
---------------------------

The :class:`~quantumdata::StateVector` and :class:`~quantumdata::DensityOperator` classes (and their non-orthogonal counterparts) exist for two main reasons:
  * As interfaces to :type:`~quantumdata::Types::StateVectorLow` and :type:`~quantumdata::Types::DensityOperatorLow`, respectively, which are more convenient to use on higher levels of the framework, e.g. in scripts, and in quantum trajectories. This is especially because while ``blitz::Array`` uses by-reference copy semantics and expression templates, these classes use the more usual by-value copy semantics and normal semantics for arithmetic operations. This means, however, that copying and arithmetics should be used judiciously, if possible only in the startup phase of simulations.

  * As an implementation of :class:`~quantumdata::LazyDensityOperator`


State vector
^^^^^^^^^^^^^^^

.. toctree::

  quantumdataStateVector


Density operator
^^^^^^^^^^^^^^^^^

.. toctree::

  quantumdataDensityOperator


Representing state vectors & density operators in non-orthogonal bases
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We rely on the :ref:`covariant-contravariant formalism <NonOrthogonalFormalism>` to represent these in the framework.

The user has considerable freedom in how he defines the metrical transformation for a given non-orthogonal basis, the only requirement being that there is a specialized instant of :class:`quantumdata::transformation::Traits` for the type representing the transformation. When compositing elementary Hilbert spaces, the framework takes care about the correct composition of the elementary metrical transformations. Of course, once there is a single non-orthogonal component, the whole composite has to be represented by non-orthogonal classes, in which case the system substitutes a :class:`quantumdata::transformation::Identity` dummy transformation for the metrics of such components as are represented in orthogonal bases.


The infrastructure for defining metrical transformations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. toctree::
   :maxdepth: 2

   The infrastructure for defining metrical transformations <quantumdataTransformations>


Non-orthogonal state vector & density operator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.. toctree::
   :maxdepth: 2

   quantumdataNonOrthogonalStateVector

   quantumdataNonOrthogonalDensityOperator






.. _ldoSlicing:

------------------------------------------------------
Slicing of ``LazyDensityOperator``
------------------------------------------------------

:ref:`Analogously to state vectors <basiSlicing>`, it is also necessary to slice :class:`~quantumdata::LazyDensityOperator`\ s because for calculating quantum expectation values of subsystem-observables (e.g. in :class:`Composite`\ s), the partial-trace density operator is needed. For the partial trace, however, only such elements of the full density operator are needed as are diagonal in the indeces *not* belonging to the given subsystem (dummy indeces). This is the reason why the tool performing the iteration is called :class:`~quantumdata::ldo::DiagonalIterator`.

Notes on implementation
^^^^^^^^^^^^^^^^^^^^^^^^^

.. toctree::
   :maxdepth: 2

   Notes on implementation <quantumdataLazyDensityOperatorSlicing>

Slicing is fully recursive, a sliced :class:`~quantumdata::LazyDensityOperator` (usually obtained by dereferencing a :class:`~quantumdata::ldo::DiagonalIterator`) can be further sliced.


Partial trace
^^^^^^^^^^^^^^

On higher levels of the framework, the iteration is performed via the function:

.. function:: const T quantumdata::partialTrace(const quantumdata::LazyDensityOperator<RANK>& matrix, F f, V v, T t)

  :ref:`template parameters <quantumdataTemplates>`: RANK, F, V, T


  ``T`` is an arithmetic type which must be default-constructible.
    The trailing function argument is only a dummy one, facilitating template-argument deduction.

    .. todo::

      (Re)consider whether it is not better to use the last argument's value as initial value in the accumulation. In this case one has to be careful e.g. that it is not correct to assign to a default-constructed ``blitz::Array``.

  ``F`` is a functor with signature
    ``const T(const typename ldo::DiagonalIterator<RANK,V>::value_type&)``

    ::

      template<int RANK, typename F, typename V, typename T>
      const T
      partialTrace(const LazyDensityOperator<RANK>& matrix, F function, V v, T)
      {
        ldo::DiagonalIterator<RANK,V> begin(matrix.begin(v));
	T init(function(*begin));

	using namespace cpputils;
	return T(
	         accumulate(++begin,
		            matrix.end(v),
		            init,
		            function,
		            cpputils::plus<T>()
		           )
	        );
      }


  **Semantics:**
    The function iterates through all the combinations of the dummy indeces, for each slice it takes the value returned by the functor ``function``, and accumulates these values. In the following we give some instructive examples of usage.

    Calculating the full partial density operator of a unary subsystem:
      Assume having a function which simply copies the content of a unary :class:`~quantumdata::LazyDensityOperator` into a :type:`structure::free::DensityOperatorLow`::

        using structure::free::DensityOperatorLow

        const DensityOperatorLow
	densityOperator(const LazyDensityOperator<1>& m)
	{
	  size_t dim=m.getDimension();
	  DensityOperatorLow res(dim,dim);
  	  for (int i=0; i<dim; i++) for (int j=0; j<dim; j++) res(i,j)=m(i,j);
  	  return res;    
	}

      Then the partial density operator of any unary subsystem of a system of any rank may be calculated as follows (e.g. from a :class:`~quantumdata::StateVector`)::

        template<int RANK, int SUBSYSTEM> // for the subsystem indexed by the index SUBSYSTEM
        const DensityOperatorLow
	partialTraceOfUnarySubsystem(const StateVector<RANK>& psi)
	{
	  return partialTrace(psi,densityOperator,tmptools::Vector<SUBSYSTEM>(),DensityOperatorLow()); 
	}

    Calculating correlations for two harmonic-oscillator modes:
      Assume having the function ``calculateADaggerB`` as defined :ref:`above <calculateADaggerB>`. Then, if these modes are embedded at index positions 3 and 1 (note that the order might be important) in a larger system of arity larger than 3, we can write the following function::

        template<int RANK> // RANK>3
	const dcomp
	calculateADaggerB_atPositions3and1(const LazyDensityOperator<RANK>& matrix)
	{
	  return partialTrace(matrix,calculateADaggerB,tmptools::Vector<3,1>(),dcomp());
	}

      This will calculate :math:`\avr{a^\dag b}` (now :math:`b` and :math:`a` being the ladder operators of the modes at position 1 and 3, respectively) for the partial density operator of the embedded system (but without explicitly calculating this partial density operator).

    Accumulating an arbitrary ensemble of quantum averages:
      An arbitrary ensemble of real or complex numbers can be stored in a :c:macro:`TTD_DARRAY(1) <TTD_DARRAY>` or :c:macro:`TTD_CARRAY(1) <TTD_CARRAY>` (or even ``std::valarray``\ s), as these types fulfill all the requirements on type ``T``. The former was used to define the abstract interface :class:`structure::Averages`, one implementation being :class:`Mode`.

      E.g. a function for a harmonic-oscillator mode calculating :math:`\avr{a^\dag a}`, :math:`\avr{\lp a^\dag a\rp^2}`, and the real and imaginary parts of :math:`\avr{a}`, may be defined as ::

        typedef TTD_DARRAY(1) Averages;

	const Averages
	calculateModeAverages(const LazyDensityOperator<1>& matrix) const
	{
  	  Averages averages(4); averages=0;

  	  for (int n=1; n<int(matrix.getDimension()); n++) {
    	    double diag=matrix(n);
    	    averages(0)+=  n*diag;
    	    averages(1)+=n*n*diag;

    	    double sqrtn=sqrt(double(n));
    	    dcomp offdiag(matrix(n,n-1));
    	    averages(2)+=sqrtn*real(offdiag);
    	    averages(3)+=sqrtn*imag(offdiag);
	  }

  	  return averages;

	}

      Then the following function will calculate these averages for a mode embedded in a larger system at index position ``MODE_POSITION``::

        template<int RANK, int MODE_POSITION> // for the subsystem indexed by the index MODE_POSITION
        const Averages
	calculateEmbeddedHarmonicOscillatorAverages(const StateVector<RANK>& psi)
	{
	  return partialTrace(psi,calculateHarmonicOscillatorAverages,tmptools::Vector<MODE_POSITION>(),Averages()); 
	}



----------------------------
Assessing entanglement
----------------------------

.. py:module:: NegPT.h
   :synopsis: Defines the function calculating the negativity of a partially transposed density operator

...

