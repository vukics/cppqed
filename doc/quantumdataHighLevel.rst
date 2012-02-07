.. rubric:: ``ArrayBase``

.. py:module:: ArrayBase.h
   :synopsis: Defines ArrayBase in namespace quantumdata

.. class:: quantumdata::ArrayBase

  ``template <int RANK>`` (cf. :ref:`template parameters <quantumdataTemplates>`), inherits privately from ``boost::noncopyable``.

  This class comprises the common functionalities of :class:`~quantumdata::StateVector` and :class:`~quantumdata::DensityOperator`.

  .. type:: ArrayLow

    ::

      typedef TTD_CArray<RANK> ArrayLow;

  .. type:: CVector

    ::

      typedef linalg::CVector CVector;

    Equivalent to :class:`TTD_CArray`\ ``<1>``.

  .. function:: explicit ArrayBase(const ArrayLow& arrayLow)

    By-reference semantics (basically the copy of a ``blitz::Array``). Otherwise, copying is not possible.

  .. function:: ArrayBase& operator=(const ArrayLow& arrayLow)
  .. function:: ArrayBase& operator=(const ArrayBase& arrayBase)

    By-value semantics (like the assignment of a ``blitz::Array``). The latter synthesised by the compiler.

  .. function:: const ArrayLow& operator()() const
 
    (also in non-const version) Returns the underlying :type:`ArrayLow`.

  .. function:: void operator+=(const ArrayBase& arrayBase)
  .. function:: void operator-=(const ArrayBase& arrayBase)

    Naive vector-space operations.

  .. function:: void operator*=(const OTHER& dc)
  .. function:: void operator/=(const OTHER& dc)

    ``template <typename OTHER>``

    Naive vector-space operations allowing also for mixed-mode arithmetic.

  .. function:: const CVector vectorView() const

    (also in non-constant version) Returns a one-dimensional view of the underlying data, created on the fly via :func:`blitzplusplus::unaryArray`.

    .. note::

      This is meant to be used only if the underlying storage is contiguous, which may not be the case since the object may reference :type:`ArrayLow`\ s of any layout. In debug mode, a non-contiguous storage is detected by the implementing function :func:`blitzplusplus::unaryArray`, and an exception of type :class:`blitzplusplus::NonContiguousStorageException` is thrown.

  .. function:: double frobeniusNorm() const

    calculates the "entrywise" norm

    .. math::

      \norm{A}=\sqrt{\sum_i \abs{A_i}^2}

    with :math:`i` running through all the multi-indices.

*************
State vector
*************

.. py:module:: StateVector.h
   :synopsis: Defines StateVector in namespace quantumdata and corresponding operations

.. class:: quantumdata::StateVector

  ``template <int RANK>`` (cf. :ref:`template parameters <quantumdataTemplates>`); inherits publicly from :class:`~quantumdata::LazyDensityOperator`\ ``<RANK>``, and privately from :class:`~quantumdata::ArrayBase`\ ``<RANK>``, and also from :class:`linalg::VectorSpace`\ ``<StateVector<RANK> >`` which adds a lot of free-standing arithmetic functions.

  .. rubric:: Types

  .. type:: Dimensions

    (inherited from :class:`DimensionsBookkeeper`)

  .. type:: StateVectorLow

    same as :type:`quantumdata::Types::StateVectorLow`
  
  .. type:: DensityOperatorLow

    same as :type:`quantumdata::Types::DensityOperatorLow`

  .. type:: Idx

    (inherited from :class:`~quantumdata::LazyDensityOperator`)

  .. rubric:: Constructors, assignment

  .. function:: StateVector(const StateVectorLow& psi, ByReference)

    Constructs the class in such a way that the underlying data will reference the same data as ``psi``.

    .. note:: 

      Since everywhere else the class represents by-value semantics, some care is needed with this constructor's use. For this reason, the tagging dummy class ``ByReference`` is introduced, to make the user conscious of what the semantics is.

  .. function:: explicit StateVector(const Dimensions& dim, bool init=true)

    Constructs the class with a newly allocated chunk of memory, which is initialized only if ``init`` is ``true``.

  .. function:: StateVector(const StateVector& psi)

    Copy constructor using by value semantics, that is, deep copy.

  .. function:: StateVector(const StateVector<RANK2>& psi1, const StateVector<RANK__MI__RANK2>& psi2)

    ``template <int RANK2>``

    Constructs the class as the direct product of ``psi1`` and ``psi2``, whose arities add up to ``RANK``.

    .. seealso:: The implementation relies on :func:`blitzplusplus::concatenateTinies` and :func:`blitzplusplus::doDirect`.

  .. function:: StateVector& operator=(const StateVector& sv)

  .. function:: StateVector& operator=(const OTHER& other)

    ``template <typename OTHER>`` ::

      template<typename OTHER> StateVector& operator=(const OTHER& other) {operator()()=other; return *this;}

    The standard assignment and the templated assignment together cover a lot of possibilities, including also assignment from a StateVectorLow, but for example also from a :class:`TTD_DArray`\ ``<RANK>``, or just a const c-number.

    By-value semantics.

  .. rubric:: Class-specific functionality

  .. function:: const linalg::CVector vectorView() const

    (inherited from :class:`~quantumdata::ArrayBase`; also in non-constant version) Returns a one-dimensional view of the underlying data.

    .. note::

      The same note applies as for :func:`quantumdata::ArrayBase::vectorView`, and here especially because of the presence of the referencing constructor. Errors are detected at runtime in debug mode.

  .. function:: const StateVectorLow& operator()() const

    (inherited from :class:`~quantumdata::ArrayBase`; also in non-constant version) Returns the underlying ``blitz::Array`` storage.

  .. function:: double norm() const

  .. function:: double renorm()

    Both functions return the norm :math:`\norm{\Psi}`, but the latter one also renormalizes. Implemented in terms of :func:`quantumdata::ArrayBase::frobeniusNorm`.

  .. function:: const DensityOperatorLow dyad(const StateVector& psi) const

  .. function:: const DensityOperatorLow dyad() const

    Both functions form a dyad, the second one with the same object::

      const DensityOperatorLow dyad() const {return dyad(*this);}

    This is a rather expensive operation, implemented in terms of :func:`blitzplusplus::doDirect`.

  .. function:: void addTo(DensityOperator<RANK>& densityOperator)

    This function adds a dyad of the present object to ``densityOperator``, without actually forming the dyad in memory (so that this is not implemented in terms of :func:`~quantumdata::StateVector::dyad`). This is important in situations when an average density operator is needed from an ensemble of state vectors, an example being :class:`quantumtrajectory::EnsembleMCWF`.

  .. rubric:: Linear algebra

  .. function:: StateVector& operator+=(const StateVector& psi)
  
  .. function:: StateVector& operator-=(const StateVector& psi)

  .. function:: const StateVector operator-() const
  
  .. function:: const StateVector operator+() const

  .. function:: StateVector& operator*=(const OTHER& dc)

  .. function:: StateVector& operator/=(const OTHER& dc)

    ``template <typename OTHER>``

    These are vector-space operations implemented in a naive way (as opposed to e.g. the expression-template mechanism of Blitz), the last two being templated to allow for mixed-mode arithmetics.

  .. rubric:: Implementing the :class:`~quantumdata::LazyDensityOperator` interface

  .. function:: const dcomp operator()(const Idx& i, const Idx& j)

    This function implements the virtual indexing function :func:`quantumdata::LazyDensityOperator::operator()` in a dyadic-product way::

      const dcomp operator()(const Idx& i, const Idx& j) const
      {
        return operator()()(i)*conj(operator()()(j));
      }


.. rubric:: Free-standing helpers

The inheritance of :class:`~quantumdata::StateVector` from :class:`linalg::VectorSpace` provides for a lot of free-standing helpers describing vector-space algebra. These are all naively based on the arithmetic member functions like :func:`~quantumdata::StateVector::operator+=`, :func:`~quantumdata::StateVector::operator*=`, etc.

There are two further free-standing helpers:

.. function:: const StateVector<RANK1__PL__RANK2> quantumdata::operator*(const StateVector<RANK1>& psi1, const StateVector<RANK2>& psi2)

  ``template <int RANK1, int RANK2>`` (cf. :ref:`template parameters <quantumdataTemplates>`)

  This function creates the direct product, relying on the direct-product constructor.

.. function:: const dcomp quantumdata::braket(const StateVector<RANK>& psi1, const StateVector<RANK>& psi2)

  ``template <int RANK>`` (cf. :ref:`template parameters <quantumdataTemplates>`)

  Calculates the inner product, relying on :func:`~quantumdata::StateVector::vectorView`.


****************
Density operator
****************

.. py:module:: DensityOperator.h
   :synopsis: Defines DensityOperator in namespace quantumdata and corresponding operations

:class:`~quantumdata::DensityOperator`\ s interface is similar to :class:`~quantumdata::StateVector`\ s with obvious differences. Here only the most important will be tackled:

.. class:: quantumdata::DensityOperator

  ``template <int RANK>`` (cf. :ref:`template parameters <quantumdataTemplates>`); inherits publicly from :class:`~quantumdata::LazyDensityOperator`\ ``<RANK>``, and privately from :class:`~quantumdata::ArrayBase`\ ``<2*RANK>``, and also from :class:`linalg::VectorSpace`\ ``<DensityOperator<RANK> >`` which adds a lot of free-standing arithmetic functions.

  .. note::

    A :class:`~quantumdata::DensityOperator`\ ``<RANK>`` represents a density operator on a Hilbert space of arity ``RANK``. This makes that the number of its indeces is actually ``2*RANK``. This is the reason why it inherits from :class:`~quantumdata::ArrayBase`\ ``<2*RANK>``.

  .. function:: explicit DensityOperator(const StateVector<RANK>& psi)

    Constructs the class as a dyadic product of ``psi``.

  .. function:: double norm() const

  .. function:: double renorm()

    Both functions return the trace "norm", but the latter one also renormalizes.

  .. function:: const linalg::CMatrix matrixView() const

    (also in non-constant version) Returns a two-dimensional view of the underlying data, created on the fly via :func:`blitzplusplus::binaryArray`.

  .. function:: const dcomp operator()(const Idx& i, const Idx& j)

    This function implements the virtual indexing function :func:`quantumdata::LazyDensityOperator::operator()` in a trivial way, simply by accessing the necessary element in memory (it must rely on :func:`blitzplusplus::concatenateTinies`, though)::

      const dcomp operator()(const Idx& i, const Idx& j) const 
      {
        return operator()()(blitzplusplus::concatenateTinies<int,int,RANK,RANK>(i,j));
      }

