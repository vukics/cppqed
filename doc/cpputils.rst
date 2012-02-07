.. _compile-time assertion: http://www.boost.org/doc/libs/1_48_0/libs/mpl/doc/refmanual/assert-msg.html

.. _cpputils:

***************************************************
Some further general-purpose modules from ``utils``
***************************************************

From ``utils`` what deserves documentation.

.. _cpputils_Evolved:

========================================
``Evolved`` and the underlying principle
========================================

.. class:: evolved::Evolved


============================
The ``trajectory`` namespace
============================

.. function:: void trajectory::runDt(trajectory::TrajectoryBase&, double time, double deltaT, bool displayInfo)

.. function:: void trajectory::run  (trajectory::Trajectory<A> &, double time, int          , bool displayInfo)

  template<typename A>


.. _cpputils_Parameters:

========
``Pars``
========

.. class:: parameters::ParameterTable

.. function:: void parameters::update(parameters::ParameterTable& table, int argc, char** argv, const std::string& mod="--")

======================
Multi-index iterator
======================

.. class:: cpputils::MultiIndexIterator


.. _cpputils_VFMSI:

=================================================
Iterating through rows/columns of multi-matrices
=================================================



==============================
Template metaprogramming tools
==============================


.. class:: tmptools::Range

  ``template <int N, int Nbeg>``, inherits publicly from ``boost::mpl::range_c<int,Nbeg,Nbeg+N>`` (Cf. `semantics <http://www.boost.org/doc/libs/1_48_0/libs/mpl/doc/refmanual/range-c.html>`_)

  Contains a `compile-time assertion`_ of ``N`` being nonnegative.

.. class:: tmptools::Ordinals

  ``template <int N>``, inherits publicly from :class:`~tmptools::Range`\ ``<N,0>``

.. class:: tmptools::numerical_contains

  ``template <typename Seq, typename ICW>``

  ``Seq`` is a Boost.MPL `forward sequence <http://www.boost.org/doc/libs/1_48_0/libs/mpl/doc/refmanual/forward-sequence.html>`_, and ``ICW`` is a compliant `integral constant wrapper <http://www.boost.org/doc/libs/1_48_0/libs/mpl/doc/refmanual/integral-constant.html>`_.

  This metaalgorithm differs from `boost::mpl::contains <http://www.boost.org/doc/libs/1_48_0/libs/mpl/doc/refmanual/contains.html>`_ in that it does not look for whether ``ICW`` as a *type* is contained in ``Seq``, but whether the *value* held (statically) by ``ICW`` is held by any element of the type sequence ``Seq``.


.. class:: tmptools::IsEvenAssert

  ``template <int N>``, inherits publicly from ``boost::mpl::int_<N/2>`` (Cf. `semantics <http://www.boost.org/doc/libs/1_48_0/libs/mpl/doc/refmanual/int.html>`_)

  Contains a `compile-time assertion`_ of ``N`` being an even number. 


.. class:: tmptools::Vector

  ``template <int V0=TMPTOOLS_VECTOR_DEFAULT_ARG, ...>``


==============
Linear algebra
==============

.. class:: linalg::VectorSpace

  ``template <typename T, typename B>``



===========================
Further extensions to Blitz
===========================

.. py:module:: BlitzArrayExtensions.h
   :synopsis: Defines some generic extensions of ``blitz::Array`` functionality.
   

.. class:: blitzplusplus::NonContiguousStorageException

   inherits publicly from :class:`cpputils::Exception`

.. function:: const blitz::Array<T,1> blitzplusplus::unaryArray(const blitz::Array<T,RANK>& array)

  ``template <typename T, int RANK>``

  Returns a unary view of ``array``.

  This is meant to be used only if the underlying storage is contiguous, because while a multi-rank array may be able to represent a view of memory of some more or less intricate structure pattern (e.g. slices), a unary array is not capable of this. In debug mode, violation is detected at runtime via the ``blitz::Array`` class's ``isStorageContiguous`` member function and an exception of type :class:`~blitzplusplus::NonContiguousStorageException` is thrown.

  The returned array does not take part in the reference-counting mechanism of ``blitz::Array``, therefore, it does not own its data, and …

  .. warning::

    … the returned array has no means of ensuring that the referenced data still exists.

.. class:: blitzplusplus::BinaryArrayOrderingErrorException

  inherits publicly from :class:`cpputils::Exception`

.. function:: const blitz::Array<T,2> blitzplusplus::binaryArray(const blitz::Array<T,TWO_TIMES_RANK>& array)

  ``template <typename T, int TWO_TIMES_RANK>``

  Returns a binary view of ``array``. ``TWO_TIMES_RANK`` must be an even number, violation is detected at compile time by :class:`tmptools::IsEvenAssert`.

  The same requirement of contiguity an the same warning applies as :func:`above <blitzplusplus::unaryArray>`, and in addition, further assumptions on the storage order must be made: The storage of the two multi-indeces must not be intertwined and must be layed out in the same way, so that e.g. for ``RANK=4``, the member function ``array.ordering()`` should return an octary tiny vector like::

    <1 3 2 0 | 5 7 6 4>

  Violation is detected at runtime, and an exception of type :class:`~blitzplusplus::BinaryArrayOrderingErrorException` is thrown.

.. py:module:: ComplexArrayExtensions.h
   :synopsis: Defines some generic extensions for complex ``blitz::Array``\ s.

.. function:: const TTD_CArray<RANK1__PL__RANK2> blitzplusplus::doDirect(const TTD_CArray<RANK1>& array1, const TTD_CArray<RANK2>& array2, boost::mpl::bool_<MULT> tag)

  ``template <int RANK1, int RANK2, bool MULT>``

  Returns the direct product (if ``tag`` is a ``true`` boolean constant wrapper) :math:`A_{i,j}=A1_i*A2_j`, or direct sum (otherwise) :math:`A_{i,j}=A1_i+A2_j` of ``array1`` and ``array2``, with :math:`i,\ j` running through all the multi-indices. The implementation is in terms of :func:`unaryArray <blitzplusplus::unaryArray>` views of the arguments. 




.. py:module:: BlitzTinyExtensions.h
   :synopsis: Defines some generic extensions for ``blitz::TinyVector``\ s.


.. function:: blitz::TinyVector<T1,RANK1__PL__RANK2> blitzplusplus::concatenateTinies(const blitz::TinyVector<T1,RANK1>& tiny1, const blitz::TinyVector<T2,RANK2>& tiny2)

  ``template <typename T1, typename T2, int RANK1, int RANK2>``

  Concatenates ``tiny1`` and ``tiny2`` with the help of the compile-time–runtime facility ``boost::mpl::for_each``. ``T2`` must be convertible to ``T1``.


.. function:: blitz::TinyVector<T,TWO_TIMES_RANK__PE__2> halfCutTiny(const blitz::TinyVector<T,TWO_TIMES_RANK>& tiny)

  ``template<typename T, int TWO_TIMES_RANK>``


