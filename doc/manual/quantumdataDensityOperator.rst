Density operator
^^^^^^^^^^^^^^^^^

.. py:module:: DensityOperator.h
   :synopsis: Defines DensityOperator in namespace quantumdata and corresponding operations

.. class:: quantumdata::DensityOperator

  :ref:`template parameters <quantumdataTemplates>`: RANK; inherits publicly from :class:`~quantumdata::LazyDensityOperator`, and privately from :class:`~quantumdata::ArrayBase`, and also from :class:`linalg::VectorSpace` which adds a lot of free-standing arithmetic functions.

  .. note::

    A ``DensityOperator<RANK>`` represents a density operator on a Hilbert space of arity ``RANK``. This makes that the number of its indeces is actually ``2*RANK``. This is the reason why it inherits from ``ArrayBase<2*RANK>``.

  The interface is similar to :class:`~quantumdata::StateVector` with self-evident differences. Here only the most important of these will be touched upon:

  .. function:: explicit DensityOperator(const StateVector<RANK>& psi)

    Constructs the class as a dyadic product of ``psi``.

  .. function:: double norm() const

  .. function:: double renorm()

    Both functions return the trace "norm", but the latter one also renormalizes.

  .. function:: const linalg::CMatrix matrixView() const

    (also in non-constant version) Returns a two-dimensional view of the underlying data.

  .. function:: const dcomp operator()(const Idx& i, const Idx& j)

    This function implements the virtual indexing function :func:`quantumdata::LazyDensityOperator::operator()` in a trivial way, simply by accessing the necessary element in memory::

      const dcomp operator()(const Idx& i, const Idx& j) const {return operator()()(blitzplusplus::concatenateTinies<int,int,RANK,RANK>(i,j));}

