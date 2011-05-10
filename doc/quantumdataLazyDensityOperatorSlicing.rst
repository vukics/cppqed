:class:`~quantumdata::ldo::DiagonalIterator`---Notes on implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Difficulty: :class:`~quantumdata::LazyDensityOperator` is an abstract interface, and the organization of its data varies by implementation.

Implemented using a classical inheritance-based strategy idiom, together with both compile-time and run-time implementation selection (similarly to :ref:`blitz array slice iterator <basiImplementation>`, a special implementation (``DI_ImplSpecial``) is needed when the size of the compile-time vector ``V`` equals ``RANK``).

.. image:: figures/ldoDiagonalIterator.png
   :width: 643

The run-time polymorphy of :class:`~quantumdata::LazyDensityOperator` necessitates run-time implementation selection. This happens when the actual implementation (``DI_SV_Impl`` or ``DI_DO_Impl``, or their non-orthogonal counterparts) is chosen depending on whether the given :class:`~quantumdata::LazyDensityOperator` is a :class:`~quantumdata::StateVector` or a :class:`~quantumdata::DensityOperator` (or their :ref:`non-orthogonal counterparts <quantumdataNonOrthogonal>`).

.. class:: quantumdata::ldo::DiagonalIterator

  ``template <int RANK, typename V>`` (cf. :ref:`template parameters <quantumdataTemplates>`)

  Model of `InputIterator <http://www.cplusplus.com/reference/std/iterator/InputIterator/>`_, implemented with the help of ``boost::input_iterator_helper`` from `Boost.Operator <http://www.boost.org/doc/libs/1_44_0/libs/utility/operators.htm#iterator>`_.

  .. type:: value_type

    equivalent to ``const LazyDensityOperator<mpl::size<V>::value>``

  .. type:: Impl

    pointer to implementation::

      typedef boost::shared_ptr<typename mpl::eval_if_c<RANK==mpl::size<V>::value,
						        mpl::identity<details::DI_ImplSpecial<V> >,
						        mpl::identity<details::DI_Impl<RANK,V> > >::type> Impl;
