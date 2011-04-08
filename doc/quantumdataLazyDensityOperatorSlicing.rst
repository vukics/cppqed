:class:`~quantumdata::ldo::DiagonalIterator`---Notes on implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Difficulty: :class:`~quantumdata::LazyDensityOperator` is an abstract interface, and the organization of its data varies by implementation.

.. class:: quantumdata::ldo::DiagonalIterator

  ``template <int RANK, typename V>`` (:ref:`template parameters <quantumdataTemplates>`)

  Model of `InputIterator <http://www.cplusplus.com/reference/std/iterator/InputIterator/>`_.

  .. type:: value_type

    equivalent to ``const LazyDensityOperator<mpl::size<V>::value>``

Implemented using a classical inheritance-based strategy idiom, together with some compile-time implementation selection.

.. image:: figures/ldoDiagonalIterator.png
   :width: 643

...
