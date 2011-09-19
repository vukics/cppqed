******************************************************************
A metaprogramming example
******************************************************************

In the following we analyse a metaprogramming example typical for the framework: how the compile-time vector ``0,3,2,6,4,5,1,9,8,7,10`` for the self-transposition :func:`above <blitzplusplus::basi::Transposer::transpose>` is prepared.

This is done by the following snippet in ``utils/include/impl/BlitzArraySliceIterator.tcc``:

.. literalinclude:: examples/multiarrayImplementationMetaProgramming.cc
  :language: c++
  :linenos:

Line 13:
  We are using the `fold <http://www.boost.org/doc/libs/1_44_0/libs/mpl/doc/refmanual/fold.html>`_ meta-algorithm from Boost.MPL. It iterates over the :type:`sequence of ordinals <tmptools::Ordinals>` between ``0`` and ``RANK-1``.

Line 14:
  The initial state for the fold algorithm is an empty `compile-time vector of integers <http://www.boost.org/doc/libs/1_44_0/libs/mpl/doc/refmanual/vector-c.html>`_ and the `iterator <http://www.boost.org/doc/libs/1_44_0/libs/mpl/doc/refmanual/begin.html>`_ pointing to the first element of the compile-time vector ``V``. These two are "zipped" into a `compile-time pair <http://www.boost.org/doc/libs/1_44_0/libs/mpl/doc/refmanual/pair.html>`_. At the end, the first element of this pair will hold the result.

Lines 15-25
  express the forward operation of the fold algorithm in the form of a `compile-time lambda expression <http://www.boost.org/doc/libs/1_44_0/libs/mpl/doc/tutorial/handling-placeholders.html>`_. After each step, the new state will again be a pair composed of

  #. (Lines 15-20:) the vector `augmented <http://www.boost.org/doc/libs/1_44_0/libs/mpl/doc/refmanual/push-back.html>`_ either by the current element in ``V`` pointed to by the iterator (Line 17), or the current ordinal (Line 18), depending on whether ``V`` :class:`~tmptools::numerical_contains` this ordinal.

  #. (Lines 21-24:) the iterator is `advanced <http://www.boost.org/doc/libs/1_44_0/libs/mpl/doc/refmanual/next.html>`_ (Line 22) if the same numerical containment criterion is met, otherwise it is left untouched for the same element to be considered again (Line 23).

(Note that in TMP, "calling" ``numerical_contains<V,mpl::_2>`` twice is not a waste because the second time no further template instantiations are needed, the compiler will use the ones already instantiated the first time.)

Lines 32-34: the `first element <http://www.boost.org/doc/libs/1_44_0/libs/mpl/doc/refmanual/trivial-metafunctions-summary.html#first>`_ of the resulting pair of the above algorithm is picked. As it is easy to verify, ``TransposerMeta<11,tmptools::Vector<3,6,1,9,7> >::type`` will be an ``mpl::vector_c<int,0,3,2,6,4,5,1,9,8,7,10>``.
  
