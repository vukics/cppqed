******************************************************************
A metaprogramming example
******************************************************************

In the following we analyse a metaprogramming example typical for the framework: how the compile-time vector ``0,3,2,6,4,5,1,9,8,7,10`` for the self-transposition :func:`above <blitzplusplus::basi::Transposer::transpose>` is prepared.

This is done by the following snippet in ``utils/include/impl/BlitzArraySliceIterator.tcc``:

.. literalinclude:: examples/multiarrayImplementationMetaProgramming.cc
  :language: c++
  :linenos:

