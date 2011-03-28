Implementing a harmonic-oscillator mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We demonstrate how to implement an element representing a pumped lossy mode in a truncated Fock space. It is described by the Hamiltonian:

.. math::

  H=\omega a^\dagger a+\lp\eta a^\dagger+\hermConj\rp


and Liouvillean:

.. math::

  \Liou\rho=\kappa\lp(n+1)\lp2a\rho a^\dagger-\comm{a^\dagger a}{\rho}_+\rp+n\lp2a^\dagger\rho a-\comm{a\,a^\dagger}{\rho}_+\rp\rp

The frequency-like parameters are :math:`\omega`, :math:`\kappa`, and :math:`\eta`, representing the mode frequency, loss rate, and pumping Rabi frequency, respectively. A dimensionless parameter is :math:`n`, the average number of thermal photons in the heat bath, and the cutoff of the Fock space.

The non-Hermitian Hamiltonian for the Monte Carlo wave-function method reads:

.. math::

  \HnH=\lp\omega-i\kappa\rp a^\dagger a+\lp\eta a^\dagger+\hermConj\rp

The element has to be represented by a class which inherits publicly from the necessary classes in the ``structure`` namespace. In this simple case, it is basically two helper functions returning :class:`quantumoperator::Tridiagonal`\ s, a constructor, and two virtual functions inherited from :class:`structure::ElementAveraged` that have to be written. Consider the file ``ExampleMode.h``:

.. literalinclude:: examples/ExampleMode.h
  :language: c++
  :linenos:
  :lines: 2-

This will suffice here. Let us look at the implementations in ``ExampleMode.cc``:

.. literalinclude:: examples/ExampleMode.cc
  :language: c++
  :linenos:

Lines 18-22:
  We construct the :class:`~structure::Free` base with the dimension of the system and the tuples for the frequency-like parameters of the system, which in this case are all complex (cf. :ref:`explanation <dynamicsBase>`). The tool from Boost.Assign facilitates the creation of lists of such tuples.

Lines 23-25:
  We construct the time-independent :class:`~structure::TridiagonalHamiltonian` base. This is greatly facilitated by the algebra and helpers of the :class:`quantumoperator::Tridiagonal` class.

Lines 26-29:
  We construct the :class:`~structure::ElementLiouvillean` base whose second template argument denotes the number of different quantum jumps, which is 2 in this case. The constructor takes the strategies for calculating the impact of a jump on a :type:`~structure::free::StateVectorLow`, and for calculating the probability from a :type:`~structure::free::LazyDensityOperator`. These strategy functions are produced from the free-standing helpers in Lines 10-14 through binding.

Lines 30-31:
  We construct the :class:`~structure::ElementAveraged` base, with parameters necessary to produce a simple key for quantum averages communicated towards the user. Here we calculate only three such averages, the expectation value of the number operator, and the real and imaginary parts of that of the ladder operator.

Lines 35-52:
  The inherited function :func:`~structure::Averaged::average` is implemented.

  Line 41:
    the expectation value of the photon number is calculated

  Lines 46-47:
    the expectation value of the ladder operator is calculated


The implementation of the helpers is also quite straightforward. It may come to a separate file ``ExampleModeImpl.cc``:

.. literalinclude:: examples/ExampleModeImpl.cc
  :language: c++
  :linenos:
