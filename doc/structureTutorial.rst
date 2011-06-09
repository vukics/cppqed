****************************************
Implementing a harmonic-oscillator mode
****************************************

We demonstrate how to implement an element representing a pumped lossy mode in a truncated Fock space. It is described by the Hamiltonian:

.. math::
  :label: harmonicOscillatorHamiltonian

  H=-\delta a^\dagger a+\lp\eta a^\dagger+\hermConj\rp,


where :math:`\delta` is the detuning between the oscillator and the pump frequencies. The Liouvillean:

.. math::

  \Liou\rho=\kappa\lp(n+1)\lp2a\rho a^\dagger-\comm{a^\dagger a}{\rho}_+\rp+n\lp2a^\dagger\rho a-\comm{a\,a^\dagger}{\rho}_+\rp\rp

The frequency-like parameters are :math:`\delta`, :math:`\kappa`, and :math:`\eta`, representing the mode frequency, loss rate, and pumping Rabi frequency, respectively. A dimensionless parameter is :math:`n`, the average number of thermal photons in the heat bath, and the cutoff of the Fock space.

Using the notation of Sec. :ref:`MCWF_method`,

.. math::

  J_0=\sqrt{2\kappa(n+1)}a

.. math::

  J_1=\sqrt{2\kappa n}a^\dagger

The non-Hermitian Hamiltonian for the Monte Carlo wave-function method reads:

.. math::
  :label: harmonicOscillatorNonHermitianHamiltonian

  \HnH=\lp-\delta-i\kappa(2n+1)\rp a^\dagger a+\lp\eta a^\dagger+\hermConj\rp

The element has to be represented by a class which inherits publicly from the necessary classes in the ``structure`` namespace. In this simple case, it is basically two helper functions returning :class:`quantumoperator::Tridiagonal`\ s, a constructor, and two virtual functions inherited from :class:`~structure::ElementAveraged` that have to be written. Consider the file ``ExampleMode.h``:

.. literalinclude:: examples/ExampleMode.h
  :language: c++
  :linenos:
  :lines: 2-5, 7-28

This will suffice here. Let us look at the implementations in ``ExampleMode.cc``:

.. literalinclude:: examples/ExampleMode.cc
  :language: c++
  :linenos:
  :lines: 1-54

Lines 18-22:
  We construct the :class:`~structure::Free` base with the dimension of the system and the tuples for the frequency-like parameters of the system, which in this case are all complex (cf. :ref:`explanation <dynamicsBase>`). The tool from Boost.Assign facilitates the creation of lists of such tuples.

Lines 23-25:
  We construct the time-independent :class:`~structure::TridiagonalHamiltonian` base. This is greatly facilitated by the algebra and helpers of the :class:`quantumoperator::Tridiagonal` class.

.. warning::

  When implementing the Hamiltonian, not :math:`H` itself but :math:`\frac Hi` has to supplied!

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


===================
Interaction picture
===================

In many situations, it pays to transfer to interaction picture defined by the first term of the Hamiltonian :eq:`harmonicOscillatorHamiltonian`. The ladder operator in interaction picture is :math:`a\Int(t)=e^{i\delta t}a`, so that the Hamiltonian reads

.. math::

  H\Int(t)=-i\kappa(2n+1) a^\dagger a+\lp\eta a^\dagger e^{-i\delta t}+\hermConj\rp

In this case, the class representing the element has to be derived from :class:`~structure::Exact` as well, which represents the transformation between the two pictures. In addition, instead of :class:`~structure::TridiagonalHamiltonian`\ ``<1,false>``, we need to derive from :class:`~structure::TridiagonalHamiltonian`\ ``<1,true>``, because the Hamiltonian is now time-dependent.

.. note:: 

  The framework requires that the jump and the averages be calculated in the normal picture also in this case (cf. explanation of classes :class:`structure::Hamiltonian`, :class:`structure::Exact`, :class:`structure::Liouvillean`, and :class:`structure::Averaged`, and furthermore Sec. :ref:`MCWF_Trajectory`). This allows for reusing the same code in both pictures. (Incidentally, here the Liouvillean remains unchanged anyway.)


.. literalinclude:: examples/ExampleMode.h
  :language: c++
  :linenos:
  :lines: 2-11, 30-


In the implementation, the only difference from the previous case will be the constructor which now also requires a :type:`~structure::free::Frequencies` object, and the virtual function :func:`~structure::FreeExact::updateU`.

.. literalinclude:: examples/ExampleMode.cc
  :language: c++
  :linenos:
  :lines: 57-83

:class:`~structure::FreeExact` assumes that the operator transforming between the two pictures is diagonal, and the factors to update are simply its diagonal elements. If this is not the case, :class:`~structure::Exact` has to be used instead.

The construction of the :type:`~structure::free::Frequencies` object is also straightforward:

.. literalinclude:: examples/ExampleModeImpl.cc
  :language: c++
  :linenos:
  :lines: 58-62

.. note::

  Since a lot of the code from the previous case can be reused here, one will usually adopt an inheritence- or class-composition-based solution to implement classes like ``PumpedLossyMode`` and ``PumpedLossyModeIP`` (cf. :ref:`generalElements_Mode`).




*******************************
Implementing an X-X interaction
*******************************

Let us consider the interaction described by the Hamiltonian

.. math::

  H_\text{X-X}=g(a+a^\dagger)(b+b^\dagger)

The class implementing this interaction has to be derived from :class:`~structure::Interaction`\ ``<2>`` because it is a binary interaction, and :class:`~structure::TridiagonalHamiltonian`\ ``<2,...>`` (note that :class:`quantumoperator::Tridiagonal` is capable to represent direct products of tridiagonal matrices).

The only thing requiring some care is that once we transform some elements into interaction picture, the whole Hamiltonian is transformed, that is, :math:`a` or :math:`b` or both may be in interaction picture. Here, for the sake of simplicity, we assume that both constituents are of the type ``PumpedLossyModeIP``.

Consider ``ExampleInteraction.h``:

.. literalinclude:: examples/ExampleInteraction.h
  :language: c++
  :linenos:
  :lines: 2-

``ExampleInteraction.cc`` then reads

.. literalinclude:: examples/ExampleInteraction.cc
  :language: c++
  :linenos:

As we see, the Hamiltonian can be written in a rather straightforward way, and the corresponding :class:`quantumoperator::Frequencies` objects are also multiplied, expressing the fact that both constituents are in interaction picture.
