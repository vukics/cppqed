.. _generalElements:

===============================
Some general-purpose elements
===============================

-------------------------------------------
Polarizable particles in optical fields
-------------------------------------------

.. _generalElements_Mode:

Mode
^^^^^^

.. math::

  H_\text{mode}&=\omega\adagger a+i\lp\eta\adagger-\hermConj\rp\\
  &\equiv-iz\adagger a+i\lp\eta\adagger-\hermConj\rp\equiv(-\delta-i\kappa)\adagger a+i\lp\eta\adagger-\hermConj\rp

where we have introduced

.. math::

  z\equiv\kappa-i\delta

.. type:: mode::SmartPtr

  ::

    typedef boost::shared_ptr<const ModeBase> SmartPtr;


.. class:: Mode


.. function:: const mode::StateVector mode::coherent(const dcomp& alpha, size_t cutoff)

.. function:: const mode::StateVector mode::init(const mode::Pars&)

.. function:: const mode::SmartPtr mode::maker(const mode::Pars           &, QM_Picture, const A&)

.. function:: const mode::SmartPtr mode::maker(const mode::ParsLossy      &, QM_Picture, const A&)

.. function:: const mode::SmartPtr mode::maker(const mode::ParsPumped     &, QM_Picture, const A&)

.. function:: const mode::SmartPtr mode::maker(const mode::ParsPumpedLossy&, QM_Picture, const A&)

  template<typename A>


.. todo::

   Clarify mfNKX_AbsSqr





.. _generalElements_Particle:

Particle in momentum space
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. class::Spatial

  discuss intimacies of discrete Fourier transform...

.. math::

  H_\text{kinetic}=\frac{p^2}{2m}\equiv\omrec k^2


.. todo::

   Particle::Averaged is wasteful since it FFTs a matrix even in the case when the LazyDensityOperator it receives is in fact a StateVector. Solution: implement an fft in LazyDensityOperator depending on its "origin". 

Interactions
^^^^^^^^^^^^^^

.. class:: JaynesCummings

.. class:: ParticleTwoModes2D


script: 1particle1mode
^^^^^^^^^^^^^^^^^^^^^^^




z cavity axis, x orthogonal direction

g pump mode function, f cavity mode function

1

.. math::

  H_\text{kinetic}^x+\sqrt{U_0\vclass}\lp g(x)\adagger+\hermConj\rp+H_\text{mode}


2

.. math::

  H_\text{kinetic}^z+U_0\abs{f(z)}^2\adagger a+\sqrt{U_0\vclass}\lp f(z)\adagger+\hermConj\rp+H_\text{mode}


3

.. math::

  H_\text{kinetic}^z+\vclass\abs{g(z)}^2+U_0\abs{f(z)}^2\adagger a+\sqrt{U_0\vclass}\lp f(z)g(z)\adagger+\hermConj\rp+H_\text{mode}


4

.. math::

  H_\text{kinetic}^z+\vclass\abs{g(z)}^2+U_0\abs{f(z)}^2\adagger a+H_\text{mode}


---------------
``MultiLevel``
---------------

.. class:: MultiLevel

  Define multi-level systems (e.g. atoms with arbitrary level schemes) with various driving and loss schemes at compile time.

------------------
Other
------------------

.. todo::

   Implement a bosonic many-body system with an arbitrary number of modes (the number known at compile time). The mode operators can be sparse matrices. This needs the system of extended quantum operators.
