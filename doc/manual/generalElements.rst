.. _generalElements:

===============================
Some general-purpose elements
===============================

-------------------------------------------
Polarizable particles in optical fields
-------------------------------------------


Mode
^^^^^^

.. math::

  H_\text{mode}&=\omega\adagger a+i\lp\eta\adagger-\hermConj\rp\\
  &\equiv-iz\adagger a+i\lp\eta\adagger-\hermConj\rp\equiv(-\delta-i\kappa)\adagger a+i\lp\eta\adagger-\hermConj\rp


.. class:: Mode

.. todo::

   Clarify mfNKX_AbsSqr


Particle
^^^^^^^^^

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



------------------
Other
------------------

.. todo::

   Implement a bosonic many-body system with an arbitrary number of modes (the number known at compile time). The mode operators can be sparse matrices. This needs the system of extended quantum operators.
