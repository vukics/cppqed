.. _quantumtrajectory:

===================================
The ``quantumtrajectory`` namespace
===================================

.. image:: figures/trajectory.png
   :width: 801
   
   
.. todo::

   Apply Boost.Parameter for facilitating call of functions with complex signature, e.g. constructors of complex classes, especially in cases like in the trajectory bundle where many sensible default parameter values could be defined. (Eg. it is extremely tedious scaleAbs needs to be specified each time.) This is in general most useful for constructors of complex classes.


.. function:: void evolve(quantumdata::StateVector<RANK>&, const structure::QuantumSystem<RANK>&, const ParsEvolution&, V)

  ``template<int RANK, typename V>``

.. _MCWF_Trajectory:

--------------------
MCWF trajectory
--------------------

In the framework, a single :ref:`Monte-Carlo wave function step <MCWF_method>` at time :math:`t` (at which point the Schrödinger and interaction pictures coincide) is implemented as a sequence of the following stages:

1. If the system time evolution has Hamiltonian part, it is evolved with an adaptive-size step (cf. :ref:`Evolved <cpputils_Evolved>`). This takes the system into :math:`t+\Delta t`.

2. The exact part (if any) of the time evolution is applied, making that the Schrödinger and interaction pictures coincide again at :math:`t+\Delta t`.

3. The state vector is renormalized.

4. If the system is Liouvillean, the possibility of a quantum jump is considered:

  #. The rates (probabilities per unit time) corresponding to all jump operators are calculated. If some rates are found negative ("special jump", cf. explanation at :func:`~structure::Liovillean::probability`), then :math:`J_\text{at}\ket\Psi` is calculated (and tabulated) instead, and the probability is calculated as :math:`\delta r_\text{at}=\norm{J_\text{at}\ket\Psi}^2`.

  #. The total jump rate :math:`\delta r` is calculated.

  #. It is randomly decided which (if any) of the jumps to perform. If it is found to be a special jump, then the tabulated :math:`J_\text{at}\ket\Psi` is taken.

5. Time-step management is performed: the adaptive-stepsize ODE stepper gives a guess for the next timestep (:math:`\Delta t_\text{next}`). If the probability :math:`\delta r\Delta t_\text{next}` is found to overshoot a given limit (usually 0.1), then the next timestep is decreased.


.. class:: quantumtrajectory::MCWF_Trajectory

  ``template <RANK>``

  .. function:: void step(double deltaT) const

Alternative probability calculation based on tabbed jumped state vectors signalled by negative jump probabilities.


--------------------------
Ensembles of trajectories
--------------------------

.. class:: quantumtrajectory::EnsembleMCWF

---------------------------
Master equation evolution
---------------------------

.. math::
  :label: masterEqInTermsOfMCWF

  \dot\rho=\frac1{i\hbar}\comm{H}{\rho}+\sum_m\lp{J_m\rho{J_m^\dag}-\frac12\comm{J_m^\dag J_m}{\rho}_+}\rp\equiv\frac1{i\hbar}\lp\HnH\rho-\rho\HnH^\dag\rp+\sum_mJ_m\rho{J_m^\dag}=2\Re\lbr\frac\HnH{i\hbar}\rho\rbr+\sum_mJ_m\lp{J_m\rho}\rp^\dag


.. namespace:: quantumtrajectory
.. class:: quantumtrajectory::Master


Alternative Liouvillean calculation based on additional virtual function signalled by something.

On the basis of Pictures.pdf find out when exactly the code reusal for the calculation of the Liouvillean can be applied.


Performance profile
^^^^^^^^^^^^^^^^^^^^

./release/1particle1mode --evol master --dc 1 --T 1 --fin 7



======================== =============== =========================
operation                n/timestep      time
======================== =============== =========================
Hamiltonian              (nRejected+1)*5 1.31s
TwoTimesRealPartOfSelf   "               0.05s
Liouvillean              (")*nJump       0.33s
Exact                    1               0.40s
Smoothing                1               0.10s
Averaging                ≤1              0.04s
======================== =============== =========================
