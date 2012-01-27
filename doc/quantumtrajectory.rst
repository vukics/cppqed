.. _quantumtrajectory:

===================================
The ``quantumtrajectory`` namespace
===================================

.. image:: figures/trajectory.png
   :width: 801
   

.. _MCWF_Trajectory:

--------------------
MCWF trajectory
--------------------



In the framework, a single :ref:`Monte-Carlo wave function step <MCWF_method>` at time :math:`t` (at which point the Schrödinger- and interaction pictures coincide) is implemented as a sequence of the following stages:

I. Coherent time development is applied:

  #. If the system time evolution has Hamiltonian part, it is evolved with an adaptive-size step (cf. :ref:`Evolved <cpputils_Evolved>`). This takes the system into :math:`t+\Delta t`.

  #. The exact part (if any) of the time evolution is applied, making that the Schrödinger and interaction pictures coincide again at :math:`t+\Delta t`.

  #. The state vector is renormalized.

II. If the system is *not* Liouvillean, the timestep ends here, reducing to a simple ODE evolution. Otherwise:

  #. The rates (probabilities per unit time) corresponding to all jump operators are calculated. If some rates are found negative ("special jump", cf. explanation at :func:`~structure::Liouvillean::probabilities`), then :math:`J_\text{at}\ket\Psi` is calculated (and tabulated) instead, and the probability is calculated as :math:`\delta r_\text{at}=\norm{J_\text{at}\ket\Psi}^2`.

  #. First, it is verified whether the total jump probability is not too big. This is performed on two levels:

    a. The total jump rate :math:`\delta r` is calculated.

    b. If :math:`\delta r\delta t>\delta p_\text{limit}'`, the step is retraced: both the state vector and the state of the ODE stepper are restored to cached values at the beginning of the timestep, and phase I. is performed anew with a smaller stepsize :math:`\delta p_\text{limit}/\delta r`. With this, we ensure that :math:`\delta p_\text{limit}` is likely not to be overshot in the next try.

      .. note:: It is assumed that :math:`\delta p_\text{limit}'>\delta p_\text{limit}`, their ratio being a parameter of the MCWF stepper.

    c. If just :math:`\delta r\delta t_\text{next}>\delta p_\text{limit}` (where :math:`\Delta t_\text{next}` is a guess for the next timestep given by the ODE stepper), the coherent step is accepted, but the timestep to try next is modified, to reduce the likeliness of overshoot: :math:`\delta t_\text{next}\longrightarrow\delta p_\text{limit}/\delta r`.

    .. seealso:: The discussion at :ref:`MCWF_method_adaptive`.

  3. After a successful coherent step resulting in an acceptable :math:`\delta r`, the possible occurence of a quantum jump is considered: It is randomly decided which (if any) of the jumps to perform. If it is found to be a special jump, then the tabulated :math:`J_\text{at}\ket\Psi` is taken.


.. note:: In phase II.2.b., another approach would be not to trace back the whole step, but make a coherent step *backwards* to an intermediate time instant found by linear interpolation. This has several drawbacks, however, the most significant being that in the ODE stepper, it is not clear what to take as the timestep to try at the point when the direction of time is reversed. (Although in :class:`Evolved` it is simply taken to be the timestep done in the last step…)


.. class:: quantumtrajectory::MCWF_Trajectory

  ``template <RANK>``

  .. function:: void step(double deltaT) const


--------------------------
Ensembles of trajectories
--------------------------

.. class:: quantumtrajectory::EnsembleMCWF

  ...

---------------------------
Master equation evolution
---------------------------

.. math::
  :label: masterEqInTermsOfMCWF

  \dot\rho=\frac1{i\hbar}\comm{H}{\rho}+\sum_m\lp{J_m\rho{J_m^\dag}-\frac12\comm{J_m^\dag J_m}{\rho}_+}\rp\equiv\frac1{i\hbar}\lp\HnH\rho-\rho\HnH^\dag\rp+\sum_mJ_m\rho{J_m^\dag}=2\Re\lbr\frac\HnH{i\hbar}\rho\rbr+\sum_mJ_m\lp{J_m\rho}\rp^\dag


.. namespace:: quantumtrajectory

.. class:: quantumtrajectory::Master

  ...

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


---------------------------------
Dispatcher
---------------------------------

.. function:: evolve(quantumdata::StateVector<RANK>& psi, const structure::QuantumSystem<RANK>& sys, const ParsEvolution& pe, V v)

  ``template<int RANK, typename V>``

