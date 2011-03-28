.. _quantumtrajectory:

===================================
The ``quantumtrajectory`` namespace
===================================

.. image:: figures/trajectory.png
   :width: 801
   
   
.. todo::

   Apply Boost.Parameter for facilitating call of functions with complex signature, e.g. constructors of complex classes, especially in cases like in the trajectory bundle where many sensible default parameter values could be defined. (Eg. it is extremely tedious scaleAbs needs to be specified each time.) This is in general most useful for constructors of complex classes.


.. function:: void evolve(quantumdata::StateVector<RANK>&, const structure::QuantumSystem<RANK>&, const ParsEvolution&, V)

  template<int RANK, typename V>

--------------------
MCWF trajectory
--------------------

.. class:: quantumtrajectory::MCWF_Trajectory

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
