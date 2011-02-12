.. _quantumtrajectory:

===================================
The ``quantumtrajectory`` namespace
===================================

.. image:: figures/trajectory.png
   :width: 801
   
   
.. todo::

   Apply Boost.Parameter for facilitating call of functions with complex signature, e.g. constructors of complex classes, especially in cases like in the trajectory bundle where many sensible default parameter values could be defined. (Eg. it is extremely tedious scaleAbs needs to be specified each time.) This is in general most useful for constructors of complex classes.

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
