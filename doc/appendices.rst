.. _appendices:

********************************
Appendices
********************************


.. _MCWF_method:

==============================
Description of the MCWF method
==============================

The MCWF method \cite{carmichael87,dalibard92,dum92,molmer93} aims at the simulation of open quantum systems based on a stochastic ("Monte Carlo") trajectory. In terms of dimensionality, this is certainly a huge advantage as compared to solving the Master equation directly. On the other hand, stochasticity requires us to run many trajectories, but the method provides an optimal sampling of the ensemble density operator so that the relative error is inversely proportional to the number of trajectories.

The optimal sampling is achieved by evolving the state vector in two steps, one deterministic and one stochastic (quantum jump). Suppose that the Master equation of the system is of the form

.. math::
  :label: mm

  \dot\rho=\frac1{i\hbar}\comm{H}\rho+\Liou\rho\equiv\frac1{i\hbar}\comm{H}\rho+\sum_m\lp J_m\rho J_m^\dag-\frac12\comm{J_m^\dag J_m}{\rho}_+\rp\equiv\frac1{i\hbar}\comm{\HnH}\rho+\sum_mJ_m\rho J_m^\dag,

the usual form in quantum optics, and, in fact, the most general (so-called Lindblad) form. The non-Hermitian Hamiltonian is defined as
 
.. math::
  :label: HamnH

  \HnH=H-\frac{i\hbar}2\sum_m J^\dag_m J_m.

At time :math:`t` the system is in a state with normalised state vector :math:`\ket{\Psi(t)}`. To obtain the state vector at time :math:`t+\delta t` up to first order in :math:`\delta t`:

#. The state vector is evolved according to the nonunitary dynamics

   .. math::
     :label: DynnU

     i\hbar\frac{d\ket{\Psi}}{dt}=\HnH\ket{\Psi}

   to obtain (up to first order in :math:`\delta t`) 

   .. math::
     :label: PsinH
   
      \ket{\Psi_{\text{nH}}(t+\delta t)}=\lp1-\frac{i\HnH\,\delta t}\hbar\rp \ket{\Psi(t)}.

   Since :math:`\HnH` is non-Hermitian, this new state vector is not normalised. The square of its norm reads 

   .. math::

     \braket{\Psi_{\text{nH}}(t+\delta t)}{\Psi_{\text{nH}}(t+\delta t)}=\bra{\Psi(t)}\lp1+\frac{iH^\dag_{\text{nH}}\,\delta t}\hbar\rp\lp1-\frac{i\HnH\,\delta t}\hbar\rp\ket{\Psi(t)}\equiv 1-\delta p,

   where :math:`\delta p` reads 

   .. math::
   
     \delta p&=\delta t\,\frac i\hbar \bra{\Psi(t)}\HnH-H^\dag_{\text{nH}}\ket{\Psi(t)}\equiv\sum_m\delta p_m,

   .. math::

     \delta p_m&=\delta t\,\bra{\Psi(t)} J^\dag_m J_m\ket{\Psi(t)}\geq 0.

   Note that the timestep :math:`\delta t` should be small enough so that this first-order calculation be valid. In particular, we require that

   .. math::
     :label: dplimit

     \delta p\ll1.

#. A possible quantum jump with total probability :math:`\delta p`. For the physical interpretation of such a jump see e.g. Refs.\ \citep{dum92,molmer93}. We choose a random number :math:`\epsilon` between 0 and 1, and if :math:`\delta p<\epsilon`, which should mostly be the case, no jump occurs and for the new normalised state vector at :math:`t+\delta t` we take

  .. math::
    :label: renorm_ordodt

    \ket{\Psi(t+\delta t)}=\frac{\ket{\Psi_{\text{nH}}(t+\delta t)}}{\sqrt{1-\delta p}}.

  If :math:`\epsilon<\delta p`, on the other hand, a quantum jump occurs, and the new normalised state vector is chosen from among the different state vectors :math:`J_m\ket{\Psi(t)}` according to the probability distribution :math:`\Pi_m=\delta p_m/\delta p`:

  .. math::
  
    \ket{\Psi(t+\delta t)}=\sqrt{\delta t}\frac{J_m\ket{\Psi(t)}}{\sqrt{\delta p_m}}.


------------------------
Refinement of the method
------------------------

.. _MCWF_method_adaptive:

An adaptive MCWF method
^^^^^^^^^^^^^^^^^^^^^^^

The method as described above has several shortcomings. Firstly, subsequent steps of no-jump evolution (Step 1) reduce to the first order (Euler) method of evolving the Schrödinger equation, which is inappropriate in most cases of interest. Instead, to perform Step 1, we use an adaptive step-size ODE routine, usually the embedded Runge-Kutta Cash-Karp algorithm \cite{numrec}. This has an intrinsic time-step management, which strives to achieve an optimal stepsize within specified (absolute and relative) error bounds.

A second problem with the original proposal is that it leaves unspecified what is to be done if *after the coherent step,* when calculating :math:`\delta p`, we find that condition :eq:`dplimit` is not fulfilled. (With the fixed stepsize Euler method this happens if the jump *rate* grows too big, but with the adaptive stepsize algorithm, it can happen also if the timestep grows too big.)

In the framework, we adopt a further heuristic in this case, introducing a tolerance interval for too big :math:`\delta p` values: If at :math:`t+\delta t`, :math:`\delta p` overshoots a certain :math:`\delta p_\text{limit}\ll1`, then from the jump rate at this time instant, a decreased stepsize is extracted to be feeded into the ODE stepper as the stepsize to try for the next timestep. On the other hand, if :math:`\delta p` overshoots a :math:`\delta p_\text{limit}'>\delta p_\text{limit}`, then the step is altogether discarded, the state vector and the state of the ODE stepper being restored to cached states at time :math:`t`.


.. note::

  In spite of the RKCK method being of order :math:`O\lp\delta t^4\rp`, the whole method remains :math:`O\lp\sqrt{\delta t}\rp`, since the treatment of jumps is essentially the same as in the original proposal. (Events of multiple jumps in one timestep are neglected.)

For the implementation cf. :ref:`MCWF_Trajectory`.


Exploiting interaction picture
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In many situations it is worth using some sort of interaction picture, which means that instead of Eq. :eq:`DynnU` we strive to solve 

.. math::

  i\hbar\frac{d\ket{\Psi_{\text{I}}}}{dt}=U^{-1}\lp\HnH U-i\hbar\frac{dU}{dt}\rp\ket{\Psi_{\text I}},

where :math:`\ket{\Psi_{\text I}}=U^{-1}\ket\Psi`. Note that :math:`U` can be nonunitary, and therefore in general :math:`U^{-1}\neq U^\dagger`. (On non-unitary transformations in quantum mechanics cf. `these notes <http://optics.szfki.kfki.hu/~vukics/Pictures.pdf>`_.) The two pictures are accorded after each timestep, i.e. before the timestep :math:`\ket{\Psi_{\text I}(t)}=\ket{\Psi(t)}` and after the timestep the transformation :math:`\ket{\Psi(t+\delta t)}=U(\delta t)\ket{\Psi_{\text I}(t+\delta t)}` is performed. This we do on one hand for convenience and for compatibility with the case when no interaction picture is used, but on the other hand also because :math:`U(t)` is nonunitary and hence for :math:`t\to\infty` some of its elements will become very large, while others very small, possibly resulting in numerical problems. It is in fact advisable to avoid evaluating :math:`U(t)` with very large :math:`t` arguments.

Cf. also the discussion at :ref:`Exact`.


.. _NonOrthogonalFormalism:

========================
Nonorthogonal basis sets
========================

---------
Formalism
---------

We adopt the so-called covariant-contravariant formalism (see also \cite{artacho91}), which in physics is primarily known from the theory of relativity. Assume we have a basis :math:`\lbr\ket{i}\rbr_{i\in\mathbb{N}}`, where the basis vectors are nonorthogonal, so that the metric tensor

.. math::

  g_{ij}\equiv\braket{i}{j}

is nondiagonal. The contravariant components of a state vector :math:`\ket\Psi` are then defined as the expansion coefficients

.. math::

  \ket\Psi\equiv\Psi^i\ket{i},

where we have adopted the convention that there is summation for indeces appearing twice. In this case one index has to be down, while the other one up, otherwise we have as it were a syntax error. This ensures that the result is independent of the choice of basis, as one would very well expect e.g. from the trace of an operator. The covariant components are the projections

.. math::

  \Psi_i\equiv\braket{i}\Psi.

The two are connected with the metric tensor: 

.. math::

  \Psi_i=g_{ij}\Psi^j,\quad\Psi^i=g^{ij}\Psi_j,

where it is easy to show that :math:`g^{ij}` is the matrix inverse of the above defined :math:`g_{ij}`:

.. math::

  g_{ik}g^{kj}\equiv g_i^j=\delta_{ij}.

For the representation of operators assume :math:`A` is an operator, then

.. math::

  \ket\Phi\equiv A\ket\Psi=\Psi^jA\ket{j},

multiplying from the left by :math:`\bra{i}` yields

.. math::

  \Phi_i=\bra iA\ket j\Psi^j\equiv A_{ij}\Psi^j.

While the definition of the up-indeces matrix of the operator can be read from

.. math::

  \Phi=\Phi^i\ket i\equiv A^{ij}\Psi_j\ket{i}=A^{ij}\braket{j}\Psi\ket{i}=\lp A^{ij}\ket{i}\bra{j}\rp\ket\Psi=A\ket\Psi.


Hermitian conjugation
^^^^^^^^^^^^^^^^^^^^^

From the above it is easy to see that

.. math::

  \lp A^\dag\rp_{ij}=\lp A_{ji}\rp^*,\;\text{and}\;\lp A^\dag\rp^{ij}=\lp A^{ji}\rp^*.

The situation, however, is more involved with the mixed-index matrices since conjugation mixes the two kinds:

.. math::

  {\lp A^\dag\rp_i}^j=\lp {A^j}_i\rp^*,\;\text{and}\;{\lp A^\dag\rp^i}_j=\lp {A_j}^i\rp^*.

Here we prove the second:

.. math::

  {A_j}^i=A_{jk}g^{ki}=\lp\lp A^\dag\rp_{kj}\rp^*g^{ki},

so that

.. math::

  \lp{A_j}^i\rp^*=\lp A^\dag\rp_{kj}\lp g^{ki}\rp^*=g^{ik}\lp A^\dag\rp_{kj}={\lp A^\dag\rp^i}_j,

where we have used the Hermiticity of the metric tensor.

.. _NonOrthogonalFormalismImplications:

------------
Implications
------------

If one is to apply the above formalism for density matrices in nonorthogonal bases, one has to face some strange implications: The well-known properties of Hermiticity and unit trace are split: :math:`\rho_{ij}` and :math:`\rho^{ij}` will be Hermitian matrices, but their trace is irrelevant, in fact, it is not conserved during time evolution, while :math:`{\rho_i}^j` and :math:`{\rho^i}_j` have unit conserved trace.

This means that e.g. the :class:`quantumtrajectory::Master` driver has to use either :math:`\rho_{ij}` or :math:`\rho^{ij}` (by convention, it uses the latter, which also means that :class:`quantumtrajectory::MCWF_Trajectory` uses :math:`\Psi^i`) to be able to calculate the evolution using the real part as in \Eqref{eq:Master-MCWF} (otherwise, it is not the real part which appears there). On the other hand, for all normalizations one index has to be pulled, and this (at least in the case of density operators) has to be done in place to avoid excessive copying, which also means that the index has to be pulled back again afterwards.

The situation is even more delicate if one wants to calculate averages or jump probabilities for a subsystem. Then, as explained in \Sref{sec:Jump} and \Sref{sec:Displayed}, the partial density operator has to be used. Imagine that we want to calculate the expectation value of an operator :math:`O` acting on a given subsystem. Imagine that we are using an MCWF trajectory, so that we dispose of the expansion coefficients :math:`\Psi^{iI}`, where :math:`i` is the index (quantum number) for the given subsystem, while :math:`I` is a multi-index, fusing the (dummy) indeces of the other subsystems, which are all up. Now the matrix of partial density operator of the given subsystem reads

\begin{equation}
\label{eq:NO_PartialTrace}
\rho^{ij}=\lp{\Psi^i}_I\rp^*\Psi^{jI},
\end{equation}

and then, since we dispose of :math:`O_{ji}` (this is most easily obtained, since only scalar products have to be calculated), the desired expectation value is obtained as

\begin{equation}
\avr O=\rho^{ij}O_{ji}.
\end{equation}

This means that somehow all the indeces *but* the index of the given subsystem have to be pulled. Then, when we move to the next subsystem (this is what \code{Composite}'s \code{Display} does cf \Sref{sec:Composite}) to calculate the expectation value of one of its operators, we have to pull up again its index, and pull down the index of the previous subsystem, and so on. This is tedious, and possibly also wasteful.

What happens in practice is that the \code{MCWF\_Trajectory} driver, when it comes to display, pulls \emph{all indeces} by calling the function \code{Pull} for the complete \code{NonOrthogonal}, and passes both the covariant and contravariant vectors to \code{Display}. The summation in \Eqref{eq:NO_PartialTrace} is then performed by \code{ElementDisplay}, while the \code{Average} of the given subsystem works with :math:`{\rho_i}^j` (check definition!!!!!!) to calculate :math:`\avr O={\rho_i}^j{O_j}^i`, so that the additional pull has to be implemented here on :math:`O`. The same procedure is repeated with \code{Jump}, since to obtain the jump probability, also an expectation value has to be calculated. The drivers \code{Master} and \code{EnsembleMCWF} working with density operators :math:`\rho^{IJ}` pull one (multi-)index in place to obtain :math:`{\rho_I}^J`, passes this for display, and afterwards pull the index back again.

A given subsystem's \code{Average} function can therefore receive four possible sets of parameters:

\begin{enumerate}
\item for orthogonal basis and state-vector evolution only a single
array of \code{double}s, representing :math:`\Psi_I=\Psi^I`.
\item for nonorthogonal basis and state-vector evolution two arrays,
representing :math:`\Psi_I` and :math:`\Psi^I` separately.
\item for orthogonal basis and density-operator evolution a single
array, representing the matrix of :math:`\rho`.
\item for nonorthogonal basis and density-operator evolution a single
array, representing :math:`{\rho_I}^J`.
\end{enumerate}

However, an actual implementor of such a function does not see the difference between these cases, since as explained in \Chref{chap:Indexing}, they are hidden behind the \code{Indexing} interface, an abstract class, which is then implemented to cover the four separate cases. In fact, if the given subsystem is represented in an orthogonal basis, it does not really have to know about the existence of nonorthogonal bases, the only thing which may come as a surprise is that the diagonal elements of :math:`{\rho_I}^J` are in general not real. This is why the diagonal function of \code{Indexing} returns complex...

For reference we finally list the Master equation for :math:`\rho^{ij}`:


================
Some conventions
================

The Schr\"odinger equation is:
\begin{equation}
\label{eq:SchEq}
i\hbar\frac{d\ket\Psi}{dt}=H\ket\Psi.
\end{equation}
It may sound strange to call this a convention, but the position
of :math:`i` \emph{is} conventional in this equation. Of course, it has to
be consistent with other postulates, in particular the canonical
commutation relations, e.g.
\begin{equation}
\comm xp=i\hbar.
\end{equation}
The consistency of these two equations can be seen e.g. from Ehrenfest's
theorem.

\newcommand\Hinter{H^{\text{int}}}
\newcommand\IAP[1]{{#1_{\text{I}}}}

From \Eqref{eq:SchEq} we can derive the necessary definitions for
interaction picture. If the Hamiltonian is of the form
\begin{equation}
H\equiv H^{(0)}+\Hinter,
\end{equation}
where :math:`H^{(0)}` is something we can solve, then the operator :math:`U`
transforming between Schr\"odinger picture and the interaction picture
defined by :math:`H^{(0)}`
\begin{equation}
\ket\Psi\equiv U(t)\ket{\IAP\Psi}
\end{equation}
reads
\begin{equation}
U(t)=e^{-\frac i\hbar H^{(0)} t}.
\end{equation}
The Schr\"odinger equation reads
\begin{equation}
\label{eq:SchEqIAP}
i\hbar\frac{d\ket{\IAP\Psi}}{dt}=U^{-1}(t)\Hinter U(t)\ket{\IAP\Psi}
\equiv\IAP\Hinter.
\end{equation}
If :math:`\mathcal O` is an operator then in interaction picture it
evolves as
\begin{equation}
\frac{d\IAP{\mathcal O}}{dt}=
\frac i\hbar\comm{H^{(0)}}{\IAP{\mathcal O}}+
U^{-1}(t)\frac{\partial\mathcal O}{\partial t}U(t)
\end{equation}

For p look at how Sakurai does
