Description of the MCWF method {#mcwftrajectory}
==============================

\tableofcontents

The MCWF method \cite carmichael87 \cite dalibard92 \cite dum92 \cite molmer93 aims at the simulation of open quantum systems based on a stochastic (“Monte Carlo”) trajectory.
In terms of dimensionality, this is certainly a huge advantage as compared to solving the Master equation directly. On the other hand, stochasticity requires us to run many trajectories,
but the method provides an optimal sampling of the ensemble density operator so that the relative error is inversely proportional to the number of trajectories.

The optimal sampling is achieved by evolving the state vector in two steps, one deterministic and one stochastic (quantum jump).
Suppose that the Master equation of the system is of the form
\f[\dot\rho=\frac1{i\hbar}\comm{H}\rho+\Liou\rho\equiv\frac1{i\hbar}\comm{H}\rho+\sum_m\lp J_m\rho J_m^\dag-\frac12\comm{J_m^\dag J_m}{\rho}_+\rp\equiv\frac1{i\hbar}\lp\HnH\rho-\rho\HnH^\dagger\rp+\sum_mJ_m\rho J_m^\dag.\f]

This is the usual form in quantum optics, and, in fact, the most general (so-called Lindblad) form.
The non-Hermitian Hamiltonian is defined as \f[\HnH=H-\frac{i\hbar}2\sum_m J^\dag_m J_m.\f]
At time \f$t\f$ the system is in a state with normalised state vector \f$\ket{\Psi(t)}\f$. To obtain the state vector at time \f$t+\delta t\f$ up to first order in \f$\delta t\f$:
-# The state vector is evolved according to the nonunitary dynamics \f[i\hbar\frac{d\ket{\Psi}}{dt}=\HnH\ket{\Psi}\f] to obtain (up to first order in \f$\delta t\f$) 
\f[\ket{\Psi_{\text{nH}}(t+\delta t)}=\lp1-\frac{i\HnH\,\delta t}\hbar\rp \ket{\Psi(t)}.\f]
Since \f$\HnH\f$ is non-Hermitian, this new state vector is not normalised. The square of its norm reads
\f[\braket{\Psi_{\text{nH}}(t+\delta t)}{\Psi_{\text{nH}}(t+\delta t)}=\bra{\Psi(t)}\lp1+\frac{iH^\dag_{\text{nH}}\,\delta t}\hbar\rp\lp1-\frac{i\HnH\,\delta t}\hbar\rp\ket{\Psi(t)}\equiv 1-\delta p,\f]
where \f$\delta p\f$ reads \f[\delta p=\delta t\,\frac i\hbar \bra{\Psi(t)}\HnH-H^\dag_{\text{nH}}\ket{\Psi(t)}\equiv\sum_m\delta p_m,\quad\delta p_m=\delta t\,\bra{\Psi(t)} J^\dag_m J_m\ket{\Psi(t)}\geq 0.\f]
Note that the timestep \f$\delta t\f$ should be small enough so that this first-order calculation be valid. In particular, we require that \f[\delta p\ll1.\f]
-# A possible quantum jump with total probability \f$\delta p\f$. For the physical interpretation of such a jump, cf. \cite dum92 \cite molmer93.
We choose a random number \f$\epsilon\f$ between 0 and 1, and if \f$\delta p<\epsilon\f$, which should mostly be the case, no jump occurs and for the new normalised state vector at
\f$t+\delta t\f$ we take \f[\ket{\Psi(t+\delta t)}=\frac{\ket{\Psi_{\text{nH}}(t+\delta t)}}{\sqrt{1-\delta p}}.\f]
If \f$\epsilon<\delta p\f$, on the other hand, a quantum jump occurs, and the new normalised state vector is chosen from among the different state vectors
\f$J_m\ket{\Psi(t)}\f$ according to the probability distribution \f$\Pi_m=\delta p_m/\delta p\f$: \f[\ket{\Psi(t+\delta t)}=\sqrt{\delta t}\frac{J_m\ket{\Psi(t)}}{\sqrt{\delta p_m}}.\f]

Refinement of the method {#refinementofthemethod}
========================

An adaptive MCWF method {#anadaptivemcwfmethod}
-----------------------

The method as described above has several shortcomings. Firstly, subsequent steps of no-jump evolution (Step 1) reduce to the first order (Euler) method of evolving the Schrödinger equation,
which is inappropriate in most cases of interest. Instead, to perform Step 1, we use an adaptive step-size ODE routine, usually the embedded Runge-Kutta Cash-Karp algorithm \cite numrec.
This has an intrinsic time-step management, which strives to achieve an optimal stepsize within specified (absolute and relative) error bounds.

A second problem with the original proposal is that it leaves unspecified what is to be done if *after the coherent step,* when calculating \f$\delta p\f$,
we find that condition \f$\delta p\ll1\f$ is not fulfilled. (With the fixed stepsize Euler method this happens if the jump *rate* grows too big, but with the adaptive stepsize algorithm,
it can happen also if the timestep grows too big.)

In the framework, we adopt a further heuristic in this case, introducing a tolerance interval for too big \f$\delta p\f$ values:
If at \f$t+\delta t\f$, \f$\delta p\f$ overshoots a certain \f$\delta p_\text{limit}\ll1\f$, then from the jump rate at this time instant,
a decreased stepsize is extracted to be feeded into the ODE stepper as the stepsize to try for the next timestep.
On the other hand, if \f$\delta p\f$ overshoots a \f$\delta p_\text{limit}'>\delta p_\text{limit}\f$, then the step is altogether discarded,
the state vector and the state of the ODE stepper being restored to cached states at time \f$t\f$.

\note In spite of the RKCK method being of order \f$O\lp\delta t^4\rp\f$, the whole method remains \f$O\lp\sqrt{\delta t}\rp\f$, since the treatment of jumps is essentially the same
as in the original proposal. (Events of multiple jumps in one timestep are neglected.)

\see quantumtrajectory::MCWF_Trajectory for the actual implementation

Exploiting interaction picture {#exploitinginteractionpicture}
------------------------------

In many situations it is worth using some sort of interaction picture, which means that instead of the non-Hermitian Schrödinger equation above, we strive to solve
\f[i\hbar\frac{d\ket{\Psi_{\text{I}}}}{dt}=U^{-1}\lp\HnH U-i\hbar\frac{dU}{dt}\rp\ket{\Psi_{\text I}},\f]
where \f$\ket{\Psi_{\text I}}=U^{-1}\ket\Psi\f$. Note that \f$U\f$ can be nonunitary, and therefore in general \f$U^{-1}\neq U^\dagger\f$.

\see On non-unitary transformations in quantum mechanics [these notes](http://optics.szfki.kfki.hu/~vukics/Pictures.pdf).

The two pictures are accorded after each timestep, i.e. before the timestep \f$\ket{\Psi_{\text I}(t)}=\ket{\Psi(t)}\f$ and after the timestep,
the transformation \f$\ket{\Psi(t+\delta t)}=U(\delta t)\ket{\Psi_{\text I}(t+\delta t)}\f$ is performed.
This we do on one hand for convenience and for compatibility with the case when no interaction picture is used,
but on the other hand also because \f$U(t)\f$ is nonunitary and hence for \f$t\to\infty\f$ some of its elements will become very large,
while others very small, possibly resulting in numerical problems. It is in fact advisable to avoid evaluating \f$U(t)\f$ with very large \f$t\f$ arguments.

\see the discussion at structure::Exact



