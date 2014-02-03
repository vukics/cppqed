Description of the Master-equation evolution {#masterequation}
============================================

\tableofcontents

The \link quantumtrajectory::Master Master equation\endlink for the evolution of a density operator reads
\f[\dot\rho=\frac1{i\hbar}\comm{H}{\rho}+\sum_m\lp{J_m\rho{J_m^\dag}-\frac12\comm{J_m^\dag J_m}{\rho}_+}\rp\equiv\frac1{i\hbar}\lp\HnH\,\rho-\rho\HnH^\dag\rp+\sum_mJ_m\rho{J_m^\dag}=2\Re\lbr\frac\HnH{i\hbar}\rho\rbr+\sum_mJ_m\lp{J_m\rho}\rp^\dag,\f]
where, after the second equality, the non-Hermitian “Hamiltonian” is defined as (cf. \ref mcwftrajectory)
\f[\HnH=H-\frac{i\hbar}2\sum_m J^\dag_m J_m.\f]
Since both the Hamiltonian and Lindblad operators are expressed in a finite discrete basis, this is simply an \link evolved::Evolved ODE evolution\endlink.

Since each \link structure::Hamiltonian Hamiltonian\endlink and/or \link structure::Liouvillean lossy\endlink physical system in the framework  provides services
to calculate the effect of \f$\HnH/(i\hbar)\f$ (cf. structure::Hamiltonian::addContribution) and the jumps \f$J_m\f$ (cf. structure::Liouvillean::actWithJ) on state vectors, 
these being indispensable for the Monte Carlo wave-function method, no separate services are needed to calculate Master-equation evolution for these systems.

The steps are as follows:
  1. \link blitzplusplus::vfmsi::Iterator Iterate\endlink through the columns of the (in general, multi-)matrix \f$\rho\f$ to get \f$\HnH\,\rho/(i\hbar)\f$.
  2. Calculate two times the real part of this matrix *in place*.
  3. For each \f$m\f$:
    1. Take *a temporary copy* of \f$\rho\f$ and iterate through its columns to get \f$J_m\rho\f$.
    2. Take the Hermitian conjugate of this matrix *in place,* and repeat step (a).
    3. Add the contribution obtained in steps (a)-(b)-(a) to the matrix obtained in step 2.

The drawback of this method is that it necessitates a temporary copy of the full density operator. 
Alternatively, systems may provide specialized code to calculate \f$J_m\rho J_m^\dag\f$ directly, hence circumventing steps (a)-(b)-(a).

\todo Implement this possibility (cf. also [this tracker](http://sourceforge.net/p/cppqed/feature-requests/4/))

Limitations {#masterequationlimitations}
===========

Consider [these notes](http://optics.szfki.kfki.hu/~vukics/Pictures.pdf), whose notation we adopt. It is easy to see that a density-operator evolution of the form
\f[\dot\rho_T=2\Re\lbr\mathcal{A}^T\rho_T\rbr+\sum_mJ_m^T\rho_T\lp J_m^T\rp^\dag\f]
can always be unravelled up to first order in \f$\delta t\f$ via an \ref mcwftrajectory "MCWF-like" two-step stochastic process as:
1. \f$\ket{\Psi_T(t+\delta t)}=\frac{1+\delta t\mathcal{A}^T}{\sqrt{1-\delta p}}\ket{\Psi_T(t)}\f$, with probability \f$1-\delta p\f$.
2. \f$\ket{\Psi_T(t+\delta t)}=\frac{\sqrt{\delta t}J_m^T}{\sqrt{\delta p_m}}\ket{\Psi_T(t)}\f$, with probability \f$\delta p_m\f$ for either \f$m\f$, such that \f$\delta p=\sum_m\delta p_m\f$.

If we require that the second possibility do not alter the norm, as in the case of the MCWF method (the first may in the case when \f$\mathcal{A}^T\f$ is non-Hermitian),
that is, if the norm may change only *in time* (jumps considered *instantaneous*), then we get the further requirement that
\f$\delta p_m=\bra{\Psi_T(t)}\lp J_m^T\rp^\dag J_m^T\ket{\Psi_T(t)}\f$.

This shows that the \link quantumtrajectory::Master Master driver\endlink will only work correctly if \f[J_m^T\rho\lp J_m^T\rp^\dag=J_m\rho J_m^\dag\quad\forall\rho,\;m,\f]
because structure::Liouvillean::actWithJ calculates the latter, while quantumtrajectory::Master expects the former.

\see structure::ExactCommon::applicableInMaster