// -*- C++ -*-
#ifndef QUANTUMTRAJECTORY_MCWF_TRAJECTORY_H_INCLUDED
#define QUANTUMTRAJECTORY_MCWF_TRAJECTORY_H_INCLUDED

#include "MCWF_TrajectoryFwd.h"

#include "StateVectorFwd.h"

#include "MCWF_TrajectoryLogger.h"
#include "Structure.h"

#include "StochasticTrajectory.h"

#include <boost/tuple/tuple.hpp>


namespace quantumtrajectory {


///////////////////////////////////////
//
// Monte Carlo wave-function trajectory
//
///////////////////////////////////////

// NEEDS_WORK factor out template-parameter independent code


#define BASE_class trajectory::Stochastic<typename quantumdata::Types<RANK>::StateVectorLow, const quantumdata::StateVector<RANK>&>


/// Implements a single Monte Carlo wave-function trajectory
/**
 * In the framework, a single \ref mcwftrajectory "Monte Carlo wave-function step" at time \f$t\f$ (at which point, if the system inherits from structure::Exact,
 * the Schrödinger- and interaction pictures coincide) is implemented as a sequence of the following stages:
 * -# Coherent time development is applied:
 *   -# If the system time evolution has Hamiltonian part, it is evolved with an adaptive-size step (cf. evolved::Evolved). This takes the system into \f$t+\Delta t\f$.
 *   -# The exact part (if any) of the time evolution is applied, making that the Schrödinger and interaction pictures coincide again at \f$t+\Delta t\f$.
 *   -# The state vector is renormalized.
 * -# If the system is *not* Liouvillean, the timestep ends here, reducing to a simple ODE evolution. Otherwise:
 *   -# The rates (probabilities per unit time) corresponding to all jump operators are calculated. If some rates are found negative
 *      (“special jump”, cf. explanation at structure::Liouvillean::probabilities, then \f$J_\text{at}\ket\Psi\f$ is calculated (and tabulated) instead,
 *      and the rate is calculated as \f$\delta r_\text{at}=\norm{J_\text{at}\ket\Psi}^2\f$. \see \ref specialjump
 *   -# First, it is verified whether the total jump probability is not too big. This is performed on two levels:
 *     -# The total jump rate \f$\delta r\f$ is calculated.
 *     -# If \f$\delta r\delta t>\delta p_\text{limit}'\f$, the step is retraced: both the state vector and the state of the ODE stepper are restored to cached values
 *        at the beginning of the timestep, and phase I. is performed anew with a smaller stepsize \f$\delta p_\text{limit}/\delta r\f$. 
 *        With this, we ensure that \f$\delta p_\text{limit}\f$ is likely not to be overshot in the next try.
 *        \note It is assumed that \f$\delta p_\text{limit}'>\delta p_\text{limit}\f$, their ratio being a parameter of the MCWF stepper.
 *     -# If just \f$\delta r\delta t_\text{next}>\delta p_\text{limit}\f$ (where \f$\Delta t_\text{next}\f$ is a guess for the next timestep given by the ODE stepper),
 *        the coherent step is accepted, but the timestep to try next is modified, to reduce the likeliness of overshoot: \f$\delta t_\text{next}\longrightarrow\delta p_\text{limit}/\delta r\f$.
 *     \see The discussion at Sec. \ref anadaptivemcwfmethod "An adaptive MCWF method".
 *   -# After a successful coherent step resulting in an acceptable \f$\delta r\f$, the possible occurence of a quantum jump is considered:
 *      It is randomly decided which (if any) of the jumps to perform. If it is found to be a special jump, then the tabulated \f$J_\text{at}\ket\Psi\f$ is taken.
 * 
 * \note In phase 2.b.ii, another approach would be not to trace back the whole step, but make a coherent step *backwards* to an intermediate time instant found by linear interpolation.
 * This has several drawbacks, however, the most significant being that in the ODE stepper, it is not clear what to take as the timestep to try at the point when the direction of time is reversed.
 * (Although in evolved::Evolved it is simply taken to be the timestep done in the last step…)
 */

template<int RANK>
class MCWF_Trajectory : public BASE_class
{
public:
  typedef structure::Exact        <RANK> Exact      ;
  typedef structure::Hamiltonian  <RANK> Hamiltonian;
  typedef structure::Liouvillean  <RANK> Liouvillean;
  typedef structure::Averaged     <RANK> Averaged   ;

  typedef quantumdata::StateVector<RANK> StateVector;

  typedef typename StateVector::StateVectorLow StateVectorLow;

  typedef structure::QuantumSystemWrapper<RANK,true> QuantumSystemWrapper;

  typedef BASE_class Base;

#undef  BASE_class

  typedef boost::tuple<int,StateVectorLow> IndexSVL_tuple;

  using Base::getEvolved; using Base::getRandomized; using Base::getDtDid; using Base::getDtTry; using Base::getTime;

  template<typename SYS>
  MCWF_Trajectory(StateVector& psi, const SYS& sys, const ParsMCWF&, const StateVectorLow& =StateVectorLow());

  void derivs(double, const StateVectorLow&, StateVectorLow&) const;

  const StateVector& getPsi() const {return psi_;} 

  const MCWF_Logger& getLogger() const {return logger_;}
  
protected:
  std::ostream&    display_v(std::ostream&, int    ) const;
  std::ostream& displayKey_v(std::ostream&, size_t&) const;
  
  const QuantumSystemWrapper getQS() const {return qs_;}

  cpputils::iarchive&  readState_v(cpputils::iarchive& iar)       {Base:: readState_v(iar) & logger_; if (qs_.getEx()) tIntPic0_=getTime(); return iar;}
  cpputils::oarchive& writeState_v(cpputils::oarchive& oar) const {return Base::writeState_v(oar) & logger_;}

  std::ostream& logOnEnd_v(std::ostream& os) const {return logger_.onEnd(os);}
  
private:
  typedef std::vector<IndexSVL_tuple> IndexSVL_tuples;
  typedef typename Liouvillean::Rates Rates;

  void step_v(double); // performs one single adaptive-stepsize MCWF step of specified maximal length

  std::ostream& displayParameters_v(std::ostream&) const;

  const StateVector& toBeAveraged_v() const {return psi_;} 

  double                coherentTimeDevelopment(                    double Dt);
  const IndexSVL_tuples calculateSpecialRates  (      Rates* rates, double  t) const;

  bool                  manageTimeStep             (const Rates& rates, evolved::TimeStepBookkeeper*, bool logControl=true);

  void                  performJump                (const Rates&, const IndexSVL_tuples&, double); // LOGICALLY non-const
  // helpers to step---we are deliberately avoiding the normal technique of defining such helpers, because in that case the whole MCWF_Trajectory has to be passed

  mutable double tIntPic0_ ; // The time instant of the beginning of the current time step.

  StateVector& psi_;

  const QuantumSystemWrapper qs_;

  const double dpLimit_, overshootTolerance_;

  mutable MCWF_Logger logger_;

};


} // quantumtrajectory


/** \page mcwftrajectory Description of the MCWF method
 * 
 * The MCWF method \cite carmichael87 \cite dalibard92 \cite dum92 \cite molmer93 aims at the simulation of open quantum systems based on a stochastic (“Monte Carlo”) trajectory.
 * In terms of dimensionality, this is certainly a huge advantage as compared to solving the Master equation directly. On the other hand, stochasticity requires us to run many trajectories,
 * but the method provides an optimal sampling of the ensemble density operator so that the relative error is inversely proportional to the number of trajectories.
 * 
 * The optimal sampling is achieved by evolving the state vector in two steps, one deterministic and one stochastic (quantum jump).
 * Suppose that the Master equation of the system is of the form
 * \f[\dot\rho=\frac1{i\hbar}\comm{H}\rho+\Liou\rho\equiv\frac1{i\hbar}\comm{H}\rho+\sum_m\lp J_m\rho J_m^\dag-\frac12\comm{J_m^\dag J_m}{\rho}_+\rp\equiv\frac1{i\hbar}\lp\HnH\rho-\rho\HnH^\dagger\rp+\sum_mJ_m\rho J_m^\dag.\f]
 * 
 * the usual form in quantum optics, and, in fact, the most general (so-called Lindblad) form. The non-Hermitian Hamiltonian is defined as \f[\HnH=H-\frac{i\hbar}2\sum_m J^\dag_m J_m.\f]
 * At time \f$t\f$ the system is in a state with normalised state vector \f$\ket{\Psi(t)}\f$. To obtain the state vector at time \f$t+\delta t\f$ up to first order in \f$\delta t\f$:
 * -# The state vector is evolved according to the nonunitary dynamics \f[i\hbar\frac{d\ket{\Psi}}{dt}=\HnH\ket{\Psi}\f] to obtain (up to first order in \f$\delta t\f$) 
 * \f[\ket{\Psi_{\text{nH}}(t+\delta t)}=\lp1-\frac{i\HnH\,\delta t}\hbar\rp \ket{\Psi(t)}.\f]
 * Since \f$\HnH\f$ is non-Hermitian, this new state vector is not normalised. The square of its norm reads
 * \f[\braket{\Psi_{\text{nH}}(t+\delta t)}{\Psi_{\text{nH}}(t+\delta t)}=\bra{\Psi(t)}\lp1+\frac{iH^\dag_{\text{nH}}\,\delta t}\hbar\rp\lp1-\frac{i\HnH\,\delta t}\hbar\rp\ket{\Psi(t)}\equiv 1-\delta p,\f]
 * where \f$\delta p\f$ reads \f[\delta p=\delta t\,\frac i\hbar \bra{\Psi(t)}\HnH-H^\dag_{\text{nH}}\ket{\Psi(t)}\equiv\sum_m\delta p_m,\quad\delta p_m=\delta t\,\bra{\Psi(t)} J^\dag_m J_m\ket{\Psi(t)}\geq 0.\f]
 * Note that the timestep \f$\delta t\f$ should be small enough so that this first-order calculation be valid. In particular, we require that \f[\delta p\ll1.\f]
 * -# A possible quantum jump with total probability \f$\delta p\f$. For the physical interpretation of such a jump, cf. \cite dum92 \cite molmer93.
 * We choose a random number \f$\epsilon\f$ between 0 and 1, and if \f$\delta p<\epsilon\f$, which should mostly be the case, no jump occurs and for the new normalised state vector at
 * \f$t+\delta t\f$ we take \f[\ket{\Psi(t+\delta t)}=\frac{\ket{\Psi_{\text{nH}}(t+\delta t)}}{\sqrt{1-\delta p}}.\f]
 * If \f$\epsilon<\delta p\f$, on the other hand, a quantum jump occurs, and the new normalised state vector is chosen from among the different state vectors
 * \f$J_m\ket{\Psi(t)}\f$ according to the probability distribution \f$\Pi_m=\delta p_m/\delta p\f$: \f[\ket{\Psi(t+\delta t)}=\sqrt{\delta t}\frac{J_m\ket{\Psi(t)}}{\sqrt{\delta p_m}}.\f]
 * 
 * ## Refinement of the method
 * 
 * \anchor anadaptivemcwfmethod
 * 
 * ### An adaptive MCWF method
 * 
 * The method as described above has several shortcomings. Firstly, subsequent steps of no-jump evolution (Step 1) reduce to the first order (Euler) method of evolving the Schrödinger equation,
 * which is inappropriate in most cases of interest. Instead, to perform Step 1, we use an adaptive step-size ODE routine, usually the embedded Runge-Kutta Cash-Karp algorithm \cite numrec.
 * This has an intrinsic time-step management, which strives to achieve an optimal stepsize within specified (absolute and relative) error bounds.
 * 
 * A second problem with the original proposal is that it leaves unspecified what is to be done if *after the coherent step,* when calculating \f$\delta p\f$,
 * we find that condition \f$\delta p\ll1\f$ is not fulfilled. (With the fixed stepsize Euler method this happens if the jump *rate* grows too big, but with the adaptive stepsize algorithm,
 * it can happen also if the timestep grows too big.)
 * 
 * In the framework, we adopt a further heuristic in this case, introducing a tolerance interval for too big \f$\delta p\f$ values:
 * If at \f$t+\delta t\f$, \f$\delta p\f$ overshoots a certain \f$\delta p_\text{limit}\ll1\f$, then from the jump rate at this time instant,
 * a decreased stepsize is extracted to be feeded into the ODE stepper as the stepsize to try for the next timestep.
 * On the other hand, if \f$\delta p\f$ overshoots a \f$\delta p_\text{limit}'>\delta p_\text{limit}\f$, then the step is altogether discarded,
 * the state vector and the state of the ODE stepper being restored to cached states at time \f$t\f$.
 * 
 * \note In spite of the RKCK method being of order \f$O\lp\delta t^4\rp\f$, the whole method remains \f$O\lp\sqrt{\delta t}\rp\f$, since the treatment of jumps is essentially the same
 * as in the original proposal. (Events of multiple jumps in one timestep are neglected.)
 * 
 * \see quantumtrajectory::MCWF_Trajectory for the actual implementation
 * 
 * ### Exploiting interaction picture
 * 
 * In many situations it is worth using some sort of interaction picture, which means that instead of the non-Hermitian Schrödinger equation above, we strive to solve
 * \f[i\hbar\frac{d\ket{\Psi_{\text{I}}}}{dt}=U^{-1}\lp\HnH U-i\hbar\frac{dU}{dt}\rp\ket{\Psi_{\text I}},\f]
 * where \f$\ket{\Psi_{\text I}}=U^{-1}\ket\Psi\f$. Note that \f$U\f$ can be nonunitary, and therefore in general \f$U^{-1}\neq U^\dagger\f$.
 * 
 * \see On non-unitary transformations in quantum mechanics [these notes](http://optics.szfki.kfki.hu/~vukics/Pictures.pdf).
 * 
 * The two pictures are accorded after each timestep, i.e. before the timestep \f$\ket{\Psi_{\text I}(t)}=\ket{\Psi(t)}\f$ and after the timestep,
 * the transformation \f$\ket{\Psi(t+\delta t)}=U(\delta t)\ket{\Psi_{\text I}(t+\delta t)}\f$ is performed.
 * This we do on one hand for convenience and for compatibility with the case when no interaction picture is used,
 * but on the other hand also because \f$U(t)\f$ is nonunitary and hence for \f$t\to\infty\f$ some of its elements will become very large,
 * while others very small, possibly resulting in numerical problems. It is in fact advisable to avoid evaluating \f$U(t)\f$ with very large \f$t\f$ arguments.
 * 
 * \see the discussion at structure::Exact
 */


#endif // QUANTUMTRAJECTORY_MCWF_TRAJECTORY_H_INCLUDED
