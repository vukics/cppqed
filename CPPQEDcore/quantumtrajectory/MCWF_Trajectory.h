// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_MCWF_TRAJECTORY_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_MCWF_TRAJECTORY_H_INCLUDED

#include "StateVector.h"

#include "MCWF_TrajectoryLogger.h"
#include "QuantumTrajectory.h"
#include "Structure.h"

#include "StochasticTrajectory.h"

#include <tuple>


namespace quantumtrajectory {

/// Auxiliary tools to MCWF_Trajectory  
namespace mcwf {


/// Aggregate of parameters pertaining to \link MCWF_Trajectory MCWF\endlink simulations
/** \copydetails trajectory::ParsRun */
template<typename RandomEngine>
struct Pars : public trajectory::ParsStochastic<RandomEngine> {
  
  double
    &dpLimit, ///< the parameter \f$\delta p_\text{limit}\f$ (cf. 2.b.ii \link MCWF_Trajectory here\endlink)
    &overshootTolerance; ///< the parameter \f$\delta p_\text{limit}'/\delta p_\text{limit}\f$ (cf. 2.b.ii \link MCWF_Trajectory here\endlink)

  size_t
    &nBins, ///< governs how many bins should be used for the histogram of jumps created by ensemble::streamLog (a zero value means a heuristic automatic determination)
    &nJumpsPerBin; ///< the average number of jumps per bin in the histogram of jumps for the case of heuristic bin-number determination

  Pars(parameters::Table& p, const std::string& mod="")
    : trajectory::ParsStochastic<RandomEngine>{p,mod},
      dpLimit(p.addTitle("MCWF_Trajectory",mod).add("dpLimit",mod,"MCWFS stepper total jump probability limit",0.01)),
      overshootTolerance(p.add("overshootTolerance",mod,"Jump probability overshoot tolerance factor",10.)),
      nBins(p.add("nBins",mod,"number of bins used for the histogram of jumps created by EnsembleMCWF",size_t(0))),
      nJumpsPerBin(p.add("nJumpsPerBin",mod,"average number of jumps per bin in the histogram of jumps for the case of heuristic bin-number determination",size_t(50)))
    {}

};


} // mcwf


#define BASE_class trajectory::Stochastic<structure::AveragedCommon::Averages,\
                                          typename quantumdata::Types<RANK>::StateVectorLow,\
                                          quantumdata::StateVector<RANK>,\
                                          RandomEngine>


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
 *        With this, we ensure that \f$\delta p_\text{limit}\f$ (the parameter mcwf::Pars::dpLimit) is likely not to be overshot in the next try.
 *        \note It is assumed that \f$\delta p_\text{limit}'>\delta p_\text{limit}\f$, their ratio being a parameter (mcwf::Pars::overshootTolerance) of the MCWF stepper.
 *     -# If just \f$\delta r\delta t_\text{next}>\delta p_\text{limit}\f$ (where \f$\Delta t_\text{next}\f$ is a guess for the next timestep given by the ODE stepper),
 *        the coherent step is accepted, but the timestep to try next is modified, to reduce the likeliness of overshoot: \f$\delta t_\text{next}\longrightarrow\delta p_\text{limit}/\delta r\f$.
 *     \see The discussion at Sec. \ref anadaptivemcwfmethod "An adaptive MCWF method".
 *   -# After a successful coherent step resulting in an acceptable \f$\delta r\f$, the possible occurence of a quantum jump is considered:
 *      It is randomly decided which (if any) of the jumps to perform. If it is found to be a special jump, then the tabulated \f$J_\text{at}\ket\Psi\f$ is taken.
 * 
 * \tparamRANK
 * 
 * \note In phase 2.b.ii, another approach would be not to trace back the whole step, but make a coherent step *backwards* to an intermediate time instant found by linear interpolation.
 * This has several drawbacks, however, the most significant being that in the ODE stepper, it is not clear what to take as the timestep to try at the point when the direction of time is reversed.
 * (Although in evolved::Evolved it is simply taken to be the timestep done in the last step…)
 * 
 */

template<int RANK, typename RandomEngine>
class MCWF_Trajectory : public QuantumTrajectory<RANK,BASE_class>
{
public:
  MCWF_Trajectory(MCWF_Trajectory&&) = default; MCWF_Trajectory& operator=(MCWF_Trajectory&&) = default;

  typedef structure::Exact        <RANK> Exact      ;
  typedef structure::Hamiltonian  <RANK> Hamiltonian;
  typedef structure::Liouvillean  <RANK> Liouvillean;
  typedef structure::Averaged     <RANK> Averaged   ;

  typedef quantumdata::StateVector<RANK> StateVector;

  typedef typename StateVector::StateVectorLow StateVectorLow;

private:
  typedef quantumtrajectory::QuantumTrajectory<RANK,BASE_class> QuantumTrajectory;
  typedef BASE_class Base;

#undef  BASE_class

  typedef std::tuple<int,StateVectorLow> IndexSVL_tuple;

public:
  /// Templated constructor with the same idea as Master::Master
  MCWF_Trajectory(StateVector&& psi, ///< the state vector to be evolved
                  typename structure::QuantumSystem<RANK>::Ptr sys, ///< object representing the quantum system
                  const mcwf::Pars<RandomEngine>& p, ///< parameters of the evolution
                  const StateVectorLow& scaleAbs=StateVectorLow() ///< has the same role as `scaleAbs` in Master::Master
                 );

  /// \name Getters
  //@{
  const mcwf::Logger& getLogger() const {return logger_;}
  //@}
  
  using Base::getTime;

protected:
  using Base::getEvolved; using Base::getDtTry;

  typename Base::StreamReturnType stream_v(std::ostream& os, int precision) const override
  {
    return structure::QuantumSystemWrapper<RANK,true>::stream(getTime(),psi_,os,precision);
  } ///< Forwards to structure::Averaged::stream
  
  std::ostream& streamKey_v(std::ostream&, size_t&) const override; ///< Forwards to structure::Averaged::streamKey

  /// Forwards to QuantumTrajectory::readStateMore_v (that involves setting \link QuantumTrajectory::getT0 `t0`\endlink) + serializes mcwf::Logger state
  cpputils::iarchive&  readStateMore_v(cpputils::iarchive& iar) override {return QuantumTrajectory::readStateMore_v(iar) & logger_;}
  /// Forwards to Base::writeStateMore_v + serializes mcwf::Logger state
  cpputils::oarchive& writeStateMore_v(cpputils::oarchive& oar) const override {return Base::writeStateMore_v(oar) & logger_;}

  std::ostream& logMoreOnEnd_v(std::ostream& os) const override {return logger_.onEnd(os);} ///< calls mcwf::Logger::onEnd
  
private:
  typedef std::vector<IndexSVL_tuple> IndexSVL_tuples;
  typedef typename Liouvillean::Rates Rates;

  void step_v(double, std::ostream&) override; // performs one single adaptive-stepsize MCWF step of specified maximal length

  std::ostream& streamParameters_v(std::ostream&) const override;

  StateVector averaged_v() const override {return psi_;}

  const std::string trajectoryID_v() const override {return "MCWF_Trajectory";}

  double coherentTimeDevelopment(double Dt, std::ostream&);
  const IndexSVL_tuples calculateSpecialRates(Rates* rates, double t) const;

  bool manageTimeStep (const Rates& rates, evolved::TimeStepBookkeeper*, std::ostream&, bool logControl=true);

  void performJump (const Rates&, const IndexSVL_tuples&, double, std::ostream&); // LOGICALLY non-const
  // helpers to step---we are deliberately avoiding the normal technique of defining such helpers, because in that case the whole MCWF_Trajectory has to be passed

  StateVector psi_;

  const double dpLimit_, overshootTolerance_;

  mutable mcwf::Logger logger_;

};


} // quantumtrajectory


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_MCWF_TRAJECTORY_H_INCLUDED
