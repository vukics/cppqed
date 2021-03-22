// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_MCWF_TRAJECTORY_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_MCWF_TRAJECTORY_H_INCLUDED

#include "StateVector.h"

#include "MCWF_TrajectoryLogger.h"
#include "QuantumTrajectory.h"
#include "Structure.h"

#include "StochasticTrajectory.h"

#include <boost/range/algorithm/find_if.hpp>

#include <tuple>


namespace quantumtrajectory {

/// Auxiliary tools to MCWF_Trajectory  
namespace mcwf {


/// Aggregate of parameters pertaining to \link MCWF_Trajectory MCWF\endlink simulations
/** \copydetails trajectory::ParsRun *//*
template<typename RandomEngine>
struct Pars : public cppqedutils::trajectory::ParsStochastic<RandomEngine> {
  
  double
    &dpLimit, ///< the parameter \f$\delta p_\text{limit}\f$ (cf. 2.b.ii \link MCWF_Trajectory here\endlink)
    &overshootTolerance; ///< the parameter \f$\delta p_\text{limit}'/\delta p_\text{limit}\f$ (cf. 2.b.ii \link MCWF_Trajectory here\endlink)

  size_t
    &nBins, ///< governs how many bins should be used for the histogram of jumps created by ensemble::streamLog (a zero value means a heuristic automatic determination)
    &nJumpsPerBin; ///< the average number of jumps per bin in the histogram of jumps for the case of heuristic bin-number determination

  Pars(parameters::Table& p, const std::string& mod="")
    : cppqedutils::trajectory::ParsStochastic<RandomEngine>{p,mod},
      dpLimit(p.addTitle("MCWF_Trajectory",mod).add("dpLimit",mod,"MCWFS stepper total jump probability limit",0.01)),
      overshootTolerance(p.add("overshootTolerance",mod,"Jump probability overshoot tolerance factor",10.)),
      nBins(p.add("nBins",mod,"number of bins used for the histogram of jumps created by EnsembleMCWF",size_t(0))),
      nJumpsPerBin(p.add("nJumpsPerBin",mod,"average number of jumps per bin in the histogram of jumps for the case of heuristic bin-number determination",size_t(50)))
    {}

};*/


} // mcwf



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

template<int RANK, typename ODE_Engine, typename RandomEngine>
class MCWF_Trajectory : private structure::QuantumSystemWrapper<RANK,true>
{
private:
  using QuantumSystemWrapper=structure::QuantumSystemWrapper<RANK,true>;
  
public:
  MCWF_Trajectory(MCWF_Trajectory&&) = default; MCWF_Trajectory& operator=(MCWF_Trajectory&&) = default;

  using StreamedArray=structure::AveragedCommon::Averages;

  typedef structure::Exact      <RANK> Exact      ;
  typedef structure::Hamiltonian<RANK> Hamiltonian;
  typedef structure::Liouvillean<RANK> Liouvillean;
  typedef structure::Averaged   <RANK> Averaged   ;

  typedef structure::QuantumSystem<RANK> QuantumSystem;

  typedef quantumdata::StateVector<RANK> StateVector;
  typedef typename StateVector::StateVectorLow StateVectorLow;

  using EnsembleAverageElement = StateVector;
  using EnsembleAverageResult = quantumdata::DensityOperator<RANK>;
  
  /// Templated constructor with the same idea as Master::Master
  template <typename SV>
  MCWF_Trajectory(typename QuantumSystem::Ptr sys, ///< object representing the quantum system
                  SV&& psi, ///< the state vector to be evolved
                  ODE_Engine ode, RandomEngine re, double dpLimit, double overshootTolerance, int logLevel)
  : QuantumSystemWrapper{sys,true}, psi_{std::forward<StateVector>(psi)},
    ode_{ode}, re_{re}, dpLimit_{dpLimit}, overshootTolerance_{overshootTolerance},
    logger_{logLevel,this->template nAvr<structure::LA_Li>()}
  {
    if (psi!=*sys) throw DimensionalityMismatchException("during QuantumTrajectory construction");
    if (t_) if(const auto li=this->getLi()) { // On startup, dpLimit should not be overshot, either.
      Rates rates(li->rates(0.,psi_)); calculateSpecialRates(&rates);
      double dtTry=ode_.getDtTry();
      manageTimeStep(rates,0,0,dtTry,std::clog,false);
      ode_.setDtTry(dtTry);
    }
  }
  
  auto getTime() const {return t_;}

  void step(double deltaT, std::ostream& logStream);
  
  auto getDtDid() const {return ode_.getDtDid();}
  
  std::ostream& streamParameters(std::ostream&) const;

  auto stream(std::ostream& os, int precision) const {return QuantumSystemWrapper::stream(getTime(),psi_,os,precision);} ///< Forwards to structure::Averaged::stream
  
  auto& readFromArrayOnlyArchive(cppqedutils::iarchive& iar) {StateVectorLow temp; iar & temp; psi_.getArray().reference(temp); return iar;}

  /** structure of MCWF_Trajectory archives:
  * metaData – array – time – ( odeStepper – odeLogger – dtDid – dtTry ) – randomEngine – logger
  */
  // state should precede time in order to be compatible with array-only archives
  auto& stateIO(cppqedutils::iarchive& iar)
  {
    StateVectorLow temp;
    ode_.stateIO(iar & temp & t_) & re_ & logger_;
    psi_.getArray().reference(temp);
    t0_=t_; // A very important step!
    return iar;
  }
  
  auto& stateIO(cppqedutils::oarchive& oar) {return ode_.stateIO(oar & psi_.getArray() & t_) & re_ & logger_;}

  std::ostream& streamKey(std::ostream& os) const {size_t i=3; return QuantumSystemWrapper::template streamKey<structure::LA_Av>(os,i);} ///< Forwards to structure::Averaged::streamKey

  std::ostream& logOnEnd(std::ostream& os) const {return logger_.onEnd(ode_.logOnEnd(os));} ///< calls mcwf::Logger::onEnd
  
private:
  typedef std::tuple<int,StateVectorLow> IndexSVL_tuple;

  typedef std::vector<IndexSVL_tuple> IndexSVL_tuples;
  typedef typename Liouvillean::Rates Rates;

  StateVector averaged() const {return psi_;}

  void coherentTimeDevelopment(double Dt, std::ostream&);
  
  const IndexSVL_tuples calculateSpecialRates(Rates* rates) const;

  bool manageTimeStep (const Rates& rates, double tCache, double dtDidCache, double & dtTryCache, std::ostream&, bool logControl=true);

  void performJump (const Rates&, const IndexSVL_tuples&, std::ostream&); // LOGICALLY non-const
  // helpers to step---we are deliberately avoiding the normal technique of defining such helpers, because in that case the whole MCWF_Trajectory has to be passed

  double t_=0., t0_=0.;
  
  StateVector psi_;

  ODE_Engine ode_;

  RandomEngine re_;
  
  const double dpLimit_, overshootTolerance_;

  mutable mcwf::Logger logger_;

};


} // quantumtrajectory


template <int RANK, typename ODE_Engine, typename RandomEngine>
struct cppqedutils::trajectory::MakeSerializationMetadata<quantumtrajectory::MCWF_Trajectory<RANK,ODE_Engine,RandomEngine>>
{
  static auto _() {return SerializationMetadata{typeid(dcomp).name(),"MCWF_Trajectory",RANK};}
};



template<int RANK, typename ODE_Engine, typename RandomEngine>
void quantumtrajectory::MCWF_Trajectory<RANK,ODE_Engine,RandomEngine>::coherentTimeDevelopment(double Dt, std::ostream& logStream)
{
  if (this->getHa()) {
    ode_.step(Dt,logStream,[this](const StateVectorLow& psi, StateVectorLow& dpsidt, double t) {
      dpsidt=0;
      this->getHa()->addContribution(t,psi,dpsidt,t0_);
    },t_,psi_.getArray());
  }
  else {
    double stepToDo=this->getLi() ? std::min(ode_.getDtTry(),Dt) : Dt; // Cf. tracker #3482771
    t_+=stepToDo;
  }
  
  // This defines three levels:
  // 1. System is Hamiltonian -> internal timestep control is used with dtTry possibly modified by Liouvillean needs (cf. manageTimeStep)
  // 2. System is not Hamiltonian, but it is Liouvillean -> dtTry is used, which is governed solely by Liouvillean in this case (cf. manageTimeStep)
  // 3. System is neither Hamiltonian nor Liouvillean (might be of Exact) -> a step of Dt is taken
  
  if (const auto ex=this->getEx()) {
    ex->actWithU(t_,psi_.getArray(),t0_);
    t0_=t_;
  }

  logger_.processNorm(psi_.renorm());

}


template<int RANK, typename ODE_Engine, typename RandomEngine>
auto quantumtrajectory::MCWF_Trajectory<RANK,ODE_Engine,RandomEngine>::calculateSpecialRates(Rates* rates) const -> const IndexSVL_tuples
{
  IndexSVL_tuples res;
  for (int i=0; i<rates->size(); i++)
    if ((*rates)(i)<0) {
      StateVector psiTemp(psi_);
      this->actWithJ(t_,psiTemp.getArray(),i);
      res.push_back(IndexSVL_tuple(i,psiTemp.getArray()));
      (*rates)(i)=cppqedutils::sqr(psiTemp.renorm());
    } // psiTemp disappears here, but its storage does not, because the ownership is taken over by the StateVectorLow in the tuple
  return res;
}


template<int RANK, typename ODE_Engine, typename RandomEngine>
bool quantumtrajectory::MCWF_Trajectory<RANK,ODE_Engine,RandomEngine>::
manageTimeStep(const Rates& rates, double tCache, double dtDidCache, double & dtTryCache, std::ostream& logStream, bool logControl)
{
  const double totalRate=boost::accumulate(rates,0.);
  const double dtDid=ode_.getDtDid(), dtTry=ode_.getDtTry();

  const double liouvilleanSuggestedDtTry=dpLimit_/totalRate;

  // Assumption: overshootTolerance_>=1 (equality is the limiting case of no tolerance)
  if (totalRate*dtDid>overshootTolerance_*dpLimit_) {
    dtTryCache=liouvilleanSuggestedDtTry;
    t_=tCache; /*ode_.setDtDid(dtDidCache);*/ ode_.setDtTry(dtTryCache);
    logger_.stepBack(logStream,totalRate*dtDid,dtDid,ode_.getDtTry(),t_,logControl);
    return true; // Step-back required.
  }
  else if ( totalRate*dtTry>dpLimit_) {
    logger_.overshot(logStream,totalRate*dtTry,dtTry,liouvilleanSuggestedDtTry,logControl);
  }

  // dtTry-adjustment for next step:
  ode_.setDtTry(this->getHa() ? std::min(ode_.getDtTry(),liouvilleanSuggestedDtTry) : liouvilleanSuggestedDtTry);

  return false; // Step-back not required.
}


template<int RANK, typename ODE_Engine, typename RandomEngine>
void quantumtrajectory::MCWF_Trajectory<RANK,ODE_Engine,RandomEngine>::performJump(const Rates& rates, const IndexSVL_tuples& specialRates, std::ostream& logStream)
{
  double random=std::uniform_real_distribution()(re_)/this->getDtDid();

  int lindbladNo=0; // TODO: this could be expressed with an iterator into rates
  for (; random>0 && lindbladNo!=rates.size(); random-=rates(lindbladNo++))
    ;

  if(random<0) { // Jump corresponding to Lindblad no. lindbladNo-1 occurs
    --lindbladNo; auto i=boost::find_if(specialRates,[=](const IndexSVL_tuple& j){return lindbladNo==std::get<0>(j);});
    if (i!=specialRates.end())
      // special jump
      psi_=std::get<1>(*i); // RHS already normalized above
    else {
      // normal  jump
      this->actWithJ(t_,psi_.getArray(),lindbladNo);
      double normFactor=sqrt(rates(lindbladNo));
      if (!boost::math::isfinite(normFactor)) throw std::runtime_error("Infinite detected in MCWF_Trajectory::performJump");
      psi_/=normFactor;
    }

    logger_.jumpOccured(logStream,t_,lindbladNo);
  }
}


template<int RANK, typename ODE_Engine, typename RandomEngine>
void quantumtrajectory::MCWF_Trajectory<RANK,ODE_Engine,RandomEngine>::step(double Dt, std::ostream& logStream)
{
  const StateVector psiCache(psi_); // deep copy
  double tCache=t_, dtDidCache=ode_.getDtDid(), dtTryCache=ode_.getDtTry();

  coherentTimeDevelopment(Dt,logStream);

  if (const auto li=this->getLi()) {

    Rates rates(li->rates(t_,psi_));
    IndexSVL_tuples specialRates=calculateSpecialRates(&rates);

    while (manageTimeStep(rates,tCache,dtDidCache,dtTryCache,logStream)) {
      psi_=psiCache;
      coherentTimeDevelopment(Dt,logStream); // the next try
      rates=li->rates(t_,psi_);
      specialRates=calculateSpecialRates(&rates);
    }

    // Jump
    performJump(rates,specialRates,logStream);

  }

  logger_.step();

}


template<int RANK, typename ODE_Engine, typename RandomEngine>
std::ostream& quantumtrajectory::MCWF_Trajectory<RANK,ODE_Engine,RandomEngine>::streamParameters(std::ostream& os) const
{
  using namespace std;
  
  this->streamCharacteristics(this->getQS()->streamParameters(
    ode_.streamParameters(os)<<"Random engine: "<<randomutils::EngineID_v<RandomEngine><<". Initial state: "<<re_
    <<"\nMCWF Trajectory Parameters: dpLimit="<<dpLimit_<<" (overshoot tolerance factor)="<<overshootTolerance_<<endl<<endl)
  )<<endl;

  if (const auto li=this->getLi()) {
    os<<"Decay channels:\n";
    {
      size_t i=0;
      li->streamKey(os,i);
    }
    os<<"Alternative Lindblads: ";
    {
      const Rates rates(li->rates(0,psi_));
      int n=0;
      for (int i=0; i<rates.size(); ++i) if (rates(i)<0) {os<<i<<' '; ++n;}
      if (!n) os<<"none";
    }
    os<<endl;
  }

  return os;
  
}



#endif // CPPQEDCORE_QUANTUMTRAJECTORY_MCWF_TRAJECTORY_H_INCLUDED
