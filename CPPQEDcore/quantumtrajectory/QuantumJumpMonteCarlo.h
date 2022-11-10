// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#pragma once

#include "StateVector.h"

#include "DensityOperatorStreamer.h"
#include "QuantumJumpMonteCarloLogger.h"
#include "QuantumTrajectory.h"
#include "Structure.h"

#include "StochasticTrajectory.h"

#include <tuple>


namespace quantumtrajectory {

/// Auxiliary tools to QuantumJumpMonteCarlo  
namespace qjmc {


/// Aggregate of parameters of QuantumJumpMonteCarlo
template<typename RandomEngine>
struct Pars : public cppqedutils::trajectory::ParsStochastic<RandomEngine> {
  
  double
    &dpLimit, ///< the parameter \f$\Delta p\f$
    &overshootTolerance; ///< the parameter \f$\Delta p'/\Delta p\f$

  size_t
    &nBins, ///< governs how many bins should be used for the histogram of jumps created by qjmc::EnsembleLogger::stream (a zero value means a heuristic automatic determination)
    &nJumpsPerBin; ///< the average number of jumps per bin in the histogram of jumps for the case of heuristic bin-number determination

  Pars(parameters::Table& p, const std::string& mod="")
    : cppqedutils::trajectory::ParsStochastic<RandomEngine>{p,mod},
      dpLimit(p.addTitle("QuantumJumpMonteCarlo",mod).add("dpLimit",mod,"QJMC stepper total jump probability limit",0.01)),
      overshootTolerance(p.add("overshootTolerance",mod,"Jump probability overshoot tolerance factor",10.)),
      nBins(p.add("nBins",mod,"number of bins used for a histogram of jumps",size_t(0))),
      nJumpsPerBin(p.add("nJumpsPerBin",mod,"average number of jumps per bin in the histogram of jumps for the case of heuristic bin-number determination",size_t(50)))
    {}

};


} // qjmc



/// Implements a single Quantum-Jump Monte Carlo trajectory
/** the stepwise adaptive algorithm is used, cf. Comp. Phys. Comm. 238:88 (2019) */
template<size_t RANK,
         ::structure::hamiltonian<RANK> HA, ::structure::liouvillian<RANK> LI, ::structure::expectation_values<RANK> EV,
         typename ODE_Engine, std::uniform_random_bit_generator RandomEngine>
class QuantumJumpMonteCarlo
{
public:
  QuantumJumpMonteCarlo(const QuantumJumpMonteCarlo&) = default; QuantumJumpMonteCarlo(QuantumJumpMonteCarlo&&) = default; QuantumJumpMonteCarlo& operator=(QuantumJumpMonteCarlo&&) = default;

  using StreamedArray=::structure::EV_Array;

  using StateVectorView=quantumdata::StateVectorView<RANK>;

  using EnsembleAverageElement = StateVectorView;
  using EnsembleAverageResult = quantumdata::DensityOperator<RANK>;
  
  QuantumJumpMonteCarlo(const HA& ha, const LI& li, const EV& ev,
                  StateVectorView<RANK> psi, ///< the state vector to be evolved
                  ODE_Engine ode, randomutils::EngineWithParameters<RandomEngine> re,
                  double dpLimit, double overshootTolerance, int logLevel)
  : ha_{ha}, li_{li}, ev_{ev},
    psi_{psi},
    ode_{ode}, re_{re}, dpLimit_{dpLimit}, overshootTolerance_{overshootTolerance},
    logger_{logLevel,::structure::nAvr<::structure::LA_Li>(sys_)}
  {
    if (!t_ && li.size() )  { // On startup, dpLimit should not be overshot, either.
      auto rates{calculateRates(li,0.,psi_)};
      manageTimeStep(rates,0,0,std::clog,false);
    }
  }
  
  auto getTime() const {return t_;}

  void step(double deltaT, std::ostream& logStream);
  
  auto getDtDid() const {return ode_.getDtDid();}
  auto getDtTry() const {return ode_.getDtTry();}
  
  std::ostream& streamParameters(std::ostream&) const;

  /// Forwards to ::structure::Averaged::stream
  auto stream(std::ostream& os, int precision) const {return ::structure::stream(sys_,getTime(),psi_,os,precision);}
  
  auto& readFromArrayOnlyArchive(cppqedutils::iarchive& iar) {StateVectorLow temp; iar & temp; psi_.getArray().reference(temp); return iar;}

  /** 
   * structure of QuantumJumpMonteCarlo archives:
   * metaData – array – time – ( odeStepper – odeLogger – dtDid – dtTry ) – randomEngine – logger
   * (state should precede time in order to be compatible with array-only archives)
   */
  auto& stateIO(cppqedutils::iarchive& iar)
  {
    StateVectorLow temp;
    ode_.stateIO(iar & temp & t_) & re_.engine & logger_;
    psi_.getArray().reference(temp);
    t0_=t_; // A very important step!
    return iar;
  }
  
  auto& stateIO(cppqedutils::oarchive& oar) {return ode_.stateIO(oar & psi_.getArray() & t_) & re_.engine & logger_;}

  /// Forwards to ::structure::Averaged::streamKey
  std::ostream& streamKey(std::ostream& os) const {size_t i=3; return ::structure::streamKey<::structure::LA_Av>(sys_,os,i);}

  std::ostream& logOnEnd(std::ostream& os) const {return logger_.onEnd(ode_.logOnEnd(os));} ///< calls qjmc::Logger::onEnd
  
  const StateVector& getStateVector() const {return psi_;}
  
  const qjmc::Logger& getLogger() const {return logger_;}
  
  void referenceNewStateVector(const StateVector& psi) {psi_.reference(psi);}
  void setODE(ODE_Engine ode) {ode_=std::move(ode);}
  
  double sampleRandom() {return distro_(re_.engine);}
  
private:
  typedef std::tuple<int,StateVectorLow> IndexSVL_tuple;

  typedef std::vector<IndexSVL_tuple> IndexSVL_tuples;

  StateVector averaged() const {return psi_;}

  void coherentTimeDevelopment(double Dt, std::ostream&);
  
  const IndexSVL_tuples calculateSpecialRates(::structure::Rates* rates) const;

  bool manageTimeStep (const ::structure::Rates& rates, double tCache, double dtDidCache, std::ostream&, bool logControl=true);

  void performJump (const ::structure::Rates&, const IndexSVL_tuples&, std::ostream&); // LOGICALLY non-const
  // helpers to step---we are deliberately avoiding the normal technique of defining such helpers, because in that case the whole QuantumJumpMonteCarlo has to be passed

  double t_=0., t0_=0.;
  
  const HA& ha_;
  const LI& li_;
  const EV& ev_;
  
  StateVectorView<RANK> psi_;

  ODE_Engine ode_;

  randomutils::EngineWithParameters<RandomEngine> re_;
  
  const double dpLimit_, overshootTolerance_;

  mutable qjmc::Logger logger_;
  
  std::uniform_real_distribution<double> distro_{};

};



namespace qjmc {

template<typename ODE_Engine, typename RandomEngine, typename SV>
auto make(::structure::QuantumSystemPtr<std::decay_t<SV>::N_RANK> sys,
          SV&& state, const Pars<RandomEngine>& p)
{
  return QuantumJumpMonteCarlo<std::decay_t<SV>::N_RANK,ODE_Engine,RandomEngine>{
    sys,std::forward<SV>(state),{initialTimeStep(sys),p},{p.seed,p.prngStream},p.dpLimit,p.overshootTolerance,p.logLevel
  };
}

/// Here, it is very important that psi is taken by const reference, since it has to be copied by value into the individual `QuantumJumpMonteCarlo`s
template<typename ODE_Engine, typename RandomEngine, typename V, typename SYS, typename SV>
auto makeEnsemble(SYS sys, const SV& psi, const Pars<RandomEngine>& p, EntanglementMeasuresSwitch ems)
{
  constexpr auto RANK=std::decay_t<SV>::N_RANK;
  using Single=QuantumJumpMonteCarlo<RANK,ODE_Engine,RandomEngine>;

  std::vector<Single> trajs;
  
  p.logLevel=(p.logLevel>0 ? 1 : p.logLevel); // reduced logging for individual trajectories in an Ensemble

  for (size_t i=0; i<p.nTraj; ++i) {
    trajs.push_back(make<ODE_Engine,RandomEngine>(sys,quantumdata::StateVector<RANK>(psi),p));
    randomutils::incrementForNextStream(p);
  }

  auto av=::structure::castAv(sys);
  
  return cppqedutils::trajectory::Ensemble{trajs,DensityOperatorStreamer<RANK,V>{av,ems},
                                           qjmc::EnsembleLogger{p.nBins,p.nJumpsPerBin},
                                           quantumdata::DensityOperator<RANK>{psi.getDimensions()}};

}
 
} // mcwd


} // quantumtrajectory


template <size_t RANK, typename ODE_Engine, typename RandomEngine>
struct cppqedutils::trajectory::MakeSerializationMetadata<quantumtrajectory::QuantumJumpMonteCarlo<RANK,ODE_Engine,RandomEngine>>
{
  static auto _() {return SerializationMetadata{"CArray","QuantumJumpMonteCarlo",RANK};}
};



template<size_t RANK, typename ODE_Engine, typename RandomEngine>
void quantumtrajectory::QuantumJumpMonteCarlo<RANK,ODE_Engine,RandomEngine>::coherentTimeDevelopment(double Dt, std::ostream& logStream)
{
  if (const auto ha=::structure::castHa(sys_)) {
    ode_.step(Dt,logStream,[this,ha](const StateVectorLow& psi, StateVectorLow& dpsidt, double t) {
      dpsidt=0;
      ha->addContribution(t,psi,dpsidt,t0_);
    },t_,psi_.getArray());
  }
  else {
    double stepToDo=::structure::castLi(sys_) ? std::min(ode_.getDtTry(),Dt) : Dt; // Cf. tracker #3482771
    t_+=stepToDo; ode_.setDtDid(stepToDo);
  }
  
  // This defines three levels:
  // 1. System is Hamiltonian -> internal timestep control is used with dtTry possibly modified by Liouvillian needs (cf. manageTimeStep)
  // 2. System is not Hamiltonian, but it is Liouvillian -> dtTry is used, which is governed solely by Liouvillian in this case (cf. manageTimeStep)
  // 3. System is neither Hamiltonian nor Liouvillian (might be of Exact) -> a step of Dt is taken
  
  if (const auto ex=::structure::castEx(sys_)) {
    ex->actWithU(t_,psi_.getArray(),t0_);
    t0_=t_;
  }

  logger_.processNorm(psi_.renorm());

}


template<size_t RANK, typename ODE_Engine, typename RandomEngine>
auto quantumtrajectory::QuantumJumpMonteCarlo<RANK,ODE_Engine,RandomEngine>::calculateSpecialRates(::structure::Rates* rates) const -> const IndexSVL_tuples
{
  IndexSVL_tuples res;
  for (int i=0; i<rates->size(); i++)
    if ((*rates)(i)<0) {
      StateVector psiTemp(psi_);
      actWithJ(sys_,t_,psiTemp.getArray(),i);
      res.push_back(IndexSVL_tuple(i,psiTemp.getArray()));
      (*rates)(i)=cppqedutils::sqr(psiTemp.renorm());
    } // psiTemp disappears here, but its storage does not, because the ownership is taken over by the StateVectorLow in the tuple
  return res;
}


template<size_t RANK, typename ODE_Engine, typename RandomEngine>
bool quantumtrajectory::QuantumJumpMonteCarlo<RANK,ODE_Engine,RandomEngine>::
manageTimeStep(const ::structure::Rates& rates, double tCache, double dtDidCache, std::ostream& logStream, bool logControl)
{
  const double totalRate=boost::accumulate(rates,0.);
  const double dtDid=getDtDid(), dtTry=ode_.getDtTry();

  const double liouvillianSuggestedDtTry=dpLimit_/totalRate;

  // Assumption: overshootTolerance_>=1 (equality is the limiting case of no tolerance)
  if (totalRate*dtDid>overshootTolerance_*dpLimit_) {
    t_=tCache; ode_.setDtDid(dtDidCache); ode_.setDtTry(liouvillianSuggestedDtTry);
    logger_.stepBack(logStream,totalRate*dtDid,dtDid,ode_.getDtTry(),t_,logControl);
    return true; // Step-back required.
  }
  else if (totalRate*dtTry>dpLimit_) {
    logger_.overshot(logStream,totalRate*dtTry,dtTry,liouvillianSuggestedDtTry,logControl);
  }

  // dtTry-adjustment for next step:
  ode_.setDtTry(::structure::castHa(sys_) ? std::min(ode_.getDtTry(),liouvillianSuggestedDtTry) : liouvillianSuggestedDtTry);

  return false; // Step-back not required.
}


template<size_t RANK, typename ODE_Engine, typename RandomEngine>
void quantumtrajectory::QuantumJumpMonteCarlo<RANK,ODE_Engine,RandomEngine>::performJump(const ::structure::Rates& rates, const IndexSVL_tuples& specialRates, std::ostream& logStream)
{
  double random=sampleRandom()/getDtDid();

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
      ::structure::actWithJ(sys_,t_,psi_.getArray(),lindbladNo);
      double normFactor=sqrt(rates(lindbladNo));
      if (!boost::math::isfinite(normFactor)) throw std::runtime_error("Infinite detected in QuantumJumpMonteCarlo::performJump");
      psi_/=normFactor;
    }

    logger_.jumpOccured(logStream,t_,lindbladNo);
  }
}


template<size_t RANK, typename ODE_Engine, typename RandomEngine>
void quantumtrajectory::QuantumJumpMonteCarlo<RANK,ODE_Engine,RandomEngine>::step(double Dt, std::ostream& logStream)
{
  const StateVector psiCache(psi_); // deep copy
  double tCache=t_, dtDidCache=getDtDid();

  coherentTimeDevelopment(Dt,logStream);

  if (const auto li=::structure::castLi(sys_)) {

    auto rates(li->rates(t_,psi_));
    IndexSVL_tuples specialRates=calculateSpecialRates(&rates);

    while (manageTimeStep(rates,tCache,dtDidCache,logStream)) {
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


template<size_t RANK, typename ODE_Engine, typename RandomEngine>
std::ostream& quantumtrajectory::QuantumJumpMonteCarlo<RANK,ODE_Engine,RandomEngine>::streamParameters(std::ostream& os) const
{
  using namespace std;
  
  ::structure::streamCharacteristics(sys_,sys_->streamParameters(
    re_.stream(ode_.streamParameters(os))<<"\nMCWF Trajectory Parameters: dpLimit="<<dpLimit_<<" (overshoot tolerance factor)="<<overshootTolerance_<<endl<<endl) )<<endl;

  if (const auto li=::structure::castLi(sys_)) {
    os<<"Decay channels:\n";
    {
      size_t i=0;
      li->streamKey(os,i);
    }
    os<<"Alternative Lindblads: ";
    {
      const auto rates(li->rates(0,psi_));
      int n=0;
      for (int i=0; i<rates.size(); ++i) if (rates(i)<0) {os<<i<<' '; ++n;}
      if (!n) os<<"none";
    }
    os<<endl;
  }

  return os;
  
}


template<size_t RANK, typename ODE_Engine, typename RandomEngine>
struct cppqedutils::trajectory::AverageTrajectories<quantumtrajectory::QuantumJumpMonteCarlo<RANK,ODE_Engine,RandomEngine>>
{
  static const auto& _(quantumdata::DensityOperator<RANK>& rho,
                       const std::vector<quantumtrajectory::QuantumJumpMonteCarlo<RANK,ODE_Engine,RandomEngine>>& trajs)
  {
    rho=trajs.begin()->getStateVector();
      
    for (auto i=trajs.begin()+1; i<trajs.end(); i++) i->getStateVector().addTo(rho);

    return rho/=size_t2Double(trajs.size());
    
  }

};


/*
 * This could be implemented in several different ways, depending on how many arrays the archive contains
 * - Single array archive: initialize all trajectories from the same array
 * - As many arrays in archive as trajectories: initialize all trajectories from its own corresponding array
 * **Most general solution**: create a DensityOperator from the available arrays (independently of their number)
 * and sample the density operator to initialize as many trajectories as needed.
 */
template<size_t RANK, typename ODE_Engine, typename RandomEngine>
struct cppqedutils::trajectory::InitializeEnsembleFromArrayOnlyArchive<quantumtrajectory::QuantumJumpMonteCarlo<RANK,ODE_Engine,RandomEngine>>
{
  static auto& _(const std::vector<quantumtrajectory::QuantumJumpMonteCarlo<RANK,ODE_Engine,RandomEngine>>&, cppqedutils::iarchive& iar)
  {
    throw std::runtime_error("InitializeEnsembleFromArrayOnlyArchive not implemented for QuantumJumpMonteCarlo");
    return iar;
  }
};


