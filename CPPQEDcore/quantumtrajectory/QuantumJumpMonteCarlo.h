// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "QJMC_Logger.h"
#include "QuantumTrajectory.h"

#include "StochasticTrajectory.h"


namespace quantumtrajectory {

/// Auxiliary tools to QuantumJumpMonteCarlo  
namespace qjmc {


/// Aggregate of parameters of QuantumJumpMonteCarlo
template<typename RandomEngine>
struct Pars : public cppqedutils::trajectory::ParsStochastic<RandomEngine> {
  
  double dpLimit; ///< the parameter \f$\Delta p\f$

  size_t
    nBins, ///< governs how many bins should be used for the histogram of jumps created by qjmc::EnsembleLogger::stream (a zero value means a heuristic automatic determination)
    nJumpsPerBin; ///< the average number of jumps per bin in the histogram of jumps for the case of heuristic bin-number determination

  Pars(popl::OptionParser& op) : cppqedutils::trajectory::ParsStochastic<RandomEngine>{op}
  {
    using ::parameters::_;
    add(op,"QuantumJumpMonteCarlo",
      _("dpLimit","QJMC stepper total jump probability limit",0.01,dpLimit),
      _("nBins","number of bins used for a histogram of jumps",0,nBins),
      _("nJumpsPerBin","average number of jumps per bin in the histogram of jumps for the case of heuristic bin-number determination",50,nJumpsPerBin)
    );
  }

};


} // qjmc



/// Implements a single Quantum-Jump Monte Carlo trajectory
/**
 * the stepwise adaptive algorithm is used, cf. Comp. Phys. Comm. 238:88 (2019)
 * \note Finite overshoot tolerance is not supported. The extent of overshoots can be found out from the logs anyway
 */
template<
  size_t RANK,
  ::structure::quantum_system_dynamics<RANK> QSD,
  ::cppqedutils::ode::engine<typename ::quantumdata::StateVector<RANK>::StorageType> OE,
  std::uniform_random_bit_generator RandomEngine >
struct QuantumJumpMonteCarlo
{
  using StateVector=quantumdata::StateVector<RANK>;

  using EnsembleAverageElement = StateVector;
  using EnsembleAverageResult = quantumdata::DensityOperator<RANK>;

  using Rates = std::vector<double>;
  
  QuantumJumpMonteCarlo(auto&& qsd, auto&& psi, auto&& oe, randomutils::EngineWithParameters<RandomEngine> re, double dpLimit)
  : qsd{std::forward<decltype(qsd)>(qsd)}, psi{std::forward<decltype(psi)>(psi)}, oe{std::forward<decltype(oe)>(oe)}, re{re}, dpLimit_{dpLimit},
    logger_{std::size(getLi(this->qsd))}
  {
    if (const auto& li{getLi(this->qsd)}; // Careful! the argument shadows the member
        !time && size(li) ) manageTimeStep( calculateRates(li) );
  }

  double time=0., time0=0.;
  QSD qsd;
  StateVector psi;
  OE oe;
  randomutils::EngineWithParameters<RandomEngine> re;

  friend double getDtDid(const QuantumJumpMonteCarlo& q) {return getDtDid(q.oe);}

  friend double getTime(const QuantumJumpMonteCarlo& q) {return q.time;}

  /// TODO: put here the system-specific things
  friend ::cppqedutils::LogTree logIntro(const QuantumJumpMonteCarlo& q)
  {
    return {{"Quantum-Jump Monte Carlo",{{"QJMC algorithm","stepwise"},{"odeEngine",logIntro(q.oe)},{"randomEngine",logIntro(q.re)},{"System","TAKE FROM SYSTEM"}}}};
  }

  friend ::cppqedutils::LogTree logOutro(const QuantumJumpMonteCarlo& q) {return {{"QJMC",q.logger_.outro()},{"ODE_Engine",logOutro(q.oe)}};}

  /// TODO: the returned log should contain information about the coherent step + the time step change as well
  friend auto step(QuantumJumpMonteCarlo& q, double deltaT)
  {
    ::cppqedutils::LogTree res;

    // Coherent time development
    step(q.oe, std::min(q.oe.dtTry,deltaT), [&] (const typename StateVector::StorageType& psiRaw, typename StateVector::StorageType& dpsidtRaw, double t)
    {
      ::structure::applyHamiltonian(getHa(q.qsd),t,
                                    ::quantumdata::StateVectorConstView<RANK>{q.psi.extents,q.psi.strides,0,psiRaw},
                                    ::quantumdata::StateVectorView<RANK>{q.psi.extents,q.psi.strides,0,
                                                                         dpsidtRaw.begin(),std::ranges::fill(dpsidtRaw,0)},
                                    q.time0);
    },q.time,q.psi.dataStorage());

    ::structure::applyPropagator(getHa(q.qsd),q.time,q.psi.mutableView(),q.time0); q.time0=q.time;

    q.logger_.processNorm(renorm(q.psi));

    // Jump
    if (const auto& li{getLi(q.qsd)}; size(li)) {

      auto rates(q.calculateRates(li));

      q.manageTimeStep(rates);

      // perform jump
      double random=q.sampleRandom()/getDtDid(q);

      size_t lindbladNo=0; // TODO: this could be expressed with an iterator into rates
#pragma GCC warning "PROBLEM: What guarantees that rates is a random-access range?"
      for (; random>0 && lindbladNo!=rates.size(); random-=rates[lindbladNo++]) ;

      if (random<0) { // Jump corresponding to Lindblad no. lindbladNo-1 occurs
        ::structure::applyJump(li[--lindbladNo].jump,q.time,q.psi.mutableView());

        double normFactor=sqrt(rates[lindbladNo]);

        if (!boost::math::isfinite(normFactor)) throw std::runtime_error("Infinite detected in QuantumJumpMonteCarlo::performJump");

        for (dcomp& v : q.psi.dataStorage()) v/=normFactor;

        res.emplace("jump",q.logger_.jumpOccured(q.time,lindbladNo));
      }

    }

    q.logger_.step();

    return res;

  }

  friend ::cppqedutils::LogTree dataStreamKey(const QuantumJumpMonteCarlo& q) {return {{"QuantumJumpMonteCarlo","TAKE FROM SYSTEM"}};}

  friend auto temporalDataPoint(const QuantumJumpMonteCarlo& q)
  {
    return ::structure::calculateAndPostprocess<RANK>( getEV(q.qsd), q.time, ::quantumdata::StateVectorConstView<RANK>(q.psi) );
  }

  friend ::cppqedutils::iarchive& readFromArrayOnlyArchive(QuantumJumpMonteCarlo& q, ::cppqedutils::iarchive& iar) {return iar & q.psi;} // MultiArray can be (de)serialized

  /** 
   * structure of QuantumJumpMonteCarlo archives:
   * metaData – array – time – ( odeStepper – odeLogger – dtDid – dtTry ) – randomEngine – logger
   * (state should precede time in order to be compatible with array-only archives)
   */
  template <typename Archive>
  friend auto& stateIO(QuantumJumpMonteCarlo& q, Archive& ar)
  {
    stateIO(q.oe, ar & q.psi & q.time) & q.re.engine & q.logger_;
    q.time0=q.time;
    return ar;
  }
  
  auto calculateRates(const ::structure::liouvillian<RANK> auto& li) const
  {
    Rates res(std::size(li));
    std::ranges::transform(li, res.begin(), [&] (const ::structure::Lindblad<RANK>& l) -> double {return ::structure::calculateRate(l.rate,time,psi); } );
    return res;
    // TODO: with std::ranges::to the following beautiful solution will be possible:
    // return li | std::views::transform([&] (const ::structure::Lindblad<RANK>& l) -> double {return ::structure::calculateRate(l.rate,time,psi); } ) | std::ranges::to<Rates>() ;

  }

  double sampleRandom() {return distro_(re.engine);}

private:
  const double dpLimit_;
  mutable qjmc::Logger logger_;

  std::uniform_real_distribution<double> distro_{};

  StateVector averaged() const {return psi;}

  auto manageTimeStep(const Rates& rates)
  {
    ::cppqedutils::LogTree res;

    const double totalRate=std::accumulate(rates.begin(),rates.end(),0.);
    const double dtDid=getDtDid(*this), dtTry=oe.dtTry;

    const double liouvillianSuggestedDtTry=dpLimit_/totalRate;

    if (totalRate*dtTry>dpLimit_) res=logger_.overshot(totalRate*dtTry,dtTry,liouvillianSuggestedDtTry);

    // dtTry-adjustment for next step:
    oe.dtTry = std::min(oe.dtTry,liouvillianSuggestedDtTry) ;

    return res;
  }

};


template<typename QSD, typename SV, typename OE, typename RandomEngineWithParameters>
QuantumJumpMonteCarlo(QSD , SV , OE , RandomEngineWithParameters , double )
-> QuantumJumpMonteCarlo< quantumdata::multiArrayRank_v<SV>, QSD, OE, typename RandomEngineWithParameters::Engine >;



namespace qjmc {

/*
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
 */
} // qjmc


} // quantumtrajectory


template <size_t RANK, typename QSD, typename ODE_Engine, typename RandomEngine>
struct cppqedutils::trajectory::MakeSerializationMetadata<quantumtrajectory::QuantumJumpMonteCarlo<RANK,QSD,ODE_Engine,RandomEngine>>
{
  static auto _() {return SerializationMetadata{"CArray","QuantumJumpMonteCarlo",RANK};}
};


/*

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

*/

/*
 * This could be implemented in several different ways, depending on how many arrays the archive contains
 * - Single array archive: initialize all trajectories from the same array
 * - As many arrays in archive as trajectories: initialize all trajectories from its own corresponding array
 * **Most general solution**: create a DensityOperator from the available arrays (independently of their number)
 * and sample the density operator to initialize as many trajectories as needed.
 */
/*template<size_t RANK, typename ODE_Engine, typename RandomEngine>
struct cppqedutils::trajectory::InitializeEnsembleFromArrayOnlyArchive<quantumtrajectory::QuantumJumpMonteCarlo<RANK,ODE_Engine,RandomEngine>>
{
  static auto& _(const std::vector<quantumtrajectory::QuantumJumpMonteCarlo<RANK,ODE_Engine,RandomEngine>>&, cppqedutils::iarchive& iar)
  {
    throw std::runtime_error("InitializeEnsembleFromArrayOnlyArchive not implemented for QuantumJumpMonteCarlo");
    return iar;
  }
};*/


