// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "QJMC_Logger.h"
#include "QuantumTrajectory.h"

#include "StochasticTrajectory.h"


namespace quantumtrajectory {

/// Auxiliary tools to QuantumJumpMonteCarlo2
namespace qjmc {


/// Aggregate of parameters of QuantumJumpMonteCarlo2
template<typename RandomEngine>
struct Pars : public trajectory::ParsStochastic<RandomEngine> {
  
  double normTol; ///< the parameter \f$\Delta p\f$

  size_t
    nBins, ///< governs how many bins should be used for the histogram of jumps created by qjmc::EnsembleLogger::stream (a zero value means a heuristic automatic determination)
    nJumpsPerBin; ///< the average number of jumps per bin in the histogram of jumps for the case of heuristic bin-number determination

  Pars(popl::OptionParser& op) : trajectory::ParsStochastic<RandomEngine>{op}
  {
    using ::parameters::_;
    add(op,"QuantumJumpMonteCarlo22",
      _("normTol","QJMC stepper total jump probability limit",1e-6/*0.001*/,normTol),
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
  quantum_system_dynamics<RANK> QSD,
  ode::engine<StorageType> OE,
  std::uniform_random_bit_generator RandomEngine >
struct QuantumJumpMonteCarlo2
{
  using EnsembleAverageElement = const StateVector<RANK>&;
  using EnsembleAverageResult = DensityOperator<RANK>;

  using Rates = std::vector<double>;
  
  QuantumJumpMonteCarlo2(QuantumJumpMonteCarlo2&&) = default;

  QuantumJumpMonteCarlo2(auto&& qsd, auto&& psi, auto&& oe, randomutils::EngineWithParameters<RandomEngine> re, double normTol)
  : qsd{std::forward<decltype(qsd)>(qsd)}, psi{std::forward<decltype(psi)>(psi)}, oe{std::forward<decltype(oe)>(oe)}, re{re}, normTol_{normTol},
    logger_{std::size(getLi(this->qsd))}, normAt_{0.1/*sampleRandom()*/}
  {}

  double time=0., time0=0.;
  QSD qsd;
  StateVector<RANK> psi;
  OE oe;
  randomutils::EngineWithParameters<RandomEngine> re;

  friend double getDtDid(const QuantumJumpMonteCarlo2& q) {return getDtDid(q.oe);}

  friend double getTime(const QuantumJumpMonteCarlo2& q) {return q.time;}

  /// TODO: put here the system-specific things
  friend LogTree logIntro(const QuantumJumpMonteCarlo2& q)
  {
    return {{"Quantum-Jump Monte Carlo",{{"QJMC algorithm","stepwise"},{"odeEngine",logIntro(q.oe)},{"randomEngine",logIntro(q.re)},{"System","TAKE FROM SYSTEM"}}}};
  }

  friend LogTree logOutro(const QuantumJumpMonteCarlo2& q) {return {{"QJMC",q.logger_.outro()},{"ODE_Engine",logOutro(q.oe)}};}

  friend double evolveNorm(QuantumJumpMonteCarlo2& q, double deltaT)
  {

    // Coherent time development
    step(q.oe, deltaT, [&] (const StorageType& psiRaw, StorageType& dpsidtRaw, double t)
    {
      applyHamiltonian(getHa(q.qsd),t,
                       StateVectorConstView<RANK>{q.psi.extents,q.psi.strides,0,psiRaw},
                       StateVectorView<RANK>{q.psi.extents,q.psi.strides,0,dpsidtRaw.begin(),std::ranges::fill(dpsidtRaw,0)},
                       q.time0);
    },q.time,q.psi.dataStorage());

    applyPropagator(getEx(q.qsd),q.time,q.psi.mutableView(),q.time0); q.time0=q.time;

    std::cout<<"# # norm: "<<norm(q.psi)<<" "<<q.normAt_<<" "<<q.time<<std::endl;

    return norm(q.psi);
  }


  friend LogTree step( QuantumJumpMonteCarlo2& q, double deltaT, std::optional<std::pair<double,double>> bisect = std::nullopt)
  {
    LogTree res;

    if (bisect) std::cout<<"# # bisecting: "<<bisect->first<<" "<<bisect->second<<std::endl;

    if (double norm{evolveNorm(q,deltaT)}; abs(norm-q.normAt_)<q.normTol_) q.performJump();
    else if ( norm < q.normAt_-q.normTol_ ) {
     std::pair newBisect{bisect ? bisect->first : q.time-getDtDid(q), q.time} ;
      step( q, (newBisect.first-newBisect.second)/2. , newBisect );
    }
    else if (bisect) {
      // in this case the norm is above the interval, which needs special treatment only when bisecting
      // otherwise it simply means that we haven’t YET reached normAt_
      std::pair newBisect{ q.time, bisect->second } ;
      step( q, (newBisect.second-newBisect.first)/2. , newBisect );
    }

    return res;
  }


  void performJump()
  {
    normAt_=sampleRandom();
    renorm(psi);

    // Jump
    if (const auto& li{getLi(qsd)}; size(li)) {

      auto rates(calculateRates(li));

      // perform jump
      double random=sampleRandom()/getDtDid(oe);

      size_t lindbladNo=0; // TODO: this could be expressed with an iterator into rates

      for (; random>0 && lindbladNo!=rates.size(); random-=rates[lindbladNo++]) ;

      if (random<0) { // Jump corresponding to Lindblad no. lindbladNo-1 occurs
        applyJump(li[--lindbladNo].jump,time,psi.mutableView());

        double normFactor=sqrt(rates[lindbladNo]);

        if (!boost::math::isfinite(normFactor)) throw std::runtime_error("Infinite detected in QuantumJumpMonteCarlo::performJump");

        for (dcomp& v : psi.dataStorage()) v/=normFactor;

        // res.emplace("jump",logger_.jumpOccured(time,lindbladNo));
      }

    }

  }


  friend LogTree dataStreamKey(const QuantumJumpMonteCarlo2& q) {return {{"QuantumJumpMonteCarlo2","TAKE FROM SYSTEM"}};}

  friend auto temporalDataPoint(const QuantumJumpMonteCarlo2& q)
  {
    // CAREFUL!!! tdp should be normalized, since the state vector is not!
    return calculateAndPostprocess<RANK>( getEV(q.qsd), q.time, LDO<StateVector,RANK>(q.psi) );
  }

  friend iarchive& readFromArrayOnlyArchive(QuantumJumpMonteCarlo2& q, iarchive& iar) {return iar & q.psi;} // MultiArray can be (de)serialized

  /** 
   * structure of QuantumJumpMonteCarlo2 archives:
   * metaData – array – time – ( odeStepper – odeLogger – dtDid – dtTry ) – randomEngine – logger
   * (state should precede time in order to be compatible with array-only archives)
   */
  template <typename Archive>
  friend auto& stateIO(QuantumJumpMonteCarlo2& q, Archive& ar)
  {
    stateIO(q.oe, ar & q.psi & q.time) & q.re.engine & q.logger_;
    q.time0=q.time;
    return ar;
  }
  
  auto calculateRates(const Liouvillian<RANK> & li) const
  {
    Rates res(std::size(li));
    std::ranges::transform(li, res.begin(), [&] (const Lindblad<RANK>& l) -> double {return calculateRate(l.rate,time,psi); } );
    return res;
    // TODO: with std::ranges::to the following beautiful solution will be possible:
    // return li | std::views::transform([&] (const Lindblad<RANK>& l) -> double {return calculateRate(l.rate,time,psi); } ) | std::ranges::to<Rates>() ;

  }

  double sampleRandom() {return distro_(re.engine);}

  friend const StateVector<RANK>& averaged(const QuantumJumpMonteCarlo2& q) {return q.psi;}

private:
  const double normTol_;
  double normAt_;
  mutable qjmc::Logger logger_;

  std::uniform_real_distribution<double> distro_{};

};


template<typename QSD, typename SV, typename OE, typename RandomEngineWithParameters>
QuantumJumpMonteCarlo2(QSD , SV , OE , RandomEngineWithParameters , double )
-> QuantumJumpMonteCarlo2< multiArrayRank_v<SV>, QSD, OE, typename RandomEngineWithParameters::Engine >;



namespace qjmc {

template <template<typename> class OE, typename QSD, typename SV, typename RandomEngine>
auto make(QSD&& qsd, SV&& state, const Pars<RandomEngine>& p)
{
  constexpr size_t RANK=multiArrayRank_v<std::decay_t<SV>>;
  using ODE=OE<StorageType>;

  double iDt=initialTimeStep(getFreqs(qsd)); // precalculate, since qsd gets forwarded (i.e., potentially moved)
  
  return QuantumJumpMonteCarlo2<RANK,QSD,ODE,RandomEngine>{
    std::forward<QSD>(qsd),
    std::forward<SV>(state),
    ODE{iDt,p.epsRel,p.epsAbs},
    randomutils::EngineWithParameters<RandomEngine>{p.seed,p.prngStream},
    p.normTol};
}


/// Here, it is very important that psi is taken by const reference, since it has to be copied by value into the individual `QuantumJumpMonteCarlo2`s
template<template<typename> class OE, /* auto retainedAxes,*/ size_t RANK,
         quantum_system_dynamics<RANK> QSD, typename RandomEngine>
auto makeEnsemble(QSD&& qsd, const StateVector<RANK>& psi, Pars<RandomEngine>& p /*, EntanglementMeasuresSwitch ems*/)
{
  using QSD_const_ref = std::add_lvalue_reference_t<std::add_const_t<QSD>>;
  // quantum_system_dynamics is stored by value in TDP_DensityOperator, and the individual trajectories borrow it from there
  using Single=QuantumJumpMonteCarlo2<RANK,QSD_const_ref,OE<StorageType>,RandomEngine>;

  trajectory::Ensemble<Single,TDP_DensityOperator<RANK,QSD>> res{
    .trajs{},
    .tdpCalculator{std::forward<QSD>(qsd)/*,ems*/},
    // qjmc::EnsembleLogger{p.nBins,p.nJumpsPerBin},
    .ensembleAverageResult{DensityOperator<RANK>{getDimensions(psi),noInit}}
  };
  
  for (size_t i=0; i<p.nTraj; ++i) {
    res.trajs.push_back(make<OE,QSD_const_ref>(res.tdpCalculator.qsd,StateVector<RANK>{copy(psi)},p));
    randomutils::incrementForNextStream(p);
  }

  return res;

}

} // qjmc


} // quantumtrajectory


template <size_t RANK, typename QSD, typename ODE_Engine, typename RandomEngine>
struct cppqedutils::trajectory::MakeSerializationMetadata<quantumtrajectory::QuantumJumpMonteCarlo2<RANK,QSD,ODE_Engine,RandomEngine>>
{
  static auto _() {return SerializationMetadata{"CArray","QuantumJumpMonteCarlo2",RANK};}
};



template<size_t RANK, typename QSD, typename ODE_Engine, typename RandomEngine>
struct cppqedutils::trajectory::AverageTrajectories<quantumtrajectory::QuantumJumpMonteCarlo2<RANK,QSD,ODE_Engine,RandomEngine>>
{
  static const auto& _(quantumdata::DensityOperator<RANK>& rho, const auto& trajs)
  {
    std::ranges::for_each(trajs | std::views::drop(1) ,
                          [&rhoInner=(rho=dyad(trajs[0].psi,trajs[0].psi))] (const auto& traj) {addTo(traj.psi,rhoInner);});
    for (dcomp& v : rho.dataStorage()) v/=dcomp(size(trajs));
    return rho;
  }

};


/*
 * This could be implemented in several different ways, depending on how many arrays the archive contains
 * - Single array archive: initialize all trajectories from the same array
 * - As many arrays in archive as trajectories: initialize all trajectories from its own corresponding array
 * **Most general solution**: create a DensityOperator from the available arrays (independently of their number)
 * and sample the density operator to initialize as many trajectories as needed.
 */
template<size_t RANK, typename QSD, typename ODE_Engine, typename RandomEngine>
struct cppqedutils::trajectory::InitializeEnsembleFromArrayOnlyArchive<quantumtrajectory::QuantumJumpMonteCarlo2<RANK,QSD,ODE_Engine,RandomEngine>>
{
  static auto& _(const std::vector<quantumtrajectory::QuantumJumpMonteCarlo2<RANK,QSD,ODE_Engine,RandomEngine>>&, iarchive& iar)
  {
    throw std::runtime_error("InitializeEnsembleFromArrayOnlyArchive not implemented for QuantumJumpMonteCarlo2");
    return iar;
  }
};




// template<size_t RANK, typename ODE_Engine, typename RandomEngine>
// std::ostream& quantumtrajectory::QuantumJumpMonteCarlo2<RANK,ODE_Engine,RandomEngine>::streamParameters(std::ostream& os) const
// {
//   using namespace std;
//
//   streamCharacteristics(sys_,sys_->streamParameters(
//     re_.stream(ode_.streamParameters(os))<<"\nMCWF Trajectory Parameters: dpLimit="<<dpLimit_<<" (overshoot tolerance factor)="<<overshootTolerance_<<endl<<endl) )<<endl;
//
//   if (const auto li=castLi(sys_)) {
//     os<<"Decay channels:\n";
//     {
//       size_t i=0;
//       li->streamKey(os,i);
//     }
//     os<<"Alternative Lindblads: ";
//     {
//       const auto rates(li->rates(0,psi_));
//       int n=0;
//       for (int i=0; i<rates.size(); ++i) if (rates(i)<0) {os<<i<<' '; ++n;}
//       if (!n) os<<"none";
//     }
//     os<<endl;
//   }
//
//   return os;
//
// }

