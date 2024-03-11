// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "QuantumTrajectory.h"

#include "StochasticTrajectory.h"


namespace quantumtrajectory {

/// Auxiliary tools to QuantumJumpMonteCarlo  
namespace qjmc {


enum class Algorithm { integrating, stepwise } ;


/// Aggregate of parameters of QuantumJumpMonteCarlo
template<typename RandomEngine>
struct ParsBase : public trajectory::ParsStochastic<RandomEngine> {
  size_t
    nBins, ///< governs how many bins should be used for the histogram of jumps created by qjmc::EnsembleLogger::stream (a zero value means a heuristic automatic determination)
    nJumpsPerBin; ///< the average number of jumps per bin in the histogram of jumps for the case of heuristic bin-number determination

  ParsBase(popl::OptionParser& op) : trajectory::ParsStochastic<RandomEngine>{op}
  {
    using ::parameters::_;
    add(op,"QuantumJumpMonteCarlo",
      _("nBins","number of bins used for a histogram of jumps",0,nBins),
      _("nJumpsPerBin","average number of jumps per bin in the histogram of jumps for the case of heuristic bin-number determination",50,nJumpsPerBin)
    );
  }

};


template<Algorithm, typename> struct Pars;

template<typename RandomEngine>
struct Pars<Algorithm::stepwise,RandomEngine> : public ParsBase<RandomEngine> {
  double dpLimit; ///< the parameter \f$\Delta p\f$
  Pars(popl::OptionParser& op) : ParsBase<RandomEngine>{op} {add(op,::parameters::_("dpLimit","QJMC stepper total jump probability limit",0.01,dpLimit));}
};


template<typename RandomEngine>
struct Pars<Algorithm::integrating,RandomEngine> : public ParsBase<RandomEngine> {

  double normTol;

  Pars(popl::OptionParser& op) : ParsBase<RandomEngine>{op}
  {
    using ::parameters::_;
    add(op,
      _("normTol","QJMC stepper total jump probability limit",1e-6/*0.001*/,normTol)
    );
  }

};


LogTree defaultLogger()
{
  return {{"Total number of MCWF steps",0z},{"Jump trajectory",json::array{}}};
}


} // qjmc


/// Implements a single Quantum-Jump Monte Carlo trajectory
template <qjmc::Algorithm, size_t RANK, quantum_system_dynamics<RANK>, ode::engine<StorageType>, std::uniform_random_bit_generator>
struct QuantumJumpMonteCarlo;


using Rates = std::vector<double>;



template<
  size_t RANK,
  quantum_system_dynamics<RANK> QSD,
  ode::engine<StorageType> OE,
  std::uniform_random_bit_generator RandomEngine >
struct QuantumJumpMonteCarloBase
{
  using EnsembleAverageElement = const StateVector<RANK>&;
  using EnsembleAverageResult = DensityOperator<RANK>;

  QuantumJumpMonteCarloBase(QuantumJumpMonteCarloBase&&) = default;

  QuantumJumpMonteCarloBase(auto&& qsd, auto&& psi, auto&& oe, randomutils::EngineWithParameters<RandomEngine> re)
  : qsd{std::forward<decltype(qsd)>(qsd)}, psi{std::forward<decltype(psi)>(psi)}, oe{std::forward<decltype(oe)>(oe)}, re{re},
    log_{qjmc::defaultLogger()},
    intro_{{"Quantum-Jump Monte Carlo",{{"odeEngine",logIntro(this->oe)},{"randomEngine",logIntro(this->re)},{"System","TAKE FROM SYSTEM"}}}}
    {}

  double time=0., time0=0.;
  QSD qsd;
  StateVector<RANK> psi;
  OE oe;
  randomutils::EngineWithParameters<RandomEngine> re;

  friend double getDtDid(const QuantumJumpMonteCarloBase& q) {return getDtDid(q.oe);}

  friend double getTime(const QuantumJumpMonteCarloBase& q) {return q.time;}

  /// TODO: put here the system-specific things
  friend LogTree logIntro(const QuantumJumpMonteCarloBase& q)
  {
    return q.intro_;
  }

  friend LogTree logOutro(const QuantumJumpMonteCarloBase& q) {return {{"QJMC",q.log_},{"ODE_Engine",logOutro(q.oe)}};}

  friend LogTree dataStreamKey(const QuantumJumpMonteCarloBase& q) {return {{"QuantumJumpMonteCarloBase","TAKE FROM SYSTEM"}};}

  friend iarchive& readFromArrayOnlyArchive(QuantumJumpMonteCarloBase& q, iarchive& iar) {return iar & q.psi;} // MultiArray can be (de)serialized

  /** 
   * structure of QuantumJumpMonteCarlo archives:
   * metaData – array – time – ( odeStepper – odeLogger – dtDid – dtTry ) – randomEngine – logger
   * (state should precede time in order to be compatible with array-only archives)
   */
  template <typename Archive>
  friend auto& stateIO(QuantumJumpMonteCarloBase& q, Archive& ar)
  {
    stateIO(q.oe, ar & q.psi & q.time) & q.re.engine & q.log_;
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

  friend const StateVector<RANK>& averaged(const QuantumJumpMonteCarloBase& q) {return q.psi;}

protected:
  mutable LogTree log_; // serialization can be solved via converting to/from string

  LogTree intro_;

  std::uniform_real_distribution<double> distro_{};

};


/// the integrating algorithm, cf. https://qutip.org/docs/latest/guide/dynamics/dynamics-monte.html
template<
  size_t RANK,
  quantum_system_dynamics<RANK> QSD,
  ode::engine<StorageType> OE,
  std::uniform_random_bit_generator RandomEngine >
struct QuantumJumpMonteCarlo<qjmc::Algorithm::integrating,RANK,QSD,OE,RandomEngine> : QuantumJumpMonteCarloBase<RANK,QSD,OE,RandomEngine>
{
  QuantumJumpMonteCarlo(QuantumJumpMonteCarlo&&) = default;

  QuantumJumpMonteCarlo(auto&& qsd, auto&& psi, auto&& oe, randomutils::EngineWithParameters<RandomEngine> re, double normTol)
  : QuantumJumpMonteCarloBase<RANK,QSD,OE,RandomEngine>{std::forward<decltype(qsd)>(qsd),std::forward<decltype(psi)>(psi),std::forward<decltype(oe)>(oe),re},
    normTol_{normTol}, normAt_{this->sampleRandom()}
  {
    this->log_["bisectMaxIter"]=0z;
    {
      auto& intro=this->intro_["Quantum-Jump Monte Carlo"].as_object();
      intro["algorithm"]="integrating"; intro["normTol"]=normTol;
    }
  }


  friend double evolveNorm(QuantumJumpMonteCarlo& q, double deltaT)
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

    return norm(q.psi);
  }


  friend LogTree step( QuantumJumpMonteCarlo& q, double deltaT, std::optional<std::tuple<double,double,long>> bisect = std::nullopt)
  {
    LogTree res;

//    if (bisect) std::cout<<"# # bisecting: "<<std::get<0>(*bisect)<<" "<<std::get<1>(*bisect)<<" "<<std::get<2>(*bisect)<<std::endl;

    if (double norm{evolveNorm(q,deltaT)}; abs(norm-q.normAt_)<q.normTol_) {
      res=q.performJump();
      if (bisect) {auto& l=q.log_["bisectMaxIter"].as_int64(); l=std::max(l,std::get<2>(*bisect));}
    }
    else if ( norm < q.normAt_-q.normTol_ ) {
      std::tuple newBisect{bisect ? std::get<0>(*bisect) : q.time-getDtDid(q), q.time, bisect ? ++std::get<2>(*bisect) : 1 } ;
      // TODO: simple binary bisecting is very dumb here, we should use at least a linear interpolation
      // alternatively, more than two points could be passed on, to use arbitrary nonlinear interpolation
      // TODO: furthermore, bisecting results in very small timesteps, the original dtTry should be restored after bisecting is over
      step( q, (std::get<0>(newBisect)-std::get<1>(newBisect))/2. , newBisect );
    }
    else if (bisect) {
      // in this case the norm is above the interval, which needs special treatment only when bisecting
      // otherwise it simply means that we haven’t YET reached normAt_
      std::tuple newBisect{ q.time, std::get<1>(*bisect), ++std::get<2>(*bisect) } ;
      step( q, (std::get<1>(newBisect)-std::get<0>(newBisect))/2. , newBisect );
    }

    return res;
  }


  LogTree performJump()
  {
    LogTree res;

    normAt_=this->sampleRandom();
    renorm(this->psi);

    // Jump
    if (const auto& li{getLi(this->qsd)}; size(li)) {

      auto rates(this->calculateRates(li));

      // perform jump
      double random=std::accumulate(rates.begin(),rates.end(),0.)*this->sampleRandom() ;

      size_t lindbladNo=0; // TODO: this could be expressed with an iterator into rates

      for (; random>0 && lindbladNo!=rates.size(); random-=rates[lindbladNo++]) ;

      if (random<0) { // Jump corresponding to Lindblad no. lindbladNo-1 occurs
        applyJump(li[--lindbladNo].jump,this->time,this->psi.mutableView());

        double normFactor=sqrt(rates[lindbladNo]);

        if (!boost::math::isfinite(normFactor)) throw std::runtime_error("Infinite detected in QuantumJumpMonteCarlo::performJump");

        for (dcomp& v : this->psi.dataStorage()) v/=normFactor;

        res.emplace("jump",LogTree{{"no.",lindbladNo},{"at time",this->time}});
        this->log_["Jump trajectory"].as_array().push_back({this->time,lindbladNo});
      }
//      else throw std::runtime_error("Jump did not occur!") ;

    }

    return res;

  }


  friend auto temporalDataPoint(const QuantumJumpMonteCarlo& q)
  {
    auto res{calculateAndPostprocess<RANK>( getEV(q.qsd), q.time, LDO<StateVector,RANK>(q.psi) )};
    renormTDP(res,norm(q.psi));
    return hana::make_tuple(res,norm(q.psi));
  }


  template <typename Archive>
  friend auto& stateIO(QuantumJumpMonteCarlo& q, Archive& ar)
  {
    return stateIO(static_cast<QuantumJumpMonteCarloBase<RANK,QSD,OE,RandomEngine>&>(q),ar) & q.normAt_;
  }


private:
  const double normTol_;
  double normAt_;

};



/**
 * the stepwise adaptive algorithm is used, cf. Comp. Phys. Comm. 238:88 (2019)
 * \note Finite overshoot tolerance is not supported. The extent of overshoots can be found out from the logs anyway
 */
template<
  size_t RANK,
  quantum_system_dynamics<RANK> QSD,
  ode::engine<StorageType> OE,
  std::uniform_random_bit_generator RandomEngine >
struct QuantumJumpMonteCarlo<qjmc::Algorithm::stepwise,RANK,QSD,OE,RandomEngine> : QuantumJumpMonteCarloBase<RANK,QSD,OE,RandomEngine>
{
  QuantumJumpMonteCarlo(QuantumJumpMonteCarlo&&) = default;

  QuantumJumpMonteCarlo(auto&& qsd, auto&& psi, auto&& oe, randomutils::EngineWithParameters<RandomEngine> re, double dpLimit)
    : QuantumJumpMonteCarloBase<RANK,QSD,OE,RandomEngine>{std::forward<decltype(qsd)>(qsd),std::forward<decltype(psi)>(psi),std::forward<decltype(oe)>(oe),re},
      dpLimit_{dpLimit}
  {
    if (const auto& li{getLi(this->qsd)}; // Careful! the argument shadows the member
        !time && size(li) ) manageTimeStep( this->calculateRates(li) );
    this->log_["dpLimit overshot"]={{"n",0z},{"max.",0.}};
     {
      auto& intro=this->intro_["Quantum-Jump Monte Carlo"].as_object();
      intro["algorithm"]="stepwise"; intro["pdLimit"]=dpLimit;
     }
  }

  /// TODO: the returned log should contain information about the coherent step + the time step change as well
  friend auto step(QuantumJumpMonteCarlo& q, double deltaT)
  {
    LogTree res;

    // Coherent time development
    step(q.oe, deltaT, [&] (const StorageType& psiRaw, StorageType& dpsidtRaw, double t)
    {
      applyHamiltonian(getHa(q.qsd),t,
                       StateVectorConstView<RANK>{q.psi.extents,q.psi.strides,0,psiRaw},
                       StateVectorView<RANK>{q.psi.extents,q.psi.strides,0,dpsidtRaw.begin(),std::ranges::fill(dpsidtRaw,0)},
                       q.time0);
    },q.time,q.psi.dataStorage());

    applyPropagator(getEx(q.qsd),q.time,q.psi.mutableView(),q.time0); q.time0=q.time;

    renorm(q.psi);

    // Jump
    if (const auto& li{getLi(q.qsd)}; size(li)) {

      auto rates(q.calculateRates(li));

      q.manageTimeStep(rates);

      // perform jump
      double random=q.sampleRandom()/getDtDid(q);

      size_t lindbladNo=0; // TODO: this could be expressed with an iterator into rates

      for (; random>0 && lindbladNo!=rates.size(); random-=rates[lindbladNo++]) ;

      if (random<0) { // Jump corresponding to Lindblad no. lindbladNo-1 occurs
        applyJump(li[--lindbladNo].jump,q.time,q.psi.mutableView());

        double normFactor=sqrt(rates[lindbladNo]);

        if (!boost::math::isfinite(normFactor)) throw std::runtime_error("Infinite detected in QuantumJumpMonteCarlo::performJump");

        for (dcomp& v : q.psi.dataStorage()) v/=normFactor;

        res.emplace("jump",LogTree{{"no.",lindbladNo},{"at time",q.time}});
        q.log_["Jump trajectory"].as_array().push_back({q.time,lindbladNo});
      }

    }

    q.log_["Total number of MCWF steps"].as_int64()+=1;

    return res;

  }

  friend auto temporalDataPoint(const QuantumJumpMonteCarlo& q)
  {
    return calculateAndPostprocess<RANK>( getEV(q.qsd), q.time, LDO<StateVector,RANK>(q.psi) );
  }

private:
  const double dpLimit_;

  auto manageTimeStep(const Rates& rates)
  {
    LogTree res;

    const double totalRate=std::accumulate(rates.begin(),rates.end(),0.);
    const double dtDid=getDtDid(*this), dtTry=this->oe.dtTry;

    const double liouvillianSuggestedDtTry=dpLimit_/totalRate;

    if (double dp=totalRate*dtTry; dp>dpLimit_) {
      auto& l=this->log_["dpLimit overshot"].as_object();
      l["n"].as_int64()+=1; l["max."].as_double()=std::max(l["max."].as_double(),dp);
      res={{"dpLimit overshot",dp},{"timestep decreased",{{"from",dtTry},{"to",liouvillianSuggestedDtTry}}}};
    }

    // dtTry-adjustment for next step:
    this->oe.dtTry = std::min(this->oe.dtTry,liouvillianSuggestedDtTry) ;

    return res;
  }

};



namespace qjmc {

template <Algorithm a, template<typename> class OE, typename QSD, typename SV, typename RandomEngine>
auto make(QSD&& qsd, SV&& state, const Pars<a,RandomEngine>& p)
{
  constexpr size_t RANK=multiArrayRank_v<std::decay_t<SV>>;
  using ODE=OE<StorageType>;

  double iDt=initialTimeStep(getFreqs(qsd)); // precalculate, since qsd gets forwarded (i.e., potentially moved)
  
  return QuantumJumpMonteCarlo<a,RANK,QSD,ODE,RandomEngine>{
    std::forward<QSD>(qsd),
    std::forward<SV>(state),
    ODE{iDt,p.epsRel,p.epsAbs},
    randomutils::EngineWithParameters<RandomEngine>{p.seed,p.prngStream},
    [&] () -> double {
      if constexpr (a==Algorithm::stepwise) return p.dpLimit;
      else return p.normTol;
    } () };
}


/// Here, it is very important that psi is taken by const reference, since it has to be copied by value into the individual `QuantumJumpMonteCarlo`s
template<Algorithm a,
         template<typename> class OE, /* auto retainedAxes,*/ size_t RANK,
         quantum_system_dynamics<RANK> QSD, typename RandomEngine>
auto makeEnsemble(QSD&& qsd, const StateVector<RANK>& psi, Pars<a,RandomEngine>& p /*, EntanglementMeasuresSwitch ems*/)
{
  using QSD_const_ref = std::add_lvalue_reference_t<std::add_const_t<QSD>>;
  // quantum_system_dynamics is stored by value in TDP_DensityOperator, and the individual trajectories borrow it from there
  using Single=QuantumJumpMonteCarlo<a,RANK,QSD_const_ref,OE<StorageType>,RandomEngine>;

  trajectory::Ensemble<Single,TDP_DensityOperator<RANK,QSD>> res{
    .trajs{},
    .tdpCalculator{std::forward<QSD>(qsd)/*,ems*/},
    // qjmc::EnsembleLogger{p.nBins,p.nJumpsPerBin},
    .ensembleAverageResult{DensityOperator<RANK>{getDimensions(psi),noInit}}
  };

  for (size_t i=0; i<p.nTraj; ++i) {
    res.trajs.push_back(make<a,OE,QSD_const_ref>(res.tdpCalculator.qsd,StateVector<RANK>{copy(psi)},p));
    randomutils::incrementForNextStream(p);
  }

  return res;

}


} // qjmc


} // quantumtrajectory


template <quantumtrajectory::qjmc::Algorithm a, size_t RANK, typename QSD, typename ODE_Engine, typename RandomEngine>
struct cppqedutils::trajectory::MakeSerializationMetadata<quantumtrajectory::QuantumJumpMonteCarlo<a,RANK,QSD,ODE_Engine,RandomEngine>>
{
  static auto _() {return SerializationMetadata{"CArray","QuantumJumpMonteCarlo",RANK};}
};



template<quantumtrajectory::qjmc::Algorithm a, size_t RANK, typename QSD, typename ODE_Engine, typename RandomEngine>
struct cppqedutils::trajectory::AverageTrajectories<quantumtrajectory::QuantumJumpMonteCarlo<a,RANK,QSD,ODE_Engine,RandomEngine>>
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
template<quantumtrajectory::qjmc::Algorithm a, size_t RANK, typename QSD, typename ODE_Engine, typename RandomEngine>
struct cppqedutils::trajectory::InitializeEnsembleFromArrayOnlyArchive<quantumtrajectory::QuantumJumpMonteCarlo<a,RANK,QSD,ODE_Engine,RandomEngine>>
{
  static auto& _(const std::vector<quantumtrajectory::QuantumJumpMonteCarlo<a,RANK,QSD,ODE_Engine,RandomEngine>>&, iarchive& iar)
  {
    throw std::runtime_error("InitializeEnsembleFromArrayOnlyArchive not implemented for QuantumJumpMonteCarlo");
    return iar;
  }
};




// template<size_t RANK, typename ODE_Engine, typename RandomEngine>
// std::ostream& quantumtrajectory::QuantumJumpMonteCarlo<RANK,ODE_Engine,RandomEngine>::streamParameters(std::ostream& os) const
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


namespace quantumtrajectory::qjmc {


struct EnsembleLogger
{
  typedef std::list<LogTree> LoggerList;

  /// Called by Ensemble<QuantumJumpMonteCarlo>::logOnEnd, it calculates a temporal histogram of total jumps
  /** \todo The different kinds of jumps should be collected into different histograms */
  static std::ostream& stream(std::ostream&, const LoggerList&, size_t nBins, size_t nJumpsPerBin);

  template <typename SingleTrajectory>
  auto& operator()(const std::vector<SingleTrajectory>& trajs, std::ostream& os) const
  {
    LoggerList loggerList;
    for (auto& traj : trajs) loggerList.push_back(traj.getLogger());

    return stream(os,loggerList,nBins_,nJumpsPerBin_);
  }

  const size_t nBins_, nJumpsPerBin_;

};


} // quantumtrajectory::qjmc


/*

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/density.hpp>
#include <boost/accumulators/statistics/stats.hpp>

ostream& quantumtrajectory::qjmc::EnsembleLogger::stream(ostream& os, const LoggerList& loggerList, size_t nBins, size_t nJumpsPerBin)
{
  using namespace boost;

#define AVERAGE_function(f) accumulate(loggerList.begin(),loggerList.end(), 0., [] (double init, const Logger& l) {return init + l.f;})/loggerList.size()
#define MAX_function(f) ranges::max_element(loggerList, [] (const Logger& a, const Logger& b) {return a.f<b.f;})->f

  os<<"\nAverage number of total steps: "<<AVERAGE_function(nMCWF_steps_)<<endl
    <<"\nOn average, dpLimit overshot: "<<AVERAGE_function(nOvershot_)<<" times, maximal overshoot: "<<MAX_function(dpMaxOvershoot_)
    <<"\nOn average, dpTolerance overshot: "<<AVERAGE_function(nToleranceOvershot_)<<" times, maximal overshoot: "<<MAX_function(dpToleranceMaxOvershoot_)<<endl
    <<"\nMaximal deviation of norm from 1: "<<MAX_function(normMaxDeviation_)<<endl;
// NEEDS_WORK: this was before factoring out logging to Evolved/Adaptive, so it could be restored on some more fundamental level:  if (loggerList.front().isHamiltonian_)
//    os<<"\nAverage number of failed ODE steps: "<<AVERAGE_function(nFailedSteps_)<<"\nAverage number of Hamiltonian calls:"<<AVERAGE_function(nHamiltonianCalls_)<<endl;

#undef  MAX_function
#undef  AVERAGE_function

  using namespace accumulators;

  size_t nTotalJumps=accumulate(loggerList.begin(), loggerList.end(), 0, [] (double init, const Logger& l) {return init + l.traj_.size();});

  // Heuristic: if nBins is not given (0), then the number of bins is determined such that the bins contain nJumpsPerBin samples on average
  if (!nBins && nTotalJumps<2*nJumpsPerBin) {
    os<<"\nToo few jumps for a histogram\n";
  }
  else {
    size_t actualNumberOfBins=nBins ? nBins : nTotalJumps/nJumpsPerBin;
    accumulator_set<double, features<tag::density> > acc( tag::density::num_bins = actualNumberOfBins , tag::density::cache_size = nTotalJumps );

    // fill accumulator
    for (auto i : loggerList) for (auto j : i.traj_) acc(j.first);

    const iterator_range<std::vector<std::pair<double, double> >::iterator> histogram=density(acc);

    os<<"\nHistogram of jumps. Number of bins="<<actualNumberOfBins+2<<endl;
    for (auto i : histogram)
      os<<i.first<<"\t"<<i.second<<endl;
  }

  os<<"\nAverage number of total jumps: "<<nTotalJumps/double(loggerList.size())<<endl;

  return os;
}

*/


