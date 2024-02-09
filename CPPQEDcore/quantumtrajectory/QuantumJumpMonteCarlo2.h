// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "QuantumJumpMonteCarlo.h"


namespace quantumtrajectory {


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
struct QuantumJumpMonteCarlo2 : QuantumJumpMonteCarloBase<RANK,QSD,OE,RandomEngine>
{
  QuantumJumpMonteCarlo2(QuantumJumpMonteCarlo2&&) = default;

  QuantumJumpMonteCarlo2(auto&& qsd, auto&& psi, auto&& oe, randomutils::EngineWithParameters<RandomEngine> re, double normTol)
  : QuantumJumpMonteCarloBase<RANK,QSD,OE,RandomEngine>{std::forward<decltype(qsd)>(qsd),std::forward<decltype(psi)>(psi),std::forward<decltype(oe)>(oe),re},
    normTol_{normTol}, normAt_{this->sampleRandom()}
  {
    this->log_["bisectMaxIter"]=0uz;
  }


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

//    std::cout<<"# # norm: "<<norm(q.psi)<<" "<<q.normAt_<<" "<<q.time<<std::endl;

    return norm(q.psi);
  }


  friend LogTree step( QuantumJumpMonteCarlo2& q, double deltaT, std::optional<std::tuple<double,double,size_t>> bisect = std::nullopt)
  {
    LogTree res;

//    if (bisect) std::cout<<"# # bisecting: "<<std::get<0>(*bisect)<<" "<<std::get<1>(*bisect)<<" "<<std::get<2>(*bisect)<<std::endl;

    if (double norm{evolveNorm(q,deltaT)}; abs(norm-q.normAt_)<q.normTol_) {
      res=q.performJump();
      if (bisect) {auto& l=q.log_["bisectMaxIter"].as_uint64(); l=std::max(l,std::get<2>(*bisect));}
    }
    else if ( norm < q.normAt_-q.normTol_ ) {
      std::tuple newBisect{bisect ? std::get<0>(*bisect) : q.time-getDtDid(q), q.time, bisect ? ++std::get<2>(*bisect) : 1 } ;
      // TODO: simple binary bisecting is very dumb here, we should use at least a linear interpolation
      // alternatively, more than two points could be passed on, to use arbitrary nonlinear interpolation
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
      double random=this->sampleRandom() ;

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

    }

    return res;

  }


  friend auto temporalDataPoint(const QuantumJumpMonteCarlo2& q)
  {
    auto res{calculateAndPostprocess<RANK>( getEV(q.qsd), q.time, LDO<StateVector,RANK>(q.psi) )};
    renormTDP(res,norm(q.psi));
    return res;
  }


  template <typename Archive>
  friend auto& stateIO(QuantumJumpMonteCarlo2& q, Archive& ar)
  {
    return stateIO(static_cast<QuantumJumpMonteCarloBase<RANK,QSD,OE,RandomEngine>&>(q),ar) & q.normAt_;
  }
  

private:
  const double normTol_;
  double normAt_;

};


template<typename QSD, typename SV, typename OE, typename RandomEngineWithParameters>
QuantumJumpMonteCarlo2(QSD , SV , OE , RandomEngineWithParameters , double )
-> QuantumJumpMonteCarlo2< multiArrayRank_v<SV>, QSD, OE, typename RandomEngineWithParameters::Engine >;



namespace qjmc {

template <template<typename> class OE, typename QSD, typename SV, typename RandomEngine>
auto make2(QSD&& qsd, SV&& state, const Pars2<RandomEngine>& p)
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



} // qjmc


} // quantumtrajectory


template <size_t RANK, typename QSD, typename ODE_Engine, typename RandomEngine>
struct cppqedutils::trajectory::MakeSerializationMetadata<quantumtrajectory::QuantumJumpMonteCarlo2<RANK,QSD,ODE_Engine,RandomEngine>>
{
  static auto _() {return SerializationMetadata{"CArray","QuantumJumpMonteCarlo2",RANK};}
};


