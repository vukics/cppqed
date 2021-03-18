// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDCORE_QUANTUMTRAJECTORY_MCWF_TRAJECTORY_TCC_INCLUDED
#define   CPPQEDCORE_QUANTUMTRAJECTORY_MCWF_TRAJECTORY_TCC_INCLUDED

#include "MCWF_Trajectory.h"

#include "QuantumTrajectory.h"

#include "FormDouble.h"

#include <boost/range/numeric.hpp>
#include <boost/range/algorithm/find_if.hpp>


//////////////////
//
// Implementations
//
//////////////////


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



#endif // CPPQEDCORE_QUANTUMTRAJECTORY_MCWF_TRAJECTORY_TCC_INCLUDED
