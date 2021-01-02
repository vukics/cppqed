// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDCORE_QUANTUMTRAJECTORY_MCWF_TRAJECTORY_TCC_INCLUDED
#define   CPPQEDCORE_QUANTUMTRAJECTORY_MCWF_TRAJECTORY_TCC_INCLUDED

#include "MCWF_Trajectory.h"

#include "ParsMCWF_Trajectory.h"

#include "StateVector.h"

#include "EvolvedGSL.tcc"
#include "QuantumTrajectory.h"

#include "FormDouble.h"

#include <boost/range/numeric.hpp>
#include <boost/range/algorithm/find_if.hpp>


namespace quantumtrajectory {

//////////////////
//
// Implementations
//
//////////////////

template<int RANK, typename RandomEngine>
void MCWF_Trajectory<RANK,RandomEngine>::derivs(double t, const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  if (const auto ha=this->getHa()) {
    dpsidt=0;

    ha->addContribution(t,psi,dpsidt,this->getT0());
  }
}


template<int RANK, typename RandomEngine>
MCWF_Trajectory<RANK,RandomEngine>::MCWF_Trajectory(SV_Ptr psi,
                                                    typename structure::QuantumSystem<RANK>::Ptr sys,
                                                    const mcwf::Pars& p,
                                                    const StateVectorLow& scaleAbs)
  : QuantumTrajectory(sys,p.noise,
         psi->getArray(),
         bind(&MCWF_Trajectory::derivs,this,_1,_2,_3),
         initialTimeStep<RANK>(sys),
         scaleAbs,
         p,
         evolved::MakerGSL<StateVectorLow>(p.sf,p.nextDtTryCorrectionFactor)),
    psi_(psi),
    dpLimit_(p.dpLimit), overshootTolerance_(p.overshootTolerance),
    logger_(p.logLevel,this->template nAvr<structure::LA_Li>())
{
  QuantumTrajectory::checkDimension(psi);
  if (!getTime()) if(const auto li=this->getLi()) { // On startup, dpLimit should not be overshot, either.
    Rates rates(li->rates(0.,*psi_)); calculateSpecialRates(&rates,0.);
    manageTimeStep(rates,getEvolved().get(),false);
  }
}


template<int RANK, typename RandomEngine>
double MCWF_Trajectory<RANK,RandomEngine>::coherentTimeDevelopment(double Dt)
{
  if (this->getHa()) {
    getEvolved()->step(Dt);
  }
  else {
    double stepToDo=this->getLi() ? std::min(getDtTry(),Dt) : Dt; // Cf. tracker #3482771
    getEvolved()->update(getTime()+stepToDo,getDtTry());
  }
  // This defines three levels:
  // 1. System is Hamiltonian -> internal timestep control is used with dtTry possibly modified by Liouvillean needs (cf. manageTimeStep)
  // 2. System is not Hamiltonian, but it is Liouvillean -> dtTry is used, which is governed solely by Liouvillean in this case (cf. manageTimeStep)
  // 3. System is neither Hamiltonian nor Liouvillean (might be of Exact) -> as step of Dt is taken
  
  double t=getTime();

  if (const auto ex=this->getEx()) {
    ex->actWithU(getTime(),psi_->getArray(),this->getT0());
    QuantumTrajectory::setT0(t);
  }

  logger_.processNorm(psi_->renorm());

  return t;
}


template<int RANK, typename RandomEngine>
auto MCWF_Trajectory<RANK,RandomEngine>::calculateSpecialRates(Rates* rates, double t) const -> const IndexSVL_tuples
{
  IndexSVL_tuples res;
  for (int i=0; i<rates->size(); i++)
    if ((*rates)(i)<0) {
      StateVector psiTemp(*psi_);
      this->actWithJ(t,psiTemp.getArray(),i);
      res.push_back(IndexSVL_tuple(i,psiTemp.getArray()));
      (*rates)(i)=mathutils::sqr(psiTemp.renorm());
    } // psiTemp disappears here, but its storage does not, because the ownership is taken over by the StateVectorLow in the tuple
  return res;
}


template<int RANK, typename RandomEngine>
bool MCWF_Trajectory<RANK,RandomEngine>::manageTimeStep(const Rates& rates, evolved::TimeStepBookkeeper* evolvedCache, bool logControl)
{
  const double totalRate=boost::accumulate(rates,0.);
  const double dtDid=this->getDtDid(), dtTry=getDtTry();

  const double liouvilleanSuggestedDtTry=dpLimit_/totalRate;

  // Assumption: overshootTolerance_>=1 (equality is the limiting case of no tolerance)
  if (totalRate*dtDid>overshootTolerance_*dpLimit_) {
    evolvedCache->setDtTry(liouvilleanSuggestedDtTry);
    (*getEvolved())=*evolvedCache;
    logger_.stepBack(this->getLogStreamDuringRun(),totalRate*dtDid,dtDid,getDtTry(),getTime(),logControl);
    return true; // Step-back required.
  }
  else if ( totalRate*dtTry>dpLimit_) {
    logger_.overshot(this->getLogStreamDuringRun(),totalRate*dtTry,dtTry,liouvilleanSuggestedDtTry,logControl);
  }

  // dtTry-adjustment for next step:
  getEvolved()->setDtTry(this->getHa() ? std::min(getDtTry(),liouvilleanSuggestedDtTry) : liouvilleanSuggestedDtTry);

  return false; // Step-back not required.
}


template<int RANK, typename RandomEngine>
void MCWF_Trajectory<RANK,RandomEngine>::performJump(const Rates& rates, const IndexSVL_tuples& specialRates, double t)
{
  double random=std::uniform_real_distribution()(this->getRandomEngine())/this->getDtDid();

  int lindbladNo=0; // TODO: this could be expressed with an iterator into rates
  for (; random>0 && lindbladNo!=rates.size(); random-=rates(lindbladNo++))
    ;

  if(random<0) { // Jump corresponding to Lindblad no. lindbladNo-1 occurs
    --lindbladNo; auto i=boost::find_if(specialRates,[=](const IndexSVL_tuple& j){return lindbladNo==std::get<0>(j);});
    if (i!=specialRates.end())
      // special jump
      *psi_=std::get<1>(*i); // RHS already normalized above
    else {
      // normal  jump
      this->actWithJ(t,psi_->getArray(),lindbladNo);
      double normFactor=sqrt(rates(lindbladNo));
      if (!boost::math::isfinite(normFactor)) throw std::runtime_error("Infinite detected in MCWF_Trajectory::performJump");
      *psi_/=normFactor;
    }

    logger_.jumpOccured(this->getLogStreamDuringRun(),t,lindbladNo);
  }
}


template<int RANK, typename RandomEngine>
void MCWF_Trajectory<RANK,RandomEngine>::step_v(double Dt)
{
  const StateVector psiCache(*psi_);
  evolved::TimeStepBookkeeper evolvedCache(*getEvolved()); // This cannot be const since dtTry might change.

  double t=coherentTimeDevelopment(Dt);

  if (const auto li=this->getLi()) {

    Rates rates(li->rates(t,*psi_));
    IndexSVL_tuples specialRates=calculateSpecialRates(&rates,t);

    while (manageTimeStep(rates,&evolvedCache)) {
      *psi_=psiCache;
      t=coherentTimeDevelopment(Dt); // the next try
      rates=li->rates(t,*psi_);
      specialRates=calculateSpecialRates(&rates,t);
    }

    // Jump
    performJump(rates,specialRates,t);

  }

  logger_.step();

}


template<int RANK, typename RandomEngine>
std::ostream& MCWF_Trajectory<RANK,RandomEngine>::streamParameters_v(std::ostream& os) const
{
  using namespace std;
  
  this->streamCharacteristics(this->getQS()->streamParameters(Base::streamParameters_v(os)<<"MCWF Trajectory Parameters: dpLimit="<<dpLimit_<<" (overshoot tolerance factor)="<<overshootTolerance_<<endl<<endl))<<endl;

  if (const auto li=this->getLi()) {
    os<<"Decay channels:\n";
    {
      size_t i=0;
      li->streamKey(os,i);
    }
    os<<"Alternative Lindblads: ";
    {
      const Rates rates(li->rates(0,*psi_));
      int n=0;
      for (int i=0; i<rates.size(); ++i) if (rates(i)<0) {os<<i<<' '; ++n;}
      if (!n) os<<"none";
    }
    os<<endl;
  }

  return os;
  
}


template<int RANK, typename RandomEngine>
std::ostream& MCWF_Trajectory<RANK,RandomEngine>::streamKey_v(std::ostream& os, size_t& i) const
{
  return this->template streamKey<structure::LA_Av>(os,i);
}


} // quantumtrajectory


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_MCWF_TRAJECTORY_TCC_INCLUDED
