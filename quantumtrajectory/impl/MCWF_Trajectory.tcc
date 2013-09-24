// -*- C++ -*-
#ifndef   QUANTUMTRAJECTORY_IMPL_MCWF_TRAJECTORY_TCC_INCLUDED
#define   QUANTUMTRAJECTORY_IMPL_MCWF_TRAJECTORY_TCC_INCLUDED

#include "MCWF_Trajectory.h"

#include "ParsMCWF_Trajectory.h"

#include "StateVector.h"
#include "impl/StochasticTrajectory.tcc"
#include "Structure.h"

#include "impl/FormDouble.tcc"
#include "SmartPtr.h"

#include <boost/range/numeric.hpp>


namespace quantumtrajectory {

//////////////////
//
// Implementations
//
//////////////////

template<int RANK>
void MCWF_Trajectory<RANK>::derivs(double t, const StateVectorLow& psi, StateVectorLow& dpsidt) const
{
  if (const typename Hamiltonian::Ptr ha=qs_.getHa()) {
    dpsidt=0;

    ha->addContribution(t,psi,dpsidt,tIntPic0_);
    logger_.hamiltonianCalled();
  }
}


template<int RANK> template<typename SYS>
MCWF_Trajectory<RANK>::MCWF_Trajectory(
                                       StateVector& psi,
                                       const SYS& sys,
                                       const ParsMCWF& p,
                                       const StateVectorLow& scaleAbs
                                       )
  : Base(psi(),
         bind(&MCWF_Trajectory::derivs,this,_1,_2,_3),
         trajectory::initialTimeStep(cpputils::sharedPointerize(sys)->highestFrequency()),
         scaleAbs,
         p,
         evolved::MakerGSL<StateVectorLow>(p.sf,p.nextDtTryCorretionFactor),
         randomized::MakerGSL()),
    tIntPic0_(0),
    psi_(psi),
    qs_(cpputils::sharedPointerize(sys),Base::noise()),
    dpLimit_(p.dpLimit), overshootTolerance_(p.overshootTolerance),
    logger_(p.logLevel,qs_.getHa()!=0,qs_.template nAvr<structure::LA_Li>())
{
  if (psi!=*qs_.getQS()) throw DimensionalityMismatchException();
  if (!getTime()) if(const typename Liouvillean::Ptr li=qs_.getLi()) { // On startup, dpLimit should not be overshot, either.
    Rates rates(li->rates(0.,psi_)); calculateSpecialRates(&rates,0.);
    manageTimeStep(rates,getEvolved().get(),false);
  }
}


template<int RANK>
std::ostream& MCWF_Trajectory<RANK>::display_v(std::ostream& os, int precision) const
{
  using namespace std;

  return qs_.display(getTime(),psi_,os,precision);

}


template<int RANK>
double MCWF_Trajectory<RANK>::coherentTimeDevelopment(double Dt)
{
  if (qs_.getHa()) {
    getEvolved()->step(Dt);
    logger_.logFailedSteps(getEvolved()->nFailedSteps());
  }
  else {
    double stepToDo=qs_.getLi() ? std::min(getDtTry(),Dt) : Dt; // Cf. tracker #3482771
    getEvolved()->update(getTime()+stepToDo,getDtTry());
  }
  // This defines three levels:
  // 1. System is Hamiltonian -> internal timestep control is used with dtTry possibly modified by Liouvillean needs (cf. manageTimeStep)
  // 2. System is not Hamiltonian, but it is Liouvillean -> dtTry is used, which is governed solely by Liouvillean in this case (cf. manageTimeStep)
  // 3. System is neither Hamiltonian nor Liouvillean (might be of Exact) -> as step of Dt is taken
  
  double t=getTime();

  if (const typename Exact::Ptr ex=qs_.getEx()) {
    ex->actWithU(getTime(),psi_(),tIntPic0_);
    tIntPic0_=t;
  }

  logger_.processNorm(psi_.renorm());

  return t;
}


template<int RANK>
const typename MCWF_Trajectory<RANK>::IndexSVL_tuples
MCWF_Trajectory<RANK>::calculateSpecialRates(Rates* rates, double t) const
{
  IndexSVL_tuples res;
  for (int i=0; i<rates->size(); i++)
    if ((*rates)(i)<0) {
      StateVector psiTemp(psi_);
      qs_.actWithJ(t,psiTemp(),i);
      res.push_back(IndexSVL_tuple(i,psiTemp()));
      (*rates)(i)=mathutils::sqr(psiTemp.renorm());
    } // psiTemp disappears here, but its storage does not, because the ownership is taken over by the StateVectorLow in the tuple
  return res;
}


template<int RANK>
bool MCWF_Trajectory<RANK>::manageTimeStep(const Rates& rates, evolved::TimeStepBookkeeper* evolvedCache, bool logControl)
{
  const double totalRate=boost::accumulate(rates,0.);
  const double dtDid=getDtDid(), dtTry=getDtTry();

  // Assumption: overshootTolerance_>=1 (equality is the limiting case of no tolerance)
  if (totalRate*dtDid>overshootTolerance_*dpLimit_) {
    evolvedCache->setDtTry(dpLimit_/totalRate);
    (*getEvolved())=*evolvedCache;
    logger_.stepBack(totalRate*dtDid,dtDid,getDtTry(),getTime(),logControl);
    return true; // Step-back required.
  }
  else if (totalRate*dtTry>dpLimit_) {
    logger_.overshot(totalRate*dtTry,dtTry,getDtTry(),logControl);
  }

  { // dtTry-adjustment for next step:
    const double liouvilleanSuggestedDtTry=dpLimit_/totalRate;
    getEvolved()->setDtTry(qs_.getHa() ? std::min(getDtTry(),liouvilleanSuggestedDtTry) : liouvilleanSuggestedDtTry);
  }

  return false; // Step-back not required.
}


template<int RANK>
void MCWF_Trajectory<RANK>::performJump(const Rates& rates, const IndexSVL_tuples& specialRates, double t)
{
  double random=(*getRandomized())()/getDtDid();

  int jumpNo=0;
  for (; random>0 && jumpNo!=rates.size(); random-=rates(jumpNo++))
    ;

  if(random<0) { // Jump No. jumpNo-1 occurs
    struct helper
    {
      static bool p(int i, IndexSVL_tuple j) {return i==j.template get<0>();} // NEEDS_WORK how to express this with lambda?
    };

    auto i=find_if(specialRates,bind(&helper::p,--jumpNo,_1)); // See whether it's a special jump
    if (i!=specialRates.end())
      // special jump
      psi_()=i->template get<1>(); // RHS already normalized above
    else {
      // normal  jump
      qs_.actWithJ(t,psi_(),jumpNo);
      double normFactor=sqrt(rates(jumpNo));
      if (!boost::math::isfinite(normFactor)) throw structure::InfiniteDetectedException();
      psi_()/=normFactor;
    }

    logger_.jumpOccured(t,jumpNo);
  }
}


template<int RANK>
void MCWF_Trajectory<RANK>::step_v(double Dt)
{
  const StateVectorLow psiCache(psi_().copy());
  evolved::TimeStepBookkeeper evolvedCache(*getEvolved()); // This cannot be const since dtTry might change.

  double t=coherentTimeDevelopment(Dt);

  if (const typename Liouvillean::Ptr li=qs_.getLi()) {

    Rates rates(li->rates(t,psi_));
    IndexSVL_tuples specialRates=calculateSpecialRates(&rates,t);

    while (manageTimeStep(rates,&evolvedCache)) {
      psi_()=psiCache;
      t=coherentTimeDevelopment(Dt); // the next try
      rates=li->rates(t,psi_);
      specialRates=calculateSpecialRates(&rates,t);
    }

    // Jump
    performJump(rates,specialRates,t);

  }

  logger_.step();

}


template<int RANK>
std::ostream& MCWF_Trajectory<RANK>::displayParameters_v(std::ostream& os) const
{
  using namespace std;
  
  qs_.displayCharacteristics(qs_.getQS()->displayParameters(Base::displayParameters_v(os)<<"# MCWF Trajectory Parameters: dpLimit="<<dpLimit_<<" (overshoot tolerance factor)="<<overshootTolerance_<<endl<<endl))<<endl;

  if (const typename Liouvillean::Ptr li=qs_.getLi()) {
    os<<"# Decay channels:\n";
    {
      size_t i=0;
      li->displayKey(os,i);
    }
    os<<"# Alternative jumps: ";
    {
      const Rates rates(li->rates(0,psi_));
      int n=0;
      for (int i=0; i<rates.size(); i++) if (rates(i)<0) {os<<i<<' '; n++;}
      if (!n) os<<"none";
    }
    os<<endl;
  }

  return os;
  
}


template<int RANK>
std::ostream& MCWF_Trajectory<RANK>::displayKey_v(std::ostream& os, size_t& i) const
{
  return qs_.template displayKey<structure::LA_Av>(os,i);
}


} // quantumtrajectory


#endif // QUANTUMTRAJECTORY_IMPL_MCWF_TRAJECTORY_TCC_INCLUDED
