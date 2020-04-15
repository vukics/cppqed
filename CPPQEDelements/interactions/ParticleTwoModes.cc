// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ParticleTwoModes_.h"

#include "BlitzArraySliceIterator.tcc"
#include "MathExtensions.h"
#include "Tridiagonal.tcc"
#include "TridiagonalHamiltonian.tcc"


using namespace std;

using namespace mathutils;

using particle::mfNKX;
using particle::mfComposition;
using quantumoperator::identity;

using namespace mode;

namespace {


const dcomp factor(double uNot0, double uNot1, double phi)
{
  if (uNot0*uNot1<0) throw particletwomodes::UnotsSignDiscrepancy(); 
  return sign(uNot0)*sqrt(uNot0*uNot1)*exp(DCOMP_I*phi);
}


}


ParticleTwoModes::ParticleTwoModes(mode::Ptr mode0, mode::Ptr mode1, particle::Ptr part, 
                                   double uNot0, double uNot1, 
                                   const ModeFunction& mf0, const ModeFunction& mf1,
                                   double phi)
  : structure::Interaction<3>({mode0,mode1,part},
                              {RF{"Unot0",uNot0,sqrt(mode0->getDimension())},RF{"Unot1",uNot1,sqrt(mode1->getDimension())}}),
                              firstH_(factor(uNot0,uNot1,phi)*aop(mode0)*aop(mode1).dagger()*mfNKX(part,mf0)/DCOMP_I), firstHT_(-firstH_.dagger()),
    secondH_(mfNKX(part,mf1).dagger()), secondHT_(secondH_.dagger()),
    isSpecialH_(abs(get<1>(mf0))==abs(get<1>(mf1))),
    specialH_(isSpecialH_ ? Tridiagonals{quantumoperator::tridiagPlusHC_overI( factor(uNot0,uNot1,phi)*aop(mode0).dagger()*aop(mode1)*mfComposition(part,mf0,mf1) )}: Tridiagonals{})
{
  getParsStream()<<"ParticleTwoModes\nphi="<<phi<<endl;
}



void ParticleTwoModes::addContribution_v(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double t0) const
{
  if (isSpecialH_){
    specialH_.addContribution(t,psi,dpsidt,t0);
    return;
  }

  using namespace blitzplusplus;
  using basi::fullRange;

  typedef tmptools::Vector<2> V2;

  {
    double dt=t-t0;
    firstH_ .propagate(dt); firstHT_. propagate(dt);
    secondH_.propagate(dt); secondHT_.propagate(dt);
  }

  StateVectorLow dpsidtTemp(psi.shape()); // NEEDS_WORK check whether putting this into class scope saves time (only one dynamic allocation)
  {
    dpsidtTemp=0;
    apply(psi,dpsidtTemp,firstH_);
    boost::for_each(fullRange<V2>(dpsidtTemp),fullRange<V2>(dpsidt),bind(quantumoperator::apply<1>,_1,_2,secondH_));
  }
  {
    dpsidtTemp=0;
    apply(psi,dpsidtTemp,firstHT_);
    boost::for_each(fullRange<V2>(dpsidtTemp),fullRange<V2>(dpsidt),bind(quantumoperator::apply<1>,_1,_2,secondHT_));
  }

}


