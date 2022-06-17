// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ParticleTwoModes_.h"

#include "MathExtensions.h"
#include "SliceIterator.tcc"
#include "Tridiagonal.tcc"
#include "TridiagonalHamiltonian.h"


using namespace std;

using namespace cppqedutils;

using particle::mfNKX;
using particle::mfComposition;
using quantumoperator::identity;


namespace {


const dcomp factor(double uNot0, double uNot1, double phi)
{
  if (uNot0*uNot1<0) throw std::domain_error("Unots sign discrepancy"); 
  return sign(uNot0)*sqrt(uNot0*uNot1)*exp(1i*phi);
}


}


ParticleTwoModes::ParticleTwoModes(mode::Ptr mode0, mode::Ptr mode1, particle::Ptr part, 
                                   double uNot0, double uNot1, 
                                   const ModeFunction& mf0, const ModeFunction& mf1,
                                   double phi)
  : structure::Interaction<3>({mode0,mode1,part},
                              {RF{"Unot0",uNot0,sqrt(mode0->getDimension())},RF{"Unot1",uNot1,sqrt(mode1->getDimension())}}),
                              firstH_(factor(uNot0,uNot1,phi)*mode::aop(mode0)*mode::aop(mode1).dagger()*mfNKX(part,mf0)/1i), firstHT_(-firstH_.dagger()),
    secondH_(mfNKX(part,mf1).dagger()), secondHT_(secondH_.dagger()),
    isSpecialH_(abs(get<1>(mf0))==abs(get<1>(mf1))),
    specialH_(isSpecialH_ ? Tridiagonals{quantumoperator::tridiagPlusHC_overI( factor(uNot0,uNot1,phi)*mode::aop(mode0).dagger()*mode::aop(mode1)*mfComposition(part,mf0,mf1) )}: Tridiagonals{})
{
  getParsStream()<<"ParticleTwoModes\nphi="<<phi<<endl;
}



void ParticleTwoModes::addContribution_v(double t, const quantumdata::StateVectorLow<3>& psi, quantumdata::StateVectorLow<3>& dpsidt, double t0) const
{
  if (isSpecialH_){
    specialH_.addContribution(t,psi,dpsidt,t0);
    return;
  }

  using cppqedutils::sliceiterator::fullRange;

  typedef tmptools::Vector<2> V2;

  {
    double dt=t-t0;
    firstH_ .propagate(dt); firstHT_. propagate(dt);
    secondH_.propagate(dt); secondHT_.propagate(dt);
  }

  quantumdata::StateVectorLow<3> dpsidtTemp(psi.shape()); // NEEDS_WORK check whether putting this into class scope saves time (only one dynamic allocation)
  {
    dpsidtTemp=0;
    apply(psi,dpsidtTemp,firstH_);
    boost::for_each(fullRange<V2>(dpsidtTemp),fullRange<V2>(dpsidt),[&](const auto& p1, auto& p2){secondH_ .apply(p1,p2);});
  }
  {
    dpsidtTemp=0;
    apply(psi,dpsidtTemp,firstHT_);
    boost::for_each(fullRange<V2>(dpsidtTemp),fullRange<V2>(dpsidt),[&](const auto& p1, auto& p2){secondHT_.apply(p1,p2);});
  }

}


