#include "ParticleTwoModes.h"

#include "Algorithm.h"
#include "BlitzArraySliceIterator.h"
#include "Range.h"

#include<boost/assign/list_of.hpp>
#include<boost/assign/list_inserter.hpp>


using namespace std;
using namespace boost::assign;

using namespace mathutils;

using particle::mfNKX;
using particle::ModeFunction;


using namespace mode;


namespace particletwomodes {


const quantumoperator::Tridiagonal<3> 
helper0(const ModeBase* mode0, const ModeBase* mode1, const ParticleBase* part, const ModeFunction& mf0)
{
  return aop(mode0)*aop(mode1).dagger()*mfNKX(part->getDimension(),mf0);
}


const quantumoperator::Frequencies<3> 
helper1(const ModeBase* mode0, const ModeBase* mode1, const ParticleBase* part, const ModeFunction& mf0)
{
  return freqs(mode0)*freqs(mode1)*freqs(part,mf0.get<1>());
}



const dcomp factor(double uNot0, double uNot1, double phi)
{
  if (uNot0*uNot1<0) throw UnotsSignDiscrepancy(); 
  return sign(uNot0)*sqrt(uNot0*uNot1)*exp(DCOMP_I*phi)/DCOMP_I;
}


Base::Base(const ModeBase* mode0, const ModeBase* mode1, const ParticleBase* part, 
	   double uNot0, double uNot1, 
	   const ModeFunction& mf0, const ModeFunction& mf1,
	   double phi)
  : structure::Interaction<3>(Frees(mode0,mode1,part),
			      tuple_list_of("Unot0",uNot0,sqrt(mode0->getDimension()))("Unot1",uNot1,sqrt(mode1->getDimension()))),
    firstH_(factor(uNot0,uNot1,phi)*helper0(mode0,mode1,part,mf0)), firstHT_(-firstH_.dagger()),
    firstF_(helper1(mode0,mode1,part,mf0)),
    secondH_(mfNKX(part->getDimension(),mf1).dagger()), secondHT_(secondH_.dagger()),
    secondF_(freqs(part,mf1.get<1>()))
{
  getParsStream()<<"# ParticleTwoModes\n# phi="<<phi<<endl;
}



void Base::addContribution(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double tIntPic0) const
{
  using namespace blitzplusplus;
  using basi::fullRange; using cpputils::for_each;

  typedef tmptools::Vector<2> V2;

  {
    double dt=t-tIntPic0;
    firstH_ .propagate( firstF_,dt); firstHT_. propagate( firstF_,dt);
    secondH_.propagate(secondF_,dt); secondHT_.propagate(secondF_,dt);
  }

  StateVectorLow dpsidtTemp(psi.shape()); // NEEDS_WORK check whether putting this into class scope saves time (only one dynamic allocation)
  {
    dpsidtTemp=0;
    apply(psi,dpsidtTemp,firstH_);
    for_each(fullRange(dpsidtTemp,V2()),basi::begin(dpsidt,V2()),bind(quantumoperator::apply<1>,_1,_2,secondH_));
  }
  {
    dpsidtTemp=0;
    apply(psi,dpsidtTemp,firstHT_);
    for_each(fullRange(dpsidtTemp,V2()),basi::begin(dpsidt,V2()),bind(quantumoperator::apply<1>,_1,_2,secondHT_));
  }

}


} // particletwomodes
