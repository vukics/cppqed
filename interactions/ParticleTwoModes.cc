#include "ParticleTwoModes_.h"

#include "impl/Algorithm.tcc"
#include "impl/BlitzArraySliceIterator.tcc"
#include "MathExtensions.h"
#include "Range.h"
#include "impl/Tridiagonal.tcc"

#include <boost/assign/list_of.hpp>
#include <boost/assign/list_inserter.hpp>


using namespace std;
using namespace boost::assign;

using namespace mathutils;

using particle::mfNKX;

using namespace mode;


namespace particletwomodes {


const quantumoperator::Tridiagonal<3> 
helper(mode::Ptr mode0, mode::Ptr mode1, particle::Ptr part, const ModeFunction& mf0)
{
  return aop(mode0)*aop(mode1).dagger()*mfNKX(part,mf0);
}


const dcomp factor(double uNot0, double uNot1, double phi)
{
  if (uNot0*uNot1<0) throw UnotsSignDiscrepancy(); 
  return sign(uNot0)*sqrt(uNot0*uNot1)*exp(DCOMP_I*phi)/DCOMP_I;
}


Base::Base(mode::Ptr mode0, mode::Ptr mode1, particle::Ptr part, 
	   double uNot0, double uNot1, 
	   const ModeFunction& mf0, const ModeFunction& mf1,
	   double phi)
  : structure::Interaction<3>(Frees(mode0,mode1,part),
			      tuple_list_of("Unot0",uNot0,sqrt(mode0->getDimension()))("Unot1",uNot1,sqrt(mode1->getDimension()))),
    firstH_(factor(uNot0,uNot1,phi)*helper(mode0,mode1,part,mf0)), firstHT_(-firstH_.dagger()),
    secondH_(mfNKX(part,mf1).dagger()), secondHT_(secondH_.dagger())
{
  getParsStream()<<"# ParticleTwoModes\n# phi="<<phi<<endl;
}



void Base::addContribution_v(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double tIntPic0) const
{
  using namespace blitzplusplus;
  using basi::fullRange;

  typedef tmptools::Vector<2> V2;

  {
    double dt=t-tIntPic0;
    firstH_ .propagate(dt); firstHT_. propagate(dt);
    secondH_.propagate(dt); secondHT_.propagate(dt);
  }

  StateVectorLow dpsidtTemp(psi.shape()); // NEEDS_WORK check whether putting this into class scope saves time (only one dynamic allocation)
  {
    dpsidtTemp=0;
    apply(psi,dpsidtTemp,firstH_);
    cpputils::for_each(fullRange<V2>(dpsidtTemp),basi::begin<V2>(dpsidt),bind(quantumoperator::apply<1>,_1,_2,secondH_));
  }
  {
    dpsidtTemp=0;
    apply(psi,dpsidtTemp,firstHT_);
    cpputils::for_each(fullRange<V2>(dpsidtTemp),basi::begin<V2>(dpsidt),bind(quantumoperator::apply<1>,_1,_2,secondHT_));
  }

}


} // particletwomodes
