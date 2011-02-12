#include "BinarySystem.h"

#include "Interaction.h"

#include "LazyDensityOperatorSliceIterator.h"

#include "Algorithm.h"
#include "BlitzArraySliceIterator.h"
#include "Range.h"


namespace {

#include "../interactions/details/BinaryHelper.h"

}


using cpputils::for_each;
using boost   ::for_each;


BinarySystem::BinarySystem(const Interaction& ia)
  : QS2(Dimensions(ia.getFrees()(0)->getDimension(),ia.getFrees()(1)->getDimension())),
    free0_(ia.getFrees()(0)), free1_(ia.getFrees()(1)), ia_(&ia)
{
} 


// NEEDS_WORK a lot of boilerplate should be eliminated!!!


double BinarySystem::highestFrequency () const
{
  using std::max;
  return max(ia_.get()->highestFrequency(),max(free0_.get()->highestFrequency(),free1_.get()->highestFrequency()));
}



void BinarySystem::displayParameters(std::ostream& os) const
{
  using namespace std;
  os<<"# Binary System\n# Dimensions: "<<getDimensions()<<". Total: "<<getTotalDimension()<<endl<<endl
    <<"# Subsystem Nr. 0\n";     free0_.get()->displayParameters(os);
  os<<"# Subsystem Nr. 1\n";     free1_.get()->displayParameters(os);
  os<<"# 0 - 1 - Interaction\n"; ia_.get()->displayParameters(os);
}



bool BinarySystem::isUnitary() const
{
  return Ex1::isUnitary(free0_.getEx()) && Ex1::isUnitary(free1_.getEx()) && Ex2::isUnitary(ia_.getEx());
}



void BinarySystem::actWithU(double dt, StateVectorLow& psi) const
{
  using namespace blitzplusplus::basi;
  if (const Ex1* ex1=free0_.getEx()) for_each(fullRange(psi,v0),bind(&Ex1::actWithU,ex1,dt,_1));
  if (const Ex1* ex1=free1_.getEx()) for_each(fullRange(psi,v1),bind(&Ex1::actWithU,ex1,dt,_1));

  Ex2::actWithU(dt,psi,ia_.getEx());
}



void BinarySystem::addContribution(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double tIntPic0) const
{
  using namespace blitzplusplus; using basi::fullRange;
  if (const Ha1* ha1=free0_.getHa()) for_each(fullRange(psi,v0),basi::begin(dpsidt,v0),bind(&Ha1::addContribution,ha1,t,_1,_2,tIntPic0));
  if (const Ha1* ha1=free1_.getHa()) for_each(fullRange(psi,v1),basi::begin(dpsidt,v1),bind(&Ha1::addContribution,ha1,t,_1,_2,tIntPic0));

  Ha2::addContribution(t,psi,dpsidt,tIntPic0,ia_.getHa());

}



size_t BinarySystem::nJumps() const
{
  return Li1::nJumps(free0_.getLi()) + Li1::nJumps(free1_.getLi()) + Li2::nJumps(ia_.getLi());
}



const BinarySystem::Probabilities BinarySystem::probabilities(const LazyDensityOperator& ldo) const
{
  using quantumdata::partialTrace;
  using boost::copy;
  
  const Probabilities 
    p0 (partialTrace(ldo,bind(&Li1::probabilities,_1,free0_.getLi(),theStaticOne),v0,defaultArray)),
    p1 (partialTrace(ldo,bind(&Li1::probabilities,_1,free1_.getLi(),theStaticOne),v1,defaultArray)),
    p01(Li2::probabilities(ldo,ia_.getLi()));

  Probabilities p(nJumps());

  copy(p01,copy(p1,copy(p0,p.begin())));

  return p;

}



void BinarySystem::actWithJ(StateVectorLow& psi, size_t i) const
{
  using namespace blitzplusplus::basi;

  const Li1 
    * li0 =free0_.getLi(),
    * li1 =free1_.getLi();
  const Li2 
    * li01=   ia_.getLi();

  size_t n=Li1::nJumps(li0);
  if (li0 && i<n) {
    for_each(fullRange(psi,v0),bind(&Li1::actWithJ,li0,_1,i));
    return;
  }

  i-=n;  
  if (li1 && i<(n=Li1::nJumps(li1))) {
    for_each(fullRange(psi,v1),bind(&Li1::actWithJ,li1,_1,i));
    return;
  }

  i-=n;
  if (i<(n=Li2::nJumps(li01)))
    Li2::actWithJ(psi,i,li01);

}



void BinarySystem::displayKey(std::ostream& os, size_t& i) const
{
  os<<"# Binary system\n";
  Av1::displayKey(os,i,free0_.getAv());
  Av1::displayKey(os,i,free1_.getAv());
  Av2::displayKey(os,i,   ia_.getAv());
}



size_t BinarySystem::nAvr() const
{
  return Av1::nAvr(free0_.getAv()) + Av1::nAvr(free1_.getAv()) + Av2::nAvr(ia_.getAv());
}



const BinarySystem::Averages BinarySystem::average(const LazyDensityOperator& ldo) const
{
  using quantumdata::partialTrace;
  using boost::copy;

  const Averages 
    a0 (partialTrace(ldo,bind(&Av1::average,_1,free0_.getAv(),theStaticOne),v0,defaultArray)),
    a1 (partialTrace(ldo,bind(&Av1::average,_1,free1_.getAv(),theStaticOne),v1,defaultArray)),
    a01(Av2::average(ldo,ia_.getAv()));

  Averages a(nAvr());

  copy(a01,copy(a1,copy(a0,a.begin())));

  return a;  
}



void BinarySystem::process(Averages& averages) const
{
  using blitz::Range;

  const Av1 
    * av0 =free0_.getAv(),
    * av1 =free1_.getAv();
  const Av2 
    * av01=   ia_.getAv();

  ptrdiff_t l=-1, u;

  if ((u=l+Av1::nAvr(av0 ))>l) {
    Averages temp(averages(Range(l+1,u)));
    Av1::process(temp,av0 );
  }
  if ((l=u+Av1::nAvr(av1 ))>u) {
    Averages temp(averages(Range(u+1,l)));
    Av1::process(temp,av1 );
  }
  if ((u=l+Av2::nAvr(av01))>l) {
    Averages temp(averages(Range(l+1,u)));
    Av2::process(temp,av01);
  }

}



void BinarySystem::display(const Averages& averages, std::ostream& os, int precision) const
{
  using blitz::Range;

  const Av1 
    * av0 =free0_.getAv(),
    * av1 =free1_.getAv();
  const Av2 
    * av01=   ia_.getAv();

  ptrdiff_t l=-1, u;

  if ((u=l+Av1::nAvr(av0 ))>l) av0 ->display(averages(Range(l+1,u)),os,precision);

  if ((l=u+Av1::nAvr(av1 ))>u) av1 ->display(averages(Range(u+1,l)),os,precision);

  if ((u=l+Av2::nAvr(av01))>l) av01->display(averages(Range(l+1,u)),os,precision);

}
