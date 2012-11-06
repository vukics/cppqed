#ifndef BINARY_SYSTEM_AVERAGED_LIOUVILLEAN_ITERATING

#include "BinarySystem.h"

#include "Interaction.h"

#include "impl/LazyDensityOperator.tcc"

#include "impl/Algorithm.tcc"
#include "BlitzArraySliceIterator.h"
#include "Range.h"

#include "BlitzTiny.h"

#include <boost/make_shared.hpp>


using composite::SubSystemFree;


using namespace structure;

namespace {

#include "../interactions/details/BinaryHelper.h"

}


using cpputils::for_each;
using boost   ::for_each;


//////////
//      //
// Base //
//      //
//////////


binary::Base::Base(Interaction::Ptr ia)
  : QuantumSystem<2>(Dimensions(ia->getFrees()(0)->getDimension(),ia->getFrees()(1)->getDimension())),
    free0_(ia->getFrees()(0)), free1_(ia->getFrees()(1)), ia_(ia)
{
} 


double binary::Base::highestFrequency_v() const
{
  using std::max;
  return max(ia_.get()->highestFrequency(),max(free0_.get()->highestFrequency(),free1_.get()->highestFrequency()));
}


void binary::Base::displayParameters_v(std::ostream& os) const
{
  using namespace std;
  os<<"# Binary System\n# Dimensions: "<<getDimensions()<<". Total: "<<getTotalDimension()<<endl<<endl
    <<"# Subsystem Nr. 0\n";     free0_.get()->displayParameters(os);
  os<<"# Subsystem Nr. 1\n";     free1_.get()->displayParameters(os);
  os<<"# 0 - 1 - Interaction\n"; ia_.get()->displayParameters(os);
}



#define SUCCESSIVE_Ranges(f0,f1,ia) ptrdiff_t l=-1, u;  \
  if ((u=l+free0_.nAvr())>l) {                          \
    PROCESS_Range( averages(blitz::Range(l+1,u)) , f0)  \
  }                                                     \
  if ((l=u+free1_.nAvr())>u) {                          \
    PROCESS_Range( averages(blitz::Range(u+1,l)) , f1)  \
  }                                                     \
  if ((u=l+ia_.nAvr())>l) {                             \
    PROCESS_Range( averages(blitz::Range(l+1,u)) , ia)  \
  }                                                     \


void binary::Base::process_v(Averages& averages) const
{
#define PROCESS_Range(av,ss) Averages temp(av); ss.process(temp);
  
  SUCCESSIVE_Ranges(free0_,free1_,ia_) ;
  
#undef  PROCESS_Range
  
}


void binary::Base::display_v(const Averages& averages, std::ostream& os, int precision) const
{
  const Av1::Ptr 
    av0 =free0_.getAv(),
    av1 =free1_.getAv();
  const Av2::Ptr
    av01=   ia_.getAv();

#define PROCESS_Range(av,ss) ss->display(av,os,precision);
  
  SUCCESSIVE_Ranges(av0,av1,av01) ;
  
#undef  PROCESS_Range

}

#undef  SUCCESSIVE_Ranges

////////////////////////////
//                        //
// Averaged - Liouvillean //
//                        //
////////////////////////////


#define BINARY_SYSTEM_AVERAGED_LIOUVILLEAN_ITERATING

#define Class      Base
#define keyName    displayAveragedKey
#define numberName nAvr
#define ArrayName  Averages
#define func       average

#include "BinarySystem.cc"


#define BINARY_SYSTEM_AVERAGED_LIOUVILLEAN_ITERATING

#define Class      Liouvillean
#define keyName    displayLiouvilleanKey
#define numberName nJumps
#define ArrayName  Probabilities
#define func       probabilities

#include "BinarySystem.cc"



///////////
//       //
// Exact //
//       //
///////////



bool binary::Exact::isUnitary_v() const
{
  return free0_.isUnitary() && free1_.isUnitary() && ia_.isUnitary();
}



void binary::Exact::actWithU_v(double dt, StateVectorLow& psi) const
{
  using namespace blitzplusplus::basi;

  if (const Ex1::Ptr ex=free0_.getEx()) for_each(fullRange<V0>(psi),bind(&Ex1::actWithU,ex,dt,_1));
  if (const Ex1::Ptr ex=free1_.getEx()) for_each(fullRange<V1>(psi),bind(&Ex1::actWithU,ex,dt,_1));

  ia_.actWithU(dt,psi);

}


/////////////////
//             //
// Hamiltonian //
//             //
/////////////////


void binary::Hamiltonian::addContribution_v(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double tIntPic0) const
{
  using namespace blitzplusplus; using basi::fullRange;

  if (const Ha1::Ptr ha=free0_.getHa()) for_each(fullRange<V0>(psi),basi::begin<V0>(dpsidt),bind(&Ha1::addContribution,ha,t,_1,_2,tIntPic0));
  if (const Ha1::Ptr ha=free1_.getHa()) for_each(fullRange<V1>(psi),basi::begin<V1>(dpsidt),bind(&Ha1::addContribution,ha,t,_1,_2,tIntPic0));

  ia_.addContribution(t,psi,dpsidt,tIntPic0);

}


/////////////////
//             //
// Liouvillean //
//             //
/////////////////


void binary::Liouvillean::actWithJ_v(double t, StateVectorLow& psi, size_t i) const
{
  using namespace blitzplusplus::basi;

  const Li1::Ptr
    li0 =free0_.getLi(),
    li1 =free1_.getLi();

  size_t n=free0_.nJumps();
  if (li0 && i<n) {
    for_each(fullRange<V0>(psi),bind(&Li1::actWithJ,li0,t,_1,i));
    return;
  }

  i-=n;  
  if (li1 && i<(n=free1_.nJumps())) {
    for_each(fullRange<V1>(psi),bind(&Li1::actWithJ,li1,t,_1,i));
    return;
  }

  i-=n;
  if (i<ia_.nJumps())
    ia_.actWithJ(t,psi,i);

}



//////////////////
//              //
// Constructors //
//              //
//////////////////


#define BASE_ctor(Class) Class##Base(getFree0(),getFree1(),getIA())


template<bool IS_EX, bool IS_HA, bool IS_LI>
BinarySystem<IS_EX,IS_HA,IS_LI>::BinarySystem(Interaction::Ptr ia) 
: binary::Base(ia),
  BASE_ctor(Exact),
  BASE_ctor(Hamiltonian),
  BASE_ctor(Liouvillean)
{
} 


#undef BASE_ctor


namespace {

using structure::SystemCharacteristics;

const SystemCharacteristics querySystemCharacteristics(binary::Interaction::Ptr ia)
{
  using namespace structure;
  
  const QuantumSystem<1>::Ptr
    free0=ia->getFrees()(0),
    free1=ia->getFrees()(1);

  return SystemCharacteristics(qse(free0) || qse(free1) || qse<2>(ia),
			       qsh(free0) || qsh(free1) || qsh<2>(ia),
			       qsl(free0) || qsl(free1) || qsl<2>(ia));
}

}


#define DISPATCHER(EX,HA,LI) (all(querySystemCharacteristics(ia)==SystemCharacteristics(EX,HA,LI))) return boost::make_shared<BinarySystem<EX,HA,LI> >(ia)


const binary::Ptr binary::doMake(Interaction::Ptr ia)
{
  if      DISPATCHER(true ,true ,true ) ;
  else if DISPATCHER(true ,true ,false) ;
  else if DISPATCHER(true ,false,true ) ;
  else if DISPATCHER(true ,false,false) ;
  else if DISPATCHER(false,true ,true ) ;
  else if DISPATCHER(false,true ,false) ;
  else if DISPATCHER(false,false,true ) ;
  else return boost::make_shared<BinarySystem<false,false,false> >(ia);
}


#undef DISPATCHER

template class BinarySystem<true ,true ,true >;
template class BinarySystem<true ,true ,false>;
template class BinarySystem<true ,false,true >;
template class BinarySystem<true ,false,false>;
template class BinarySystem<false,true ,true >;
template class BinarySystem<false,true ,false>;
template class BinarySystem<false,false,true >;
template class BinarySystem<false,false,false>;


#else // BINARY_SYSTEM_AVERAGED_LIOUVILLEAN_ITERATING


void binary::Class::displayKey_v(std::ostream& os, size_t& i) const
{
  os<<"# Binary system\n";
  free0_.keyName(os,i);
  free1_.keyName(os,i);
  ia_   .keyName(os,i);
}


size_t binary::Class::BOOST_PP_CAT(numberName,_v)() const
{
  return free0_.numberName() + free1_.numberName() + ia_.numberName();
}


const binary::Class::ArrayName binary::Class::BOOST_PP_CAT(func,_v)(double t, const LazyDensityOperator& ldo) const
{
  using quantumdata::partialTrace;
  using boost::copy;

  const ArrayName
    a0 (partialTrace<V0,ArrayName>(ldo,bind(&SubSystemFree::func,free0_,t,_1))),
    a1 (partialTrace<V1,ArrayName>(ldo,bind(&SubSystemFree::func,free1_,t,_1))),
    a01(ia_.func(t,ldo));

  ArrayName a(numberName());

  copy(a01,copy(a1,copy(a0,a.begin())));

  return a;
}

#undef Class
#undef keyName
#undef numberName
#undef ArrayName
#undef func

#undef BINARY_SYSTEM_AVERAGED_LIOUVILLEAN_ITERATING

#endif // BINARY_SYSTEM_AVERAGED_LIOUVILLEAN_ITERATING

