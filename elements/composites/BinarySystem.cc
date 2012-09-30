#include "BinarySystem.h"

#include "Interaction.h"

#include "impl/LazyDensityOperator.tcc"

#include "impl/Algorithm.tcc"
#include "BlitzArraySliceIterator.h"
#include "Range.h"

#include "BlitzTiny.h"

#include <boost/make_shared.hpp>


using composite::SubSystemFree;


#define DISPLAY_KEY(Class,Aux) void binary::Class::displayKey_v(std::ostream& os, size_t& i) const \
  {									\
    os<<"# Binary system\n";						\
    free0_.display##Aux##Key(os,i);					\
    free1_.display##Aux##Key(os,i);					\
    ia_   .display##Aux##Key(os,i);					\
  }									\


#define ADD_UP_N(Class,Func) size_t binary::Class::n##Func##_v() const	\
  {									\
    return free0_.n##Func() + free1_.n##Func() + ia_.n##Func(); \
  }									\


#define CONCATENATE_ARRAYS(Class,ArrayName,func,NumberName) const binary::Class::ArrayName binary::Class::func##_v(double t, const LazyDensityOperator& ldo) const \
  {									\
  using quantumdata::partialTrace;					\
  using boost::copy;							\
									\
  const ArrayName							\
    a0 (partialTrace<V0,ArrayName>(ldo,bind(&SubSystemFree::func,free0_,t,_1))), \
    a1 (partialTrace<V1,ArrayName>(ldo,bind(&SubSystemFree::func,free1_,t,_1))), \
    a01(ia_.func(t,ldo));						\
									\
  ArrayName a(n##NumberName());						\
									\
  copy(a01,copy(a1,copy(a0,a.begin())));				\
									\
  return a;								\
  }									\



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


DISPLAY_KEY(Base,Averaged);

ADD_UP_N(Base,Avr);

CONCATENATE_ARRAYS(Base,Averages,average,Avr);


void binary::Base::process_v(Averages& averages) const
{
  using blitz::Range;

  ptrdiff_t l=-1, u;

  if ((u=l+free0_.nAvr())>l) {
    Averages temp(averages(Range(l+1,u)));
    free0_.process(temp);
  }
  if ((l=u+free1_.nAvr())>u) {
    Averages temp(averages(Range(u+1,l)));
    free1_.process(temp);
  }
  if ((u=l+ia_.nAvr())>l) {
    Averages temp(averages(Range(l+1,u)));
    ia_.process(temp);
  }

}



void binary::Base::display_v(const Averages& averages, std::ostream& os, int precision) const
{
  using blitz::Range;

  const Av1::Ptr 
    av0 =free0_.getAv(),
    av1 =free1_.getAv();
  const Av2::Ptr
    av01=   ia_.getAv();

  ptrdiff_t l=-1, u;

  if ((u=l+free0_.nAvr())>l) av0 ->display(averages(Range(l+1,u)),os,precision);

  if ((l=u+free1_.nAvr())>u) av1 ->display(averages(Range(u+1,l)),os,precision);

  if ((u=l+   ia_.nAvr())>l) av01->display(averages(Range(l+1,u)),os,precision);

}



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


ADD_UP_N(Liouvillean,Jumps)

CONCATENATE_ARRAYS(Liouvillean,Probabilities,probabilities,Jumps);


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


DISPLAY_KEY(Liouvillean,Liouvillean)



//////////////////
//              //
// Constructors //
//              //
//////////////////


#define BASE_CTOR(Class) Class##Base(getFree0(),getFree1(),getIA())


template<bool IS_EX, bool IS_HA, bool IS_LI>
BinarySystem<IS_EX,IS_HA,IS_LI>::BinarySystem(Interaction::Ptr ia) 
: binary::Base(ia),
  BASE_CTOR(Exact),
  BASE_CTOR(Hamiltonian),
  BASE_CTOR(Liouvillean)
{
} 


#undef BASE_CTOR
#undef CONCATENATE_ARRAYS
#undef ADD_UP_N
#undef DISPLAY_KEY


namespace {

typedef blitz::TinyVector<bool,3> SystemCharacteristics;

const SystemCharacteristics querySystemCharacteristics(binary::Interaction::Ptr ia)
{
  using namespace structure;
  
  QuantumSystem<1>::Ptr
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
