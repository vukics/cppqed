#include "BinarySystem.h"

#include "Interaction.h"

#include "impl/LazyDensityOperator.tcc"

#include "impl/Algorithm.tcc"
#include "BlitzArraySliceIterator.h"
#include "Range.h"

#include "BlitzTiny.h"

#include <boost/make_shared.hpp>


#define DISPLAY_KEY(Class,Aux) void binary::Class::displayKey(std::ostream& os, size_t& i) const \
  {									\
    os<<"# Binary system\n";						\
    Aux##1::displayKey(os,i,free0_.get##Aux());				\
    Aux##1::displayKey(os,i,free1_.get##Aux());				\
    Aux##2::displayKey(os,i,   ia_.get##Aux());				\
  }									\


#define ADD_UP_N(Class,Aux,Func) size_t binary::Class::n##Func() const	\
  {									\
    return Aux##1::n##Func(free0_.get##Aux()) + Aux##1::n##Func(free1_.get##Aux()) + Aux##2::n##Func(ia_.get##Aux()); \
  }									\


#define CONCATENATE_ARRAYS(Class,Aux,ArrayName,func,NumberName) const binary::Class::ArrayName binary::Class::func(double t, const LazyDensityOperator& ldo) const \
  {									\
  using quantumdata::partialTrace;					\
  using boost::copy;							\
									\
  const ArrayName a0 (partialTrace<V0,ArrayName>(ldo,bind(&Aux##1::func,t,_1,free0_.get##Aux(),theStaticOne))), \
    a1 (partialTrace<V1,ArrayName>(ldo,bind(&Aux##1::func,t,_1,free1_.get##Aux(),theStaticOne))), \
    a01(Aux##2::func(t,ldo,ia_.get##Aux()));				\
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


double binary::Base::highestFrequency () const
{
  using std::max;
  return max(ia_.get()->highestFrequency(),max(free0_.get()->highestFrequency(),free1_.get()->highestFrequency()));
}


void binary::Base::displayParameters(std::ostream& os) const
{
  using namespace std;
  os<<"# Binary System\n# Dimensions: "<<getDimensions()<<". Total: "<<getTotalDimension()<<endl<<endl
    <<"# Subsystem Nr. 0\n";     free0_.get()->displayParameters(os);
  os<<"# Subsystem Nr. 1\n";     free1_.get()->displayParameters(os);
  os<<"# 0 - 1 - Interaction\n"; ia_.get()->displayParameters(os);
}


DISPLAY_KEY(Base,Av);

ADD_UP_N(Base,Av,Avr);

CONCATENATE_ARRAYS(Base,Av,Averages,average,Avr);


void binary::Base::process(Averages& averages) const
{
  using blitz::Range;

  const Av1::Ptr 
    av0 =free0_.getAv(),
    av1 =free1_.getAv();
  const Av2::Ptr 
    av01=   ia_.getAv();

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



void binary::Base::display(const Averages& averages, std::ostream& os, int precision) const
{
  using blitz::Range;

  const Av1::Ptr 
    av0 =free0_.getAv(),
    av1 =free1_.getAv();
  const Av2::Ptr
    av01=   ia_.getAv();

  ptrdiff_t l=-1, u;

  if ((u=l+Av1::nAvr(av0 ))>l) av0 ->display(averages(Range(l+1,u)),os,precision);

  if ((l=u+Av1::nAvr(av1 ))>u) av1 ->display(averages(Range(u+1,l)),os,precision);

  if ((u=l+Av2::nAvr(av01))>l) av01->display(averages(Range(l+1,u)),os,precision);

}



///////////
//       //
// Exact //
//       //
///////////



bool binary::Exact::isUnitary() const
{
  return Ex1::isUnitary(free0_.getEx()) && Ex1::isUnitary(free1_.getEx()) && Ex2::isUnitary(ia_.getEx());
}



void binary::Exact::actWithU(double dt, StateVectorLow& psi) const
{
  using namespace blitzplusplus::basi;

  if (const Ex1::Ptr ex1=free0_.getEx()) for_each(fullRange<V0>(psi),bind(&Ex1::actWithU,ex1,dt,_1));
  if (const Ex1::Ptr ex1=free1_.getEx()) for_each(fullRange<V1>(psi),bind(&Ex1::actWithU,ex1,dt,_1));

  Ex2::actWithU(dt,psi,ia_.getEx());

}


/////////////////
//             //
// Hamiltonian //
//             //
/////////////////


void binary::Hamiltonian::addContribution(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double tIntPic0) const
{
  using namespace blitzplusplus; using basi::fullRange;

  if (const Ha1::Ptr ha1=free0_.getHa()) for_each(fullRange<V0>(psi),basi::begin<V0>(dpsidt),bind(&Ha1::addContribution,ha1,t,_1,_2,tIntPic0));
  if (const Ha1::Ptr ha1=free1_.getHa()) for_each(fullRange<V1>(psi),basi::begin<V1>(dpsidt),bind(&Ha1::addContribution,ha1,t,_1,_2,tIntPic0));

  Ha2::addContribution(t,psi,dpsidt,tIntPic0,ia_.getHa());

}


/////////////////
//             //
// Liouvillean //
//             //
/////////////////


ADD_UP_N(Liouvillean,Li,Jumps)

CONCATENATE_ARRAYS(Liouvillean,Li,Probabilities,probabilities,Jumps);


void binary::Liouvillean::actWithJ(double t, StateVectorLow& psi, size_t i) const
{
  using namespace blitzplusplus::basi;

  const Li1::Ptr
    li0 =free0_.getLi(),
    li1 =free1_.getLi();
  const Li2::Ptr
    li01=   ia_.getLi();

  size_t n=Li1::nJumps(li0);
  if (li0 && i<n) {
    for_each(fullRange<V0>(psi),bind(&Li1::actWithJ,li0,t,_1,i));
    return;
  }

  i-=n;  
  if (li1 && i<(n=Li1::nJumps(li1))) {
    for_each(fullRange<V1>(psi),bind(&Li1::actWithJ,li1,t,_1,i));
    return;
  }

  i-=n;
  if (i<(n=Li2::nJumps(li01)))
    Li2::actWithJ(t,psi,i,li01);

}


DISPLAY_KEY(Liouvillean,Li)



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

  return SystemCharacteristics(/* NEED_TO_UNDERSTAND Why can't the compiler deduce RANK here for frees? */
			       qse<1>(free0) || qse<1>(free1) || qse<2>(ia),
			       qsh<1>(free0) || qsh<1>(free1) || qsh<2>(ia),
			       qsl<1>(free0) || qsl<1>(free1) || qsl<2>(ia));
}

}


#define DISPATCHER(EX,HA,LI) (all(querySystemCharacteristics(ia)==SystemCharacteristics(EX,HA,LI))) return boost::make_shared<BinarySystem<EX,HA,LI> >(ia)


const binary::SmartPtr binary::make(Interaction::Ptr ia)
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
