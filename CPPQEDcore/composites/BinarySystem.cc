// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "BinarySystem.h"

#include "DensityOperator.h"
#include "Interaction.h"
#include "StateVector.h"

#include "Algorithm.h"
#include "SliceIterator.tcc"

#include "BlitzTiny.h"

#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/range/algorithm_ext/for_each.hpp>


using composite::SubSystemFree;


using namespace structure;

using cpputils::sliceiterator::fullRange;


namespace {

#include "details_BinaryHelper.h"

}


//////////
//      //
// Base //
//      //
//////////


binary::Base::Base(Interaction::Ptr ia)
  : QuantumSystem<2>(Dimensions(ia->getFrees()[0]->getDimension(),ia->getFrees()[1]->getDimension())),
    free0_(ia->getFrees()[0]), free1_(ia->getFrees()[1]), ia_(ia)
{
} 


double binary::Base::highestFrequency_v() const
{
  using std::max;
  return max(ia_.get()->highestFrequency(),max(free0_.get()->highestFrequency(),free1_.get()->highestFrequency()));
}


std::ostream& binary::Base::displayParameters_v(std::ostream& os) const
{
  using namespace std;
  os<<"Binary System\nDimensions: "<<getDimensions()<<". Total: "<<getTotalDimension()
    <<"\n\nSubsystem Nr. 0\n";     free0_.get()->displayParameters(os);
  os<<    "Subsystem Nr. 1\n";     free1_.get()->displayParameters(os);
  os<<"0 - 1 - Interaction\n";
  return ia_.get()->displayParameters(os);
}



#define SUCCESSIVE_Ranges(f0,f1,ia) ptrdiff_t l=-1, u;  \
  if ((u=l+free0_.nAvr<LA_Av>())>l) {                          \
    PROCESS_Range( averages(blitz::Range(l+1,u)) , f0)  \
  }                                                     \
  if ((l=u+free1_.nAvr<LA_Av>())>u) {                          \
    PROCESS_Range( averages(blitz::Range(u+1,l)) , f1)  \
  }                                                     \
  if ((u=l+ia_.nAvr<LA_Av>())>l) {                             \
    PROCESS_Range( averages(blitz::Range(l+1,u)) , ia)  \
  }                                                     \


void binary::Base::process_v(Averages& averages) const
{
#define PROCESS_Range(av,ss) Averages temp(av); ss.process(temp);
  
  SUCCESSIVE_Ranges(free0_,free1_,ia_) ;
  
#undef  PROCESS_Range
  
}


std::ostream& binary::Base::display_v(const Averages& averages, std::ostream& os, int precision) const
{
  const Av1::Ptr 
    av0 =free0_.getAv(),
    av1 =free1_.getAv();
  const Av2::Ptr
    av01=   ia_.getAv();

#define PROCESS_Range(av,ss) ss->display(av,os,precision);
  
  SUCCESSIVE_Ranges(av0,av1,av01) ;
  
#undef  PROCESS_Range

  return os;

}

#undef  SUCCESSIVE_Ranges

////////////////////////////
//                        //
// Averaged - Liouvillean //
//                        //
////////////////////////////


namespace binary {

// These possibilities get instantiated through compilation, so that explicit instantiation is not necessary here.
  
template<LiouvilleanAveragedTag LA>
std::ostream& displayKey(std::ostream& os, size_t& i, const SSF& free0, const SSF& free1, const SSI& ia)
{
  os<<"Binary system\n";
  free0.displayKey<LA>(os,i);
  free1.displayKey<LA>(os,i);
  return ia.displayKey<LA>(os,i);
}
template std::ostream& displayKey<LA_Li>(std::ostream& os, size_t& i, const SSF& free0, const SSF& free1, const SSI& ia);
template std::ostream& displayKey<LA_Av>(std::ostream& os, size_t& i, const SSF& free0, const SSF& free1, const SSI& ia);


template<LiouvilleanAveragedTag LA>
size_t nAvr(const SSF& free0, const SSF& free1, const SSI& ia)
{
  return free0.nAvr<LA>() + free1.nAvr<LA>() + ia.nAvr<LA>();
}
template size_t nAvr< LA_Li >(const SSF& free0, const SSF& free1, const SSI& ia); // explicit instantiation
template size_t nAvr< LA_Av >(const SSF& free0, const SSF& free1, const SSI& ia); // explicit instantiation


template<LiouvilleanAveragedTag LA>
const Base::Averages average(double t, const Base::LazyDensityOperator& ldo, const SSF& free0, const SSF& free1, const SSI& ia, size_t numberAvr)
{
  typedef Base::Averages ArrayName;
  
  using boost::copy;

  const ArrayName
    a0 {quantumdata::partialTrace<V0>(ldo,[&](const auto& m){return free0.average<LA>(t,m);})/*bind(&SubSystemFree::average<LA>,free0,t,std::placeholders::_1)*/},
    a1 {quantumdata::partialTrace<V1>(ldo,[&](const auto& m){return free1.average<LA>(t,m);})},
    a01{ia.average<LA>(t,ldo)};

  ArrayName a(numberAvr);

  copy(a01,copy(a1,copy(a0,a.begin())));

  return a;
}


} // binary


///////////
//       //
// Exact //
//       //
///////////



bool binary::Exact::applicableInMaster_v() const
{
  return free0_.applicableInMaster() && free1_.applicableInMaster() && ia_.applicableInMaster();
}



void binary::Exact::actWithU_v(double t, StateVectorLow& psi, double t0) const
{
  if (const auto ex=free0_.getEx()) for(auto& psiS : fullRange<V0>(psi)) ex->actWithU(t,psiS,t0);
  if (const auto ex=free1_.getEx()) for(auto& psiS : fullRange<V1>(psi)) ex->actWithU(t,psiS,t0);

  ia_.actWithU(t,psi,t0);

}


/////////////////
//             //
// Hamiltonian //
//             //
/////////////////


void binary::Hamiltonian::addContribution_v(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double t0) const
{
  using namespace std::placeholders;
  
  if (const auto ha=free0_.getHa()) boost::range::for_each(fullRange<V0>(psi),fullRange<V0>(dpsidt),bind(&Ha1::addContribution,ha,t,_1,_2,t0));
  if (const auto ha=free1_.getHa()) boost::range::for_each(fullRange<V1>(psi),fullRange<V1>(dpsidt),bind(&Ha1::addContribution,ha,t,_1,_2,t0));

  ia_.addContribution(t,psi,dpsidt,t0);

}


/////////////////
//             //
// Liouvillean //
//             //
/////////////////


void binary::Liouvillean::actWithJ_v(double t, StateVectorLow& psi, size_t i) const
{
  const Li1::Ptr
    li0 =free0_.getLi(),
    li1 =free1_.getLi();

  size_t n=free0_.nAvr<LA_Li>();
  if (li0 && i<n) {
    for(auto& psiS : fullRange<V0>(psi)) li0->actWithJ(t,psiS,i);
    return;
  }

  i-=n;  
  if (li1 && i<(n=free1_.nAvr<LA_Li>())) {
    for(auto& psiS : fullRange<V1>(psi)) li1->actWithJ(t,psiS,i);
    return;
  }

  i-=n;
  if (i<ia_.nAvr<LA_Li>())
    ia_.actWithJ(t,psi,i);

}


void binary::Liouvillean::actWithSuperoperator_v(double t, const DensityOperatorLow& rho, DensityOperatorLow& drhodt, size_t i) const
{
  using namespace std::placeholders;

  typedef tmptools::Vector<0,2> V0;
  typedef tmptools::Vector<1,3> V1;

  const Li1::Ptr
    li0 =free0_.getLi(),
    li1 =free1_.getLi();

  size_t n=free0_.nAvr<LA_Li>();
  if (li0 && i<n) {
    boost::range::for_each(fullRange<V0>(rho),fullRange<V0>(drhodt),bind(&Li1::actWithSuperoperator,li0,t,_1,_2,i));
    return;
  }

  i-=n;  
  if (li1 && i<(n=free1_.nAvr<LA_Li>())) {
    boost::range::for_each(fullRange<V1>(rho),fullRange<V1>(drhodt),bind(&Li1::actWithSuperoperator,li1,t,_1,_2,i));
    return;
  }

  i-=n;
  if (i<ia_.nAvr<LA_Li>())
    ia_.actWithSuperoperator(t,rho,drhodt,i);

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
    free0=ia->getFrees()[0],
    free1=ia->getFrees()[1];

  return SystemCharacteristics{qse<1>(free0) || qse<1>(free1) || qse<2>(ia),
                               qsh<1>(free0) || qsh<1>(free1) || qsh<2>(ia),
                               qsl<1>(free0) || qsl<1>(free1) || qsl<2>(ia)};
}

}


#define DISPATCHER(EX,HA,LI) (querySystemCharacteristics(ia)==SystemCharacteristics{EX,HA,LI}) return std::make_shared<BinarySystem<EX,HA,LI> >(ia)


const binary::Ptr binary::make(Interaction::Ptr ia)
{
  if      DISPATCHER(true ,true ,true ) ;
  else if DISPATCHER(true ,true ,false) ;
  else if DISPATCHER(true ,false,true ) ;
  else if DISPATCHER(true ,false,false) ;
  else if DISPATCHER(false,true ,true ) ;
  else if DISPATCHER(false,true ,false) ;
  else if DISPATCHER(false,false,true ) ;
  else return std::make_shared<BinarySystem<false,false,false> >(ia);
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

