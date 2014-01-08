#include "BinarySystem.h"

#include "Interaction.h"

#include "LazyDensityOperator.tcc"

#include "Algorithm.h"
#include "BlitzArraySliceIterator.h"

#include "BlitzTiny.h"

#include <boost/range/algorithm/for_each.hpp>
#include <boost/range/algorithm_ext/for_each.hpp>

#include <boost/make_shared.hpp>


using composite::SubSystemFree;


using namespace structure;

namespace {

#include "details_BinaryHelper.h"

}


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


std::ostream& binary::Base::displayParameters_v(std::ostream& os) const
{
  using namespace std;
  os<<"# Binary System\n# Dimensions: "<<getDimensions()<<". Total: "<<getTotalDimension()<<endl<<endl
    <<"# Subsystem Nr. 0\n";     free0_.get()->displayParameters(os);
  os<<"# Subsystem Nr. 1\n";     free1_.get()->displayParameters(os);
  os<<"# 0 - 1 - Interaction\n";
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
  os<<"# Binary system\n";
  free0.displayKey<LA>(os,i);
  free1.displayKey<LA>(os,i);
  return ia.displayKey<LA>(os,i);
}


template<LiouvilleanAveragedTag LA>
size_t nAvr(const SSF& free0, const SSF& free1, const SSI& ia)
{
  return free0.nAvr<LA>() + free1.nAvr<LA>() + ia.nAvr<LA>();
}


template<LiouvilleanAveragedTag LA>
const Base::Averages average(double t, const Base::LazyDensityOperator& ldo, const SSF& free0, const SSF& free1, const SSI& ia, size_t numberAvr)
{
  typedef Base::Averages ArrayName;
  
  using quantumdata::partialTrace;
  using boost::copy;

  const ArrayName
    a0 (partialTrace<V0,ArrayName>(ldo,bind(&SubSystemFree::average<LA>,free0,t,_1))),
    a1 (partialTrace<V1,ArrayName>(ldo,bind(&SubSystemFree::average<LA>,free1,t,_1))),
    a01(ia.average<LA>(t,ldo));

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



bool binary::Exact::isUnitary_v() const
{
  return free0_.isUnitary() && free1_.isUnitary() && ia_.isUnitary();
}



void binary::Exact::actWithU_v(double t, StateVectorLow& psi, double t0) const
{
  using namespace blitzplusplus::basi;

  if (const auto ex=free0_.getEx()) for_each(fullRange<V0>(psi),bind(&Ex1::actWithU,ex,t,_1,t0));
  if (const auto ex=free1_.getEx()) for_each(fullRange<V1>(psi),bind(&Ex1::actWithU,ex,t,_1,t0));

  ia_.actWithU(t,psi,t0);

}


/////////////////
//             //
// Hamiltonian //
//             //
/////////////////


void binary::Hamiltonian::addContribution_v(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double tIntPic0) const
{
  using namespace blitzplusplus; using basi::fullRange;

  if (const auto ha=free0_.getHa()) for_each(fullRange<V0>(psi),fullRange<V0>(dpsidt),bind(&Ha1::addContribution,ha,t,_1,_2,tIntPic0));
  if (const auto ha=free1_.getHa()) for_each(fullRange<V1>(psi),fullRange<V1>(dpsidt),bind(&Ha1::addContribution,ha,t,_1,_2,tIntPic0));

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

  size_t n=free0_.nAvr<LA_Li>();
  if (li0 && i<n) {
    for_each(fullRange<V0>(psi),bind(&Li1::actWithJ,li0,t,_1,i));
    return;
  }

  i-=n;  
  if (li1 && i<(n=free1_.nAvr<LA_Li>())) {
    for_each(fullRange<V1>(psi),bind(&Li1::actWithJ,li1,t,_1,i));
    return;
  }

  i-=n;
  if (i<ia_.nAvr<LA_Li>())
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

