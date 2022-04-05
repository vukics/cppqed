// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "BinarySystem.h"

#include "DensityOperator.h"
#include "StateVector.h"

#include "Algorithm.h"
#include "SliceIterator.h"

#include "BlitzTiny.h"


using namespace structure;

using cppqedutils::sliceiterator::fullRange;


//////////
//      //
// Base //
//      //
//////////




double binary::Base::highestFrequency_v() const
{
  return std::max(ia_->highestFrequency(),std::max((*ia_)[0]->highestFrequency(),(*ia_)[1]->highestFrequency()));
}


std::ostream& binary::Base::streamParameters_v(std::ostream& os) const
{
  return ia_->streamParameters((*ia_)[1]->streamParameters(
    (*ia_)[0]->streamParameters(os<<"Binary System\nDimensions: "<<getDimensions()<<". Total: "<<getTotalDimension()<<"\n\nSubsystem Nr. 0\n")
    <<"Subsystem Nr. 1\n")<<"0 - 1 - Interaction\n");
}



void binary::Base::process_v(Averages& averages) const
{
  ptrdiff_t l=-1, u;
  
  if ( ( u = l + ::structure::nAvr<LA_Av>((*ia_)[0]) ) > l ) {Averages temp(averages(blitz::Range(l+1,u))); ::structure::process((*ia_)[0],temp);}
  if ( ( l = u + ::structure::nAvr<LA_Av>((*ia_)[1]) ) > u ) {Averages temp(averages(blitz::Range(u+1,l))); ::structure::process((*ia_)[1],temp);}
  if ( ( u = l + ::structure::nAvr<LA_Av>(ia_) ) > l ) {Averages temp(averages(blitz::Range(l+1,u))); ::structure::process(ia_,temp);}
  
}


std::ostream& binary::Base::stream_v(const Averages& averages, std::ostream& os, int precision) const
{
  ptrdiff_t l=-1, u;
  
  if ( const auto av=::structure::castAv((*ia_)[0]); av && ( u = l + ::structure::nAvr<LA_Av>((*ia_)[0]) ) > l ) av->stream(averages(blitz::Range(l+1,u)),os,precision);
  if ( const auto av=::structure::castAv((*ia_)[1]); av && ( l = u + ::structure::nAvr<LA_Av>((*ia_)[1]) ) > u ) av->stream(averages(blitz::Range(u+1,l)),os,precision);
  if ( const auto av=::structure::castAv(ia_); av && ( u = l + ::structure::nAvr<LA_Av>(ia_) ) > l ) av->stream(averages(blitz::Range(l+1,u)),os,precision);

  return os;

}


///////////
//       //
// Exact //
//       //
///////////



bool binary::Exact::applicableInMaster_v() const
{
  return ::structure::applicableInMaster((*ia_)[0]) && ::structure::applicableInMaster((*ia_)[1]) && ::structure::applicableInMaster(ia_);
}



void binary::Exact::actWithU_v(double t, StateVectorLow& psi, double t0) const
{
  if (const auto ex=castEx((*ia_)[0])) for(auto& psiS : fullRange<V0>(psi)) ex->actWithU(t,psiS,t0);
  if (const auto ex=castEx((*ia_)[1])) for(auto& psiS : fullRange<V1>(psi)) ex->actWithU(t,psiS,t0);

  ::structure::actWithU(ia_,t,psi,t0);

}


/////////////////
//             //
// Hamiltonian //
//             //
/////////////////


void binary::Hamiltonian::addContribution_v(double t, const StateVectorLow& psi, StateVectorLow& dpsidt, double t0) const
{
  const auto lambda=[=](auto ha) {
    return [=](const auto& psiS, auto& dpsidtS) {
      ha->addContribution(t,psiS,dpsidtS,t0);
    };
  };
  
  if (const auto ha=castHa((*ia_)[0])) std::ranges::for_each(fullRange<V0>(psi),fullRange<V0>(dpsidt),lambda(ha));
  if (const auto ha=castHa((*ia_)[1])) std::ranges::for_each(fullRange<V1>(psi),fullRange<V1>(dpsidt),lambda(ha));

  ::structure::addContribution(ia_,t,psi,dpsidt,t0);

}


/////////////////
//             //
// Liouvillean //
//             //
/////////////////


void binary::Liouvillean::actWithJ_v(double t, StateVectorLow& psi, size_t i) const
{
  size_t n=::structure::nAvr<LA_Li>((*ia_)[0]);
  
  if ( const auto li0=castLi((*ia_)[0]); li0 && i<n) {
    for(auto& psiS : fullRange<V0>(psi)) li0->actWithJ(t,psiS,i);
    return;
  }

  i-=n;  
  if ( const auto li1=castLi((*ia_)[1]); li1 && i < ( n = ::structure::nAvr<LA_Li>((*ia_)[1])) ) {
    for(auto& psiS : fullRange<V1>(psi)) li1->actWithJ(t,psiS,i);
    return;
  }

  i-=n;
  if ( i < ::structure::nAvr<LA_Li>(ia_) )
    ::structure::actWithJ(ia_,t,psi,i);

}


void binary::Liouvillean::actWithSuperoperator_v(double t, const DensityOperatorLow& rho, DensityOperatorLow& drhodt, size_t i) const
{
  const auto lambda=[=,&i](auto li) {
    return [=,&i](const auto& rhoS, auto& drhodtS) {
      li->actWithSuperoperator(t,rhoS,drhodtS,i);
    };
  };

  using V0=tmptools::Vector<0,2>;
  using V1=tmptools::Vector<1,3>;

  size_t n=::structure::nAvr<LA_Li>((*ia_)[0]);
  
  if ( const auto li0=castLi((*ia_)[0]);  li0 && i<n ) {
    std::ranges::for_each(fullRange<V0>(rho),fullRange<V0>(drhodt),lambda(li0));
    return;
  }

  i-=n;
  if ( const auto li1=castLi((*ia_)[1]);  li1 && i < ( n = ::structure::nAvr<LA_Li>((*ia_)[1]) ) ) {
    std::ranges::for_each(fullRange<V1>(rho),fullRange<V1>(drhodt),lambda(li1));
    return;
  }

  i-=n;
  if ( i < ::structure::nAvr<LA_Li>(ia_) ) ::structure::actWithSuperoperator(ia_,t,rho,drhodt,i);

}


//////////////////
//              //
// Constructors //
//              //
//////////////////


#define BASE_ctor(Class) Class##Base(ia)


template<bool IS_EX, bool IS_HA, bool IS_LI>
BinarySystem<IS_EX,IS_HA,IS_LI>::BinarySystem(binary::InteractionPtr ia) 
: binary::Base(ia),
  BASE_ctor(Exact),
  BASE_ctor(Hamiltonian),
  BASE_ctor(Liouvillean)
{
} 


#undef BASE_ctor


namespace {

using structure::SystemCharacteristics;

const SystemCharacteristics querySystemCharacteristics(binary::InteractionPtr ia)
{
  using namespace structure;
  
  const auto free0=(*ia)[0], free1=(*ia)[1];

  return SystemCharacteristics{
    castEx(free0) || castEx(free1) || castEx(ia),
    castHa(free0) || castHa(free1) || castHa(ia),
    castLi(free0) || castLi(free1) || castLi(ia)};
}

}


#define DISPATCHER(EX,HA,LI) (all(querySystemCharacteristics(ia)==SystemCharacteristics{EX,HA,LI})) return std::make_shared<BinarySystem<EX,HA,LI> >(ia)


const binary::Ptr binary::make(InteractionPtr ia)
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

