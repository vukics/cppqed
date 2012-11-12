// -*- C++ -*-
#ifndef   ELEMENTS_INTERACTIONS_NX_COUPLEDMODES_H_INCLUDED
#define   ELEMENTS_INTERACTIONS_NX_COUPLEDMODES_H_INCLUDED

#include "NX_CoupledModesFwd.h"

#include "Mode_.h"

#include "Interaction.h"

#include <boost/make_shared.hpp>


namespace nxcoupledmodes {


typedef boost::shared_ptr<const Base<false> > Ptr;

template<>
class Base<false> : public structure::Interaction<2>
{
public:
  Base(mode::Ptr, mode::Ptr, double);

};

template<>
class Base<true> : public Base<false>, public structure::TridiagonalHamiltonian<2,true>
{
public:
  Base(mode::Ptr, mode::Ptr, double u);

};


} // nxcoupledmodes



#define BIG_NAMESPACE_NAME                      nxcoupledmodes
#define BIG_CLASS_NAME                          NX_CoupledModes
#define BIG_ADDITIONAL_PARAMETERS               , double u
#define BIG_ADDITIONAL_PARAMETERS_PASS          ,u
#define BIG_ADDITIONAL_TEMPLATE_PARAMETERS      bool IS_HA,
#define BIG_ADDITIONAL_TEMPLATE_PARAMETERS_PASS <IS_HA>

#include "details/BinaryInteractionGenerator.h"


namespace nxcoupledmodes {


template<typename A, typename F1, typename F2>
const Ptr make(const F1& f1, const F2& f2, double u, const A& =A())
{
  if (u) return boost::make_shared<NX_CoupledModes<true ,A> >(f1,f2,u);
  else   return boost::make_shared<NX_CoupledModes<false,A> >(f1,f2,0);
}


template<typename F1, typename F2>
const Ptr make(const F1& f1, const F2& f2, double u)
{
  return make<EmptyAveragingBaseForInteractions>(f1,f2,u);  
}


} // nxcoupledmodes


#endif // ELEMENTS_INTERACTIONS_NX_COUPLEDMODES_H_INCLUDED
