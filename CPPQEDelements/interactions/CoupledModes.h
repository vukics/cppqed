// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_INTERACTIONS_COUPLEDMODES_H_INCLUDED
#define   CPPQEDELEMENTS_INTERACTIONS_COUPLEDMODES_H_INCLUDED

#include "Mode_.h"

#include "Interaction.h"


namespace coupledmodes {


enum Coupling {
  CM_NX,
  CM_XX,
  CM_XX_RWA
};


template<bool IS_HA, Coupling C=CM_NX>
class Base;


typedef std::shared_ptr<const Base<false> > Ptr;

template<>
class Base<false> : public structure::Interaction<2>
{
public:
  Base(mode::Ptr, mode::Ptr, dcomp);

};

// NOTE! In case when C is CM_NX or CM_XX, only the real part of u is taken into account!
template<Coupling C>
class Base<true,C> : public Base<false>, public quantumoperator::TridiagonalHamiltonian<2,true>
{
public:
  Base(mode::Ptr, mode::Ptr, dcomp u);

};


} // coupledmodes



#define BIG_NAMESPACE_NAME                      coupledmodes
#define BIG_CLASS_NAME                          CoupledModes
#define BIG_ADDITIONAL_PARAMETERS               , dcomp u
#define BIG_ADDITIONAL_PARAMETERS_PASS          ,u
#define BIG_ADDITIONAL_TEMPLATE_PARAMETERS      bool IS_HA, coupledmodes::Coupling C,
#define BIG_ADDITIONAL_TEMPLATE_PARAMETERS_PASS <IS_HA,C>

#include "details_BinaryInteractionGenerator.h"


namespace coupledmodes {


template<Coupling C, typename A, typename F1, typename F2>
const Ptr make(const F1& f1, const F2& f2, dcomp u, const A& =A())
{
  if (abs(u)) return std::make_shared<CoupledModes<true ,C    ,A> >(f1,f2,u);
  else        return std::make_shared<CoupledModes<false,CM_NX,A> >(f1,f2,0);
}


template<Coupling C, typename F1, typename F2>
const Ptr make(const F1& f1, const F2& f2, dcomp u)
{
  return make<C,EmptyAveragingBaseForInteractions>(f1,f2,u);  
}


} // coupledmodes


#endif // CPPQEDELEMENTS_INTERACTIONS_COUPLEDMODES_H_INCLUDED
