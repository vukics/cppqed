// -*- C++ -*-
#ifndef   MULTI_LEVEL_JAYNES_CUMMINGS_INCLUDED
#define   MULTI_LEVEL_JAYNES_CUMMINGS_INCLUDED

#include "MLJCFwd.h"

#include "Mode_.h"
#include "MultiLevel.h"

#include "Interaction.h"

#include "details/DispatchFreeType.h"

#include<boost/fusion/algorithm/transformation/transform.hpp>


template<int N1, int N2>
class Coupling : public multilevel::Storage<dcomp>, public tmptools::pair_c<N1,N2>
{
public:
  typedef multilevel::Storage<dcomp> Base;

  Coupling(const dcomp& value) : Base(value) {}
  Coupling() : Base(dcomp()) {}

};


namespace mljc {

template<int NL, typename VC>
class Base : public structure::Interaction<2>, public structure::Hamiltonian<2>
{
public:
  typedef mode::Tridiagonal Tridiagonal;

  Base(const MultiLevelBase<NL>*, const ModeBase*, const VC&);

private:
  class ElementaryCoupling;

  void addContribution(double, const StateVectorLow&, StateVectorLow&, double) const; 

  template<int,int>
  class ModeDynamics;


  class CouplingToModeDynamics
  {
  public:
    CouplingToModeDynamics(const MultiLevelBase<NL>& ml, const ModeBase& mode) : ml_(ml), mode_(mode) {}
    
    template<typename> struct result;

    template<typename T, typename C>
    struct result<T(const C&)> : mpl::identity<ModeDynamics<C::first,C::second> > {};

    template<typename C>
    const ModeDynamics<C::first,C::second>
    operator()(const C& coupling) const 
    {
      return ModeDynamics<C::first,C::second>(ml_,mode_,coupling.get());
    }
    
  private:
    const MultiLevelBase<NL>& ml_  ;
    const ModeBase&           mode_;
    
  };


  typedef typename boost::fusion::result_of::transform<VC const,CouplingToModeDynamics>::type ModeDynamicss;
  
  const ModeDynamicss mds_;

};


} // mljc




template<int NL, typename VC>
class MLJC : public mljc::Base<NL,VC>
{
public:
  typedef mljc::Base<NL,VC> Base;

  template<typename MLB, typename MB>
  MLJC(const MLB& ml, const MB& m, const mljc::Pars<VC>& p)
    : Base(dispatchFreeType(ml),dispatchFreeType(m),p.gs) {}
};


#endif // MULTI_LEVEL_JAYNES_CUMMINGS_INCLUDED
