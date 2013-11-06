// -*- C++ -*-
#ifndef   INTERACTIONS_MLJC_H_INCLUDED
#define   INTERACTIONS_MLJC_H_INCLUDED

#include "MLJCFwd.h"

#include "Mode_.h"
#include "MultiLevel_.h"

#include "Interaction.h"

#include <boost/fusion/algorithm/transformation/transform.hpp>


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

  typedef typename MultiLevelBase<NL>::Ptr MultiLevelPtr;

  Base(MultiLevelPtr, mode::Ptr, const VC&);

private:
  class ElementaryCoupling;

  void addContribution_v(double, const StateVectorLow&, StateVectorLow&, double) const; 

  template<int,int>
  struct ModeDynamics;


  class CouplingToModeDynamics
  {
  public:
    CouplingToModeDynamics(MultiLevelPtr ml, mode::Ptr mode) : ml_(ml), mode_(mode) {}
    
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
    const MultiLevelPtr ml_  ;
    const mode::Ptr     mode_;
    
  };


  typedef typename boost::fusion::result_of::transform<VC const,CouplingToModeDynamics>::type ModeDynamicss;
  
  const ModeDynamicss mds_;

};


} // mljc



#define BIG_NAMESPACE_NAME                      mljc
#define BIG_CLASS_NAME                          MLJC
#define BIG_ADDITIONAL_PARAMETERS               , const mljc::Pars<VC>& p
#define BIG_ADDITIONAL_PARAMETERS_PASS          ,p.gs
#define BIG_ADDITIONAL_TEMPLATE_PARAMETERS      int NL, typename VC,
#define BIG_ADDITIONAL_TEMPLATE_PARAMETERS_PASS <NL,VC>

#include "details_BinaryInteractionGenerator.h"


#endif // INTERACTIONS_MLJC_H_INCLUDED
