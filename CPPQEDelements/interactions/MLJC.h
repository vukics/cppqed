/// \briefFile{Defines interaction elements of the \ref multilevelbundle "MultiLevel bundle" }
// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDELEMENTS_INTERACTIONS_MLJC_H_INCLUDED
#define   CPPQEDELEMENTS_INTERACTIONS_MLJC_H_INCLUDED

#include "MLJCFwd.h"

#include "Mode_.h"
#include "MultiLevel_.h"

#include "Interaction.h"

#include <boost/fusion/algorithm/transformation/transform.hpp>


/// Class representing an elementary coupling term (a \f$g_{ij}\f$ \ref multilevelactualHamiltonian "here") with a compile-time pair \f$i,j\f$ and a runtime complex value
template<int I, int J>
class Coupling : public multilevel::Storage<dcomp>, public tmptools::pair_c<I,J>
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


/** \cond */

#define BIG_NAMESPACE_NAME                      mljc
#define BIG_CLASS_NAME                          MLJC
#define BIG_ADDITIONAL_PARAMETERS               , const mljc::Pars<VC>& p
#define BIG_ADDITIONAL_PARAMETERS_PASS          ,p.gs
#define BIG_ADDITIONAL_TEMPLATE_PARAMETERS      int NL, typename VC,
#define BIG_ADDITIONAL_TEMPLATE_PARAMETERS_PASS <NL,VC>

#include "details_BinaryInteractionGenerator.h"

/** \endcond */

#endif // CPPQEDELEMENTS_INTERACTIONS_MLJC_H_INCLUDED
