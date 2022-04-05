/// \briefFile{Defines interaction elements of the \ref multilevelbundle "MultiLevel bundle" }
// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_INTERACTIONS_MLJC_H_INCLUDED
#define   CPPQEDELEMENTS_INTERACTIONS_MLJC_H_INCLUDED

#include "Mode_.h"
#include "MultiLevel.h"
#include "ParsMLJC.h"

#include "Interaction.h"

#include "Sigma.h"


/// Class representing an elementary coupling term (a \f$g_{ij}\f$ \ref multilevelactualHamiltonian "here") with a compile-time pair \f$i,j\f$ and a runtime complex value
template<int I, int J> using Coupling = multilevel::DynamicsPair<dcomp,I,J>;


namespace mljc {

template<int NL, typename VC>
class Base : public structure::Interaction<2>, public structure::Hamiltonian<2>
{
public:
  typedef mode::Tridiagonal Tridiagonal;

  typedef typename MultiLevelBase<NL>::Ptr MultiLevelPtr;

  Base(MultiLevelPtr ml, mode::Ptr mode, const VC& gs) 
    : structure::Interaction<2>({ml,mode},{},[&] {
        structure::DynamicsBase::ComplexFreqs res;
        for_each(gs,multilevel::elementaryComplexFreqs(res,"g"));
        return res;
      } () ),
      mds_(modeDynamicsMaker(gs,ml,mode))
  {
    getParsStream()<<"Multi-Level Jaynes-Cummings\n";
  }

private:
  void addContribution_v(double t, const structure::StateVectorLow<2>& psi, structure::StateVectorLow<2>& dpsidt, double t0) const override
  {
    boost::fusion::for_each(mds_, [&] (const auto& modeDynamics) {
      using MD=std::decay_t<decltype(modeDynamics)>;
      
      quantumoperator::Sigma<MD::first,MD::second> sigma;

      Tridiagonal& 
        a      (modeDynamics.a      .propagate(t-t0)),
        adagger(modeDynamics.adagger.propagate(t-t0));

      // Note that the multiplication with -g and conj(g) has already
      // been taken care of by ModeDynamics above
      
      (sigma         *adagger).apply(psi,dpsidt);
      (sigma.dagger()*a      ).apply(psi,dpsidt);
    });
  }

  template<int N1, int N2>
  struct ModeDynamics : public tmptools::pair_c<N1,N2>
  {
    ModeDynamics(MultiLevelPtr ml, mode::Ptr mode, dcomp g)
      : a      (-g*mode::aop(mode)),
        adagger((g*mode::aop(mode)).dagger())
    {
      if (dynamic_cast<const multilevel::Exact<NL>*>(ml.get())) throw std::runtime_error("multilevel::Exact not implemented");
    }

    mutable Tridiagonal a, adagger;

  };


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

  static auto modeDynamicsMaker(const VC& gs, MultiLevelPtr ml, mode::Ptr mode)
  {
    return boost::fusion::transform(gs,CouplingToModeDynamics(ml,mode));
    /*Here, the generic lambda doesn’t seem to be a good enough solution, since we need the result of fusion::transform as well (cf. below)
     * [=] (const auto& coupling) -> ModeDynamics<std::decay_t<decltype(coupling)>::first,std::decay_t<decltype(coupling)>::second> {
      using C=std::decay_t<decltype(coupling)>;
      return ModeDynamics<C::first,C::second>(ml,mode,coupling.get());
    }*/
  }  
  
  const std::invoke_result_t<decltype(&Base::modeDynamicsMaker), const VC&, MultiLevelPtr, mode::Ptr> mds_;

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
