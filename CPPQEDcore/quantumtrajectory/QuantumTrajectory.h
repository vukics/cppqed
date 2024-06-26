// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "ExpectationValues.h"
#include "Hamiltonian.h"

#include "EntanglementMeasures.h"

#include "Trajectory.h"



namespace quantumtrajectory {

using namespace structure;


using EntanglementMeasuresSwitch = std::bitset<3>;

  
template <size_t RANK>
double initialTimeStep(const hamiltonian<RANK> auto& h, Dimensions<RANK> d)
{
  StateVector<RANK> dpsidt{d,zeroInit}, psi{d,noInit}; for (dcomp& v : psi.mutableView().dataView) v=1.;
  applyHamiltonian(h,0.,psi,dpsidt.mutableView(),0.);
  return .1/std::ranges::max(dpsidt.dataView | std::views::transform( [] (dcomp v) {return abs(v);} ) ) ; // tenth of the inverse of the largest frequency
}


/*
template<size_t RANK>
std::ostream& streamCharacteristics(const quantum_system_dynamics<RANK> auto& qsd, std::ostream& os)
{
  return
    hana::fold(getLi(qsd),
    hana::fold(getHa(qsd),os<<"System characteristics:\nHamiltonian terms: ", [] (std::ostream& os, const auto& he) {
      {
        using T=decltype(termOrPropagator(he));
        if constexpr (one_time_dependent_propagator<T,RANK> || two_time_dependent_propagator<T,RANK>)
          os<<"*";
      }
      return os<<label(he)<<", ";
    })<<"Liouvillian ", [&]<typename LI> (std::ostream& os, const LI& lindblad) {
      if constexpr(lindblad_with_rate<LI,RANK>) os<<"*";
      if constexpr(lindblad_with_superoperator<LI,RANK>) os<<"!";
      return os<<label(lindblad);
    });
}
*/

/// Wraps common functionality of Master & EnsembleMCWF concerning calculating temporal data points on the basis of density operators
/**
 * This comprises
 * - keeping a QSD instant and calling Averaged::stream
 * - performing calculation of entanglement measures if needed
 * - extending key with entanglement measures when needed
 */
template<size_t RANK,
         expectation_values<RANK> EV/*,
         auto axesOfSubsystem */ // the axes belonging to one of the subsystems defined for entanglementMeasuresCalculation
         >
class TDP_DensityOperator
{
public:
  TDP_DensityOperator(auto&& ev /*, EntanglementMeasuresSwitch ems*/) : ev{std::forward<decltype(ev)>(ev)} {}

  EV ev;
//  const EntanglementMeasuresSwitch ems_;


  auto operator()(double t, const DensityOperator<RANK>& rho) const
  {
    return calculateExpectationValues<RANK>(ev,t,LDO<DensityOperator,RANK>(rho));
/*    auto & averages{std::get<1>(res)};
    if constexpr ( !isV_empty ) {
      if (ems_[0]) {
        auto n{negPT(rho,V{})};
        os<<'\t'<<FormDouble(precision)(n);
        averages.resizeAndPreserve(averages.size()+1); averages(averages.ubound(0))=n;
      }
      if (ems_[1]) {
        auto mi{mutualInformation(rho,V{})};
        os<<'\t'<<FormDouble(precision)(mi);
        averages.resizeAndPreserve(averages.size()+1); averages(averages.ubound(0))=mi;
      }
      if (ems_[2]) {
        auto p{purityOfPartialTrace(rho,V{})};
        os<<'\t'<<FormDouble(precision)(p);
        averages.resizeAndPreserve(averages.size()+1); averages(averages.ubound(0))=p;
      }

    }
    return {os,averages};*/
  }

};

} // quantumtrajectory

