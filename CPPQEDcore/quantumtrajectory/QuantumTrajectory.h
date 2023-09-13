// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#pragma once

#include "QuantumSystemDynamics.h"

#include "EntanglementMeasures.h"

#include "Trajectory.h"



namespace quantumtrajectory {

using EntanglementMeasuresSwitch = std::bitset<3>;

  
/// Forwards to trajectory::initialTimeStep, with the highest frequency of the system taken as structure::QuantumSystem::highestFrequency
double initialTimeStep(const ::structure::system_frequency_store auto& freqs)
{
  return ::cppqedutils::trajectory::initialTimeStep(::structure::highestFrequency(freqs));
}


/*
template<size_t RANK>
std::ostream& streamCharacteristics(const ::structure::quantum_system_dynamics<RANK> auto& qsd, std::ostream& os)
{
  return
    hana::fold(getLi(qsd),
    hana::fold(getHa(qsd),os<<"System characteristics:\nHamiltonian terms: ", [] (std::ostream& os, const auto& he) {
      {
        using T=decltype(termOrPropagator(he));
        if constexpr (::structure::one_time_dependent_propagator<T,RANK> || ::structure::two_time_dependent_propagator<T,RANK>)
          os<<"*";
      }
      return os<<label(he)<<", ";
    })<<"Liouvillian ", [&]<typename LI> (std::ostream& os, const LI& lindblad) {
      if constexpr(::structure::lindblad_with_rate<LI,RANK>) os<<"*";
      if constexpr(::structure::lindblad_with_superoperator<LI,RANK>) os<<"!";
      return os<<label(lindblad);
    });
}
*/

/// Wraps common functionality of Master & EnsembleMCWF concerning calculating temporal data points on the basis of density operators
/**
 * This comprises
 * - keeping a structure::QSD instant and calling structure::Averaged::stream
 * - performing calculation of entanglement measures if needed
 * - extending key with entanglement measures when needed
 */
template<size_t RANK,
         ::structure::quantum_system_dynamics<RANK> QSD,
         auto axesOfSubsystem // the axes belonging to one of the subsystems defined for entanglementMeasuresCalculation
         >
class TDP_DensityOperator/*
{
public:
  using DensityOperator=quantumdata::DensityOperator<RANK>;

  TDP_DensityOperator(::structure::AveragedPtr<RANK> av, EntanglementMeasuresSwitch ems) : av_(av), ems_(ems) {}

  StreamReturnType operator()(double t, const DensityOperator& rho, std::ostream& os, int precision) const 
  {
    auto res{structure::stream(av_,t,rho,os,precision)};
    auto & averages{std::get<1>(res)};
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
    return {os,averages};
  }

  std::ostream& streamKey(std::ostream& os, size_t& i) const
  {
    if (av_) av_->streamKey(os,i); 
    if constexpr ( !isV_empty ) {
      if (ems_.any()) os<<"Trajectory\n";
      if (ems_[0]) os<<i++<<". negativity\n";
      if (ems_[1]) os<<i++<<". mutual information\n";
      if (ems_[2]) os<<i++<<". purity of partial trace\n";
    }
    return os;
  }

private:
  const ::structure::AveragedPtr<RANK> av_ ;

  const EntanglementMeasuresSwitch ems_;

  static constexpr bool isV_empty=boost::mpl::empty<V>::value;
  
}*/;

} // quantumtrajectory

