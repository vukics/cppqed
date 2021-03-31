// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Highest level driver functions for quantum trajectories}
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION__H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION__H_INCLUDED

/// The evolve methods use these engines, but they can still be modified in a script by defining these macros before including Evolution.h
#ifndef QUANTUM_EVOLUTION_DEFAULT_RANDOM_ENGINE
#define QUANTUM_EVOLUTION_DEFAULT_RANDOM_ENGINE pcg64 // XoshiroCpp::Xoshiro256PlusPlus // cppqedutils::GSL_RandomEngine
#endif

#ifndef QUANTUM_EVOLUTION_DEFAULT_ODE_ENGINE
#define QUANTUM_EVOLUTION_DEFAULT_ODE_ENGINE cppqedutils::ODE_EngineBoost
#endif

#include "EvolutionMethod.h"
#include "Master.h"
#include "MCWF_Trajectory.h"

#include "TMP_Tools.h"

#include <iosfwd>


using namespace quantumtrajectory;


/// Auxiliary tools for the evolve functions
namespace evolution {


/// Aggregate of parameters pertaining to the highest level driver functions for quantum trajectories
/** \copydetails trajectory::ParsRun */
template <typename Base=mcwf::Pars<QUANTUM_EVOLUTION_DEFAULT_RANDOM_ENGINE>>
struct Pars : public cppqedutils::trajectory::Pars<Base> {

  Method &evol; ///< the method of evolution
  bool &negativity; ///< governs whether entanglement should be calculated in the case of #ENSEMBLE and #MASTER, cf. quantumtrajectory::stream_densityoperator::_, quantumdata::negPT

  Pars(parameters::Table& p, const std::string& mod="") 
    : cppqedutils::trajectory::Pars<Base>(p,mod),
      evol(p.addTitle("Evolution",mod).add("evol",mod,"Evolution mode (single, ensemble, master)",SINGLE)),
      negativity(p.add("negativity",mod,"Calculates negativity in ensemble & master",false))
      {}

};


using TemporalStreamedArray=cppqedutils::trajectory::TemporalStreamedArray<DArray<1>>;

} // evolution


/// The prototype function to evolve a Master trajectory from a DensityOperator initial condition
template <template <typename StateType> class ODE_Engine, typename V, typename DO, typename SYS, typename Parameters>
auto
evolveMaster(DO&& rho, ///<[in/out] density operator initial condition
             SYS&& sys, ///<[in] the simulated \link structure::QuantumSystem quantum system\endlink
             const Parameters& p, ///<[in] parameters of the evolution
             bool doStreaming=true, bool returnStreamedArray=false)
{
  static constexpr auto RANK=std::decay_t<DO>::N_RANK;
  
  return cppqedutils::run(master::make<ODE_Engine<quantumdata::DensityOperatorLow<RANK>>,V>
    (std::forward<SYS>(sys),
     std::forward<quantumdata::DensityOperator<RANK>>(rho),
     p,p.negativity),
     p,doStreaming,returnStreamedArray);
}


template<template <typename StateType> class ODE_Engine, int... V, typename DO, typename SYS, typename Parameters>
auto
evolveMaster(DO&& rho, SYS&& sys, const Parameters& p, bool doStreaming=true, bool returnStreamedArray=false)
{
  return evolveMaster<ODE_Engine,tmptools::Vector<V...>>(std::forward<DO>(rho),std::forward<SYS>(sys),p,doStreaming,returnStreamedArray);
}


namespace evolution {


namespace details {
  

template<auto RANK,
         template <typename StateType> class ODE_Engine,
         typename RandomEngine,
         typename V,
         typename SV_OR_DO,
         typename Parameters>
std::enable_if_t<std::is_same_v<std::decay_t<SV_OR_DO>,quantumdata::DensityOperator<RANK>> ||
                 std::is_same_v<std::decay_t<SV_OR_DO>,quantumdata::StateVector<RANK>>,
                 TemporalStreamedArray>
_(SV_OR_DO&& initial, ///<[in/out] pure state-vector initial condition
  structure::QuantumSystemPtr<RANK> sys, ///<[in] the simulated \link structure::QuantumSystem quantum system\endlink
  const Parameters& pe, ///<[in] parameters of the evolution
  bool doStreaming=true, bool returnStreamedArray=false)
{
  using OE=ODE_Engine<quantumdata::StateVectorLow<RANK>>;
  
  // DensityOperator initial condition
  if constexpr (std::is_same_v<std::decay_t<SV_OR_DO>,quantumdata::DensityOperator<RANK>>) {
    if (pe.evol==evolution::SINGLE) {
      throw std::runtime_error("Single MCWF trajectory from DensityOperator initial condition");
    }
    else if (pe.evol==evolution::ENSEMBLE) {
      /// here, it’s intentional that `initial` is converted into an lvalue
      return cppqedutils::run(mcwf::makeEnsemble<OE,RandomEngine,V>(sys,initial,pe,pe.negativity),pe,doStreaming,returnStreamedArray);
    }
    else {
      return ::evolveMaster<ODE_Engine,V>(std::forward<quantumdata::DensityOperator<RANK>>(initial),sys,pe,doStreaming,returnStreamedArray);
    }
  }
  // StateVector initial condition
  else {
    if (pe.evol==evolution::SINGLE) {
      return cppqedutils::run(mcwf::make<OE,RandomEngine>(sys,std::forward<quantumdata::StateVector<RANK>>(initial),pe),
                              pe,doStreaming,returnStreamedArray);
    }
    else if (pe.evol==evolution::ENSEMBLE) {
      /// here, it’s intentional that `initial` is converted into an lvalue
      return cppqedutils::run(mcwf::makeEnsemble<OE,RandomEngine,V>(sys,initial,pe,pe.negativity),pe,doStreaming,returnStreamedArray);
    }
    else {
      return ::evolveMaster<ODE_Engine,V>(quantumdata::DensityOperator<RANK>(initial),sys,pe,doStreaming,returnStreamedArray);
    }    
  }

}


} // details


template<template <typename StateType> class ODE_Engine, typename RandomEngine,
         typename V, typename SV_OR_DO, typename SYS, typename Parameters>
auto _(SV_OR_DO&& initial, SYS&& sys, const Parameters& p, bool doStreaming=true, bool returnStreamedArray=false)
{
  static constexpr auto RANK=std::decay_t<SV_OR_DO>::N_RANK;

  return details::_<RANK,ODE_Engine,RandomEngine,V>(std::forward<SV_OR_DO>(initial),
                                         std::forward<SYS>(sys),
                                         p,doStreaming,returnStreamedArray);
}


template<template <typename StateType> class ODE_Engine, typename RandomEngine,
         int... V, typename SV_OR_DO, typename SYS, typename Parameters>
auto _(SV_OR_DO&& initial, SYS&& sys, const Parameters& p, bool doStreaming=true, bool returnStreamedArray=false)
{
  return _<ODE_Engine,RandomEngine,tmptools::Vector<V...>>(
    std::forward<SV_OR_DO>(initial),
    std::forward<SYS>(sys),
    p,doStreaming,returnStreamedArray);
}


} // evolution


/// The prototype function to evolve a quantumtrajectory from a pure StateVector initial condition
/**
 * Basically a dispatcher invoking trajectory::run, with the difference that it also creates and stores the necessary \link #quantumtrajectory quantum trajectory\endlink
 * 
 * \tparam V has the same function as the template parameter `V` in quantumdata::negPT (cannot be inferred)
 * \tparamRANK (inferred from the 1st function argument)
 *
 * \return the final state of the evolution as a quantumdata::LazyDensityOperator
 * (for evolution::SINGLE this will be a *copy* of the evolved \link quantumdata::StateVector state vector\endlink,
 * while for evolution::MASTER and evolution::ENSEMBLE a \link quantumdata::DensityOperator density operator\endlink)
 *
 * \note The evolved state can always be accessed if the initial state is an lvalue.
 * In the case of Ensemble evolution, cppqedutils::run has to be used to access the evolved state.
 */
template<typename V, typename SV_OR_DO, typename SYS, typename Parameters>
auto
evolve(SV_OR_DO&& initial, ///<[in/out] pure state-vector initial condition
       SYS sys, ///<[in] the simulated \link structure::QuantumSystem quantum system\endlink
       const Parameters& p, ///<[in] parameters of the evolution
       bool doStreaming=true, bool returnStreamedArray=false)
{
  return evolution::_<QUANTUM_EVOLUTION_DEFAULT_ODE_ENGINE,QUANTUM_EVOLUTION_DEFAULT_RANDOM_ENGINE,V>(
    std::forward<SV_OR_DO>(initial),
    std::forward<SYS>(sys),
    p,doStreaming,returnStreamedArray);
}


/// \overload
/**
 * Adds also the syntactic sugar of being able to write
 *
 *     evolve<2,0,4>(psi,sys,p)
 *
 * instead of
 *
 *     evolve<tmptools::Vector<2,0,4> >(psi,sys,p)
 *
 */
template<int... V, typename SV_OR_DO, typename SYS, typename Parameters>
auto
evolve(SV_OR_DO&& initial,
       SYS&& sys,
       const Parameters& p,
       bool doStreaming=true, bool returnStreamedArray=false)
{
  return evolve<tmptools::Vector<V...>>(std::forward<SV_OR_DO>(initial),
                                        std::forward<SYS>(sys),
                                        p,doStreaming,returnStreamedArray);
}


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION__H_INCLUDED
