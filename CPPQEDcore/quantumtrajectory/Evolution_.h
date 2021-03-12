// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Highest level driver functions for quantum trajectories}
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION__H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION__H_INCLUDED

#ifndef EVOLUTION_DEFAULT_RANDOM_ENGINE
#define EVOLUTION_DEFAULT_RANDOM_ENGINE pcg64 // XoshiroCpp::Xoshiro256PlusPlus // cppqedutils::GSL_RandomEngine
#endif

#include "EnsembleMCWF.h"
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
template <typename Base=mcwf::Pars<EVOLUTION_DEFAULT_RANDOM_ENGINE>>
struct Pars : public trajectory::ParsRun, public Base {

  Method &evol; ///< the method of evolution
  bool &negativity; ///< governs whether entanglement should be calculated in the case of #ENSEMBLE and #MASTER, cf. quantumtrajectory::stream_densityoperator::_, quantumdata::negPT

  Pars(parameters::Table& p, const std::string& mod="") 
    : ParsRun(p,mod),
      Base(p,mod),
      evol(p.addTitle("Evolution",mod).add("evol",mod,"Evolution mode (single, ensemble, master)",SINGLE)),
      negativity(p.add("negativity",mod,"Calculates negativity in ensemble & master",false))
      {}

};


template<auto RANK>
using LDO_Ptr=typename quantumdata::LazyDensityOperator<RANK>::Ptr;

using TemporalStreamedArray=trajectory::TemporalStreamedArray<DArray<1>>;

template<auto RANK>
using ReturnType=std::tuple<evolution::LDO_Ptr<RANK>,TemporalStreamedArray>;

} // evolution


/// The prototype function to evolve a Master trajectory from a DensityOperator initial condition
template<typename V, typename DO, typename SYS, typename Parameters>
auto
evolveMaster(DO&& rho, ///<[in/out] density operator initial condition
             SYS&& sys, ///<[in] the simulated \link structure::QuantumSystem quantum system\endlink
             const Parameters& p, ///<[in] parameters of the evolution
             bool doStreaming=true, bool returnStreamedArray=false)
{
  static constexpr auto RANK=std::decay_t<DO>::N_RANK;
  
  Master<RANK,V> traj(std::forward<quantumdata::DensityOperator<RANK>>(rho),
                      std::forward<SYS>(sys),
                      p,p.negativity);

  return evolution::ReturnType<RANK>{std::make_shared<quantumdata::DensityOperator<RANK>>(std::move(rho)),
                                     trajectory::run(static_cast<typename Master<RANK,V>::TrajectoryBase&>(traj),p,doStreaming,returnStreamedArray)};
                                     // this ugly static_cast is needed because the member function `stream` is found in two separate bases of Master
}


template<int... V, typename DO, typename SYS, typename Parameters>
auto
evolveMaster(DO&& rho, SYS&& sys, const Parameters& p, bool doStreaming=true, bool returnStreamedArray=false)
{
  return evolveMaster(std::forward<DO>(rho),std::forward<SYS>(sys),p,doStreaming,returnStreamedArray);
}


namespace evolution {


namespace details {
  

template<auto RANK, typename RandomEngine, typename V, typename SV_OR_DO, typename Parameters>
std::enable_if_t<std::is_same_v<std::decay_t<SV_OR_DO>,quantumdata::DensityOperator<RANK>> || std::is_same_v<std::decay_t<SV_OR_DO>,quantumdata::StateVector<RANK>>,
                 ReturnType<RANK>>
_(SV_OR_DO&& initial, ///<[in/out] pure state-vector initial condition
  typename structure::QuantumSystem<RANK>::Ptr sys, ///<[in] the simulated \link structure::QuantumSystem quantum system\endlink
  const Parameters& pe, ///<[in] parameters of the evolution
  bool doStreaming=true, bool returnStreamedArray=false)
{
  // DensityOperator initial condition
  if constexpr (std::is_same_v<std::decay_t<SV_OR_DO>,quantumdata::DensityOperator<RANK>>) {
    if (pe.evol==evolution::SINGLE) {
      throw std::runtime_error("Single MCWF trajectory from DensityOperator initial condition");
    }
    else if (pe.evol==evolution::ENSEMBLE) {
      EnsembleMCWF<RANK,RandomEngine,V> traj{initial,sys,pe,pe.negativity};
      return {traj.averaged(),
              trajectory::run(traj,pe,doStreaming,returnStreamedArray)};
    }
    else {
      return ::evolveMaster<V>(std::forward<quantumdata::DensityOperator<RANK>>(initial),sys,pe,
                               doStreaming,returnStreamedArray);
    }
  }
  // StateVector initial condition
  else {
    if (pe.evol==evolution::SINGLE) {
      return {std::make_shared<quantumdata::StateVector<RANK>>(std::move(initial)),
              trajectory::run(static_cast<typename MCWF_Trajectory<RANK,RandomEngine>::TrajectoryBase&&>(
                                MCWF_Trajectory<RANK,RandomEngine>(std::forward<quantumdata::StateVector<RANK>>(initial),sys,pe)),
                              pe,doStreaming,returnStreamedArray)};
                              // this ugly static cast is needed because the member function `stream` is found in two separate bases of MCWF_Trajectory
    }
    else if (pe.evol==evolution::ENSEMBLE) {
      EnsembleMCWF<RANK,RandomEngine,V> traj{initial,sys,pe,pe.negativity};
      return {std::make_shared<quantumdata::DensityOperator<RANK>>(std::move(traj.averaged())),
              trajectory::run(traj,pe,doStreaming,returnStreamedArray)};
    }
    else {
      return ::evolveMaster<V>(quantumdata::DensityOperator<RANK>(initial),sys,pe,
                               doStreaming,returnStreamedArray);
    }    
  }

}


}

template<typename RandomEngine, typename V, typename SV_OR_DO, typename SYS, typename Parameters>
auto
_(SV_OR_DO&& initial, SYS&& sys, const Parameters& p, bool doStreaming=true, bool returnStreamedArray=false)
{
  static constexpr auto RANK=std::decay_t<SV_OR_DO>::N_RANK;

  return details::_<RANK,RandomEngine,V>(std::forward<SV_OR_DO>(initial),
                                         std::forward<SYS>(sys),
                                         p,doStreaming,returnStreamedArray);
}


template<typename RandomEngine, int... V, typename SV_OR_DO, typename SYS, typename Parameters>
auto
_(SV_OR_DO&& initial, SYS&& sys, const Parameters& p, bool doStreaming=true, bool returnStreamedArray=false)
{
  return _<RandomEngine,tmptools::Vector<V...>>(std::forward<SV_OR_DO>(initial),
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
 */
template<typename V, typename SV_OR_DO, typename SYS, typename Parameters>
auto
evolve(SV_OR_DO&& initial, ///<[in/out] pure state-vector initial condition
       SYS sys, ///<[in] the simulated \link structure::QuantumSystem quantum system\endlink
       const Parameters& p, ///<[in] parameters of the evolution
       bool doStreaming=true, bool returnStreamedArray=false)
{
  return evolution::_<EVOLUTION_DEFAULT_RANDOM_ENGINE,V>(std::forward<SV_OR_DO>(initial),
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
