// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Highest level driver functions for quantum trajectories}
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION__H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION__H_INCLUDED

#ifndef EVOLUTION_DEFAULT_RANDOM_ENGINE
#define EVOLUTION_DEFAULT_RANDOM_ENGINE cpputils::GSL_RandomEngine
#endif

#include "EnsembleMCWF.h"
#include "Master.tcc"
#include "MCWF_Trajectory.h"

#include "TMP_Tools.h"
#include "Trajectory.tcc"

#include <iosfwd>


using namespace quantumtrajectory;


/// Auxiliary tools for the evolve functions
namespace evolution {

/// Method of evolution for a quantum system
enum Method {
  SINGLE, ///< single \link quantumtrajectory::MCWF_Trajectory MCWF trajectory\endlink
  ENSEMBLE, ///< \link quantumtrajectory::EnsembleMCWF ensemble\endlink of MCWF trajectories
  MASTER ///< Master equation with \link quantumtrajectory::master::Base normal iteration\endlink
};

std::ostream& operator<<(std::ostream&, Method ); ///< output streaming for Method
std::istream& operator>>(std::istream&, Method&); ///< input streaming for Method


/// Aggregate of parameters pertaining to the highest level driver functions for quantum trajectories
/** \copydetails trajectory::ParsRun */
template <typename Base=mcwf::Pars>
struct Pars : public trajectory::ParsRun, public Base {

  Method &evol; ///< the method of evolution
  bool
    &negativity, ///< governs whether entanglement should be calculated in the case of #ENSEMBLE and #MASTER, cf. quantumtrajectory::stream_densityoperator::_, quantumdata::negPT
    &timeAverage; ///< governs whether in the case of #SINGLE, time averaging should be performed (by using quantumtrajectory::TimeAveragingMCWF_Trajectory instead of quantumtrajectory::MCWF_Trajectory)
  double &relaxationTime; ///< the relaxation time in the case when time averaging is desired

  Pars(parameters::Table& p, const std::string& mod="") 
    : ParsRun(p,mod),
      Base(p,mod),
      evol(p.addTitle("Evolution",mod).add("evol",mod,"Evolution mode (single, ensemble, master)",SINGLE)),
      negativity(p.add("negativity",mod,"Calculates negativity in ensemble & master",false)),
      timeAverage(p.add("timeAverage",mod,"Calculates time averages in MCWF trajectory",false)),
      relaxationTime(p.add("relaxationTime",mod,"Relaxation time for time averaging",0.)) {}    
};


template<int RANK>
using LDO_Ptr=typename quantumdata::LazyDensityOperator<RANK>::Ptr;

using StreamedArray=trajectory::StreamedArray<DArray<1>>;

template<int RANK>
using EvolutionReturnType=std::tuple<evolution::LDO_Ptr<RANK>,StreamedArray>;


template<typename V, int RANK>
EvolutionReturnType<RANK>
_(std::shared_ptr<quantumdata::DensityOperator<RANK>> rhoPtr,
  typename structure::QuantumSystem<RANK>::Ptr sys,
  const evolution::Pars<>& pe,
  bool doStreaming, bool returnStreamedArray)
{
  Master<RANK,V> traj(rhoPtr,sys,pe,pe.negativity);

  return {rhoPtr,trajectory::run(traj,pe,doStreaming,returnStreamedArray)};

}


namespace details {

template<typename RandomEngine, typename V, int RANK>
EvolutionReturnType<RANK>
_(std::shared_ptr<quantumdata::StateVector<RANK>> psiPtr,
  typename structure::QuantumSystem<RANK>::Ptr sys,
  const evolution::Pars<>& pe,
  bool doStreaming, bool returnStreamedArray)
{
  if      (pe.evol==evolution::SINGLE) {
    return {psiPtr,trajectory::run(*std::make_shared<MCWF_Trajectory<RANK,RandomEngine>>(psiPtr,sys,pe),pe,doStreaming,returnStreamedArray)};
  }
  else if (pe.evol==evolution::ENSEMBLE) {
    EnsembleMCWF<RANK,RandomEngine,V> traj(psiPtr,sys,pe,pe.negativity);
    return {traj.averaged(),trajectory::run(traj,pe,doStreaming,returnStreamedArray)};
  }
  else {
    return _<V>(std::make_shared<quantumdata::DensityOperator<RANK>>(*psiPtr),sys,pe,doStreaming,returnStreamedArray);
  }

}

} // details


template<typename RandomEngine, typename V, int RANK>
const auto
_(quantumdata::StateVector<RANK>&& psi, ///<[in/out] pure state-vector initial condition
  typename structure::QuantumSystem<RANK>::Ptr sys, ///<[in] the simulated \link structure::QuantumSystem quantum system\endlink
  const evolution::Pars<>& p, ///<[in] parameters of the evolution
  bool doStreaming=true, bool returnStreamedArray=false)
{
  return details::_<RandomEngine,V>(std::make_shared<quantumdata::StateVector<RANK>>(std::move(psi)),sys,p,doStreaming,returnStreamedArray);
}


template<typename RandomEngine, typename V, int RANK>
const auto
_(quantumdata::StateVector<RANK>& psi, typename structure::QuantumSystem<RANK>::Ptr sys, const evolution::Pars<>& p, bool doStreaming=true, bool returnStreamedArray=false)
{
  return details::_<RandomEngine,V>(std::make_shared<quantumdata::StateVector<RANK>>(psi.getArray(),quantumdata::byReference),sys,p,doStreaming,returnStreamedArray);
}


template<typename RandomEngine, int... V, typename SV_OR_DO, typename SYS>
const auto
_(SV_OR_DO&& initial, SYS sys, const evolution::Pars<>& p, bool doStreaming=true, bool returnStreamedArray=false)
{
  return _<RandomEngine,tmptools::Vector<V...>>(std::forward<SV_OR_DO>(initial),sys,p,doStreaming,returnStreamedArray);
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
template<typename V, int RANK>
const auto
evolve(quantumdata::StateVector<RANK>&& psi, ///<[in/out] pure state-vector initial condition
       typename structure::QuantumSystem<RANK>::Ptr sys, ///<[in] the simulated \link structure::QuantumSystem quantum system\endlink
       const evolution::Pars<>& p, ///<[in] parameters of the evolution
       bool doStreaming=true, bool returnStreamedArray=false)
{
  return evolution::_<EVOLUTION_DEFAULT_RANDOM_ENGINE,V,RANK>(std::forward<quantumdata::StateVector<RANK>>(psi),sys,p,doStreaming,returnStreamedArray);
}


/// \overload
template<typename V, int RANK>
const auto
evolve(quantumdata::StateVector<RANK>& psi, typename structure::QuantumSystem<RANK>::Ptr sys, const evolution::Pars<>& p, bool doStreaming=true, bool returnStreamedArray=false)
{
  return evolution::_<EVOLUTION_DEFAULT_RANDOM_ENGINE,V>(std::forward<quantumdata::StateVector<RANK>>(psi),sys,p,doStreaming,returnStreamedArray);
}


/// The prototype function to evolve a Master trajectory from a DensityOperator initial condition
template<typename V, int RANK>
const auto
evolve(quantumdata::DensityOperator<RANK>&& rho, ///<[in/out] density operator initial condition
       typename structure::QuantumSystem<RANK>::Ptr sys, ///<[in] the simulated \link structure::QuantumSystem quantum system\endlink
       const evolution::Pars<>& p, ///<[in] parameters of the evolution
       bool doStreaming=true, bool returnStreamedArray=false);


/// \overload
template<typename V, int RANK>
const auto
evolve(quantumdata::DensityOperator<RANK>& rho, typename structure::QuantumSystem<RANK>::Ptr sys, const evolution::Pars<>& p,
       bool doStreaming=true, bool returnStreamedArray=false)
{
  return evolution::_<V>(std::make_shared<quantumdata::DensityOperator<RANK>>(rho.getArray(),quantumdata::byReference),sys,p,doStreaming,returnStreamedArray);
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
template<int... V, typename SV_OR_DO, typename SYS>
const auto
evolve(SV_OR_DO&& initial,
       SYS sys,
       const evolution::Pars<>& p,
       bool doStreaming=true, bool returnStreamedArray=false)
{
  return evolve<tmptools::Vector<V...>>(std::forward<SV_OR_DO>(initial),sys,p,doStreaming,returnStreamedArray);
}


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION__H_INCLUDED
