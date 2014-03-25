/// \briefFile{Highest level driver functions for quantum trajectories}
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION__H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION__H_INCLUDED

#include "Evolution_Fwd.h"

#include "LazyDensityOperator.h"
#include "DensityOperatorFwd.h"
#include "StateVectorFwd.h"
#include "MCWF_TrajectoryFwd.h"
#include "ParsMCWF_Trajectory.h"
#include "QuantumSystemFwd.h"

#include "SmartPtr.h"
#include "TMP_Tools.h"

#include <iosfwd>



using namespace quantumtrajectory;


/// Auxiliary tools for the evolve functions
namespace evolution {

/// Method of evolution for a quantum system
enum Method {
  SINGLE, ///< single \link quantumtrajectory::MCWF_Trajectory MCWF trajectory\endlink
  ENSEMBLE, ///< \link quantumtrajectory::EnsembleMCWF ensemble\endlink of MCWF trajectories
  MASTER, ///< Master equation with \link quantumtrajectory::master::Base normal iteration\endlink
  MASTER_FAST ///< Master equation with \link quantumtrajectory::master::BaseFast “fast” iteration\endlink
};

std::ostream& operator<<(std::ostream&, Method ); ///< output streaming for Method
std::istream& operator>>(std::istream&, Method&); ///< input streaming for Method


/// Aggregate of parameters pertaining to the highest level driver functions for quantum trajectories
/** \copydetails trajectory::ParsRun */
struct Pars : public trajectory::ParsRun, public mcwf::Pars {

  Method &evol; ///< the method of evolution
  bool
    &negativity, ///< governs whether entanglement should be calculated in the case of #ENSEMBLE, #MASTER, and #MASTER_FAST, cf. quantumtrajectory::display_densityoperator::_, quantumdata::negPT
    &timeAverage; ///< governs whether in the case of #SINGLE, time averaging should be performed (by using quantumtrajectory::TimeAveragingMCWF_Trajectory instead of quantumtrajectory::MCWF_Trajectory)
  double &relaxationTime; ///< the relaxation time in the case when time averaging is desired

  Pars(parameters::ParameterTable& p, const std::string& mod="");

};


/// Dispatcher returning a quantumtrajectory::MCWF_Trajectory or quantumtrajectory::TimeAveragingMCWF_Trajectory instant, depending on the last argument (cf. Pars::timeAverage)
/**
 * \tparamRANK
 * \tparam SYS the object representing the quantum system to be simulated (similar idea as in quantumtrajectory::Master::Master)
 */
template<int RANK, typename SYS>
const boost::shared_ptr<MCWF_Trajectory<RANK> > makeMCWF(quantumdata::StateVector<RANK>&, const SYS&, const Pars&);

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
const typename quantumdata::LazyDensityOperator<RANK>::Ptr
evolve(quantumdata::StateVector<RANK>& psi, ///<[in/out] pure state-vector initial condition
       typename structure::QuantumSystem<RANK>::Ptr sys, ///<[in] the simulated \link structure::QuantumSystem quantum system\endlink
       const evolution::Pars& p ///<[in] parameters of the evolution
       );


/// The prototype function to evolve a Master trajectory from a DensityOperator initial condition
template<typename V, int RANK>
const typename quantumdata::LazyDensityOperator<RANK>::Ptr
evolve(quantumdata::DensityOperator<RANK>& rho, ///<[in/out] density operator initial condition
       typename structure::QuantumSystem<RANK>::Ptr sys, ///<[in] the simulated \link structure::QuantumSystem quantum system\endlink
       const evolution::Pars& p ///<[in] parameters of the evolution
       );


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
inline
const typename quantumdata::LazyDensityOperator<SV_OR_DO::N_RANK>::Ptr
evolve(SV_OR_DO& initial,
       const SYS& sys,
       const evolution::Pars& p)
{
  return evolve<tmptools::Vector<V...> >(initial,cpputils::sharedPointerize(sys),p);
}


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION__H_INCLUDED
