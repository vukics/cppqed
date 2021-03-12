// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_QUANTUMTRAJECTORY_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_QUANTUMTRAJECTORY_H_INCLUDED

#include "QuantumSystem.h"
#include "Structure.h"

#include "Trajectory.h"


/// Comprises modules representing trajectory drivers for simulating quantum systems
namespace quantumtrajectory {


using StreamReturnType=std::tuple<std::ostream&,typename structure::AveragedCommon::Averages>;

  
/// Forwards to trajectory::initialTimeStep, with the highest frequency of the system taken as structure::QuantumSystem::highestFrequency
template<int RANK>
inline double initialTimeStep(typename structure::QuantumSystem<RANK>::Ptr qs)
{
  return cppqedutils::trajectory::initialTimeStep(qs->highestFrequency());
}

/// Class hosting common code of MCWF_Trajectory and Master.
/**
 * It contains the structure::QuantumSystemWrapper and a trajectory::Adaptive::readStateMore_v implementation which sets t0_.
 * Because this class has to be inserted in the inheritence chain of both MCWF_Trajectory and Master, the base class of QuantumTrajectory is passed
 * as a template parameter `BASE` and the constructor uses perfect forwarding semantics.
 * 
 * \tparamRANK
 * \tparam BASE practically either a trajectory::Adaptive (for Master) or a trajectory::Stochastic (for MCWF_Trajectory)
 */
template<int RANK, typename BASE>
class QuantumTrajectory : public BASE, public structure::QuantumSystemWrapper<RANK,true>
{
public:
  using TrajectoryBase=BASE;
  
protected:
  QuantumTrajectory(QuantumTrajectory&&) = default; QuantumTrajectory& operator=(QuantumTrajectory&&) = default;

  /// Constructor forwarding to `BASE` and QuantumSystemWrapper
  template<typename... ArgumentPack>
  QuantumTrajectory(typename structure::QuantumSystem<RANK>::Ptr qs, bool isNoisy, ArgumentPack&&... argumentPack) 
    : BASE(std::forward<ArgumentPack>(argumentPack)...),
      structure::QuantumSystemWrapper<RANK,true>(qs,isNoisy),
      t0_(0) {}

  /// Forwards to `BASE`, but also sets \link getT0 `t0`\endlink
  cppqedutils::iarchive&  readStateMore_v(cppqedutils::iarchive &iar) override
    { BASE::readStateMore_v(iar); if (this->getEx()) setT0(); return iar; }

  /// The time instant of the beginning of the current time step
  /**
   * If the simulated quantum system is derived from structure::Exact, this is a very crucial piece of data,
   * since this is the time instant when the interaction picture and the normal picture coincide.
   * 
   * Hence, careful bookkeeping is necessary to keep this correct through any change of the trajectory state (e.g. state i/o).
   * 
   * \see the parameter \f$t_0\f$ in structure::Exact::actWithU
   * 
   */
  double getT0() const {return t0_;}
  
  /// \copydoc getT0
  void setT0(double t0) const {t0_=t0;}
  
  /// \copybrief getT0 Sets to the current time.
  void setT0() const {t0_=BASE::getTime();}

  /// Check the dimensions of the stored quantum system against `construct`
  /**
   * \tparam CONSTRUCT typically either shared_ptr to a quantumdata::StateVector (as in MCWF_Trajectory) or a quantumdata::DensityOperator (as in Master)
   */
  template<typename CONSTRUCT>
  void checkDimension(CONSTRUCT&& c) const {if (c!=*(this->getQS())) throw DimensionalityMismatchException("during QuantumTrajectory construction");}
  
private:
  mutable double t0_;

};



} // quantumtrajectory

#endif // CPPQEDCORE_QUANTUMTRAJECTORY_QUANTUMTRAJECTORY_H_INCLUDED
