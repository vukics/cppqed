// Copyright András Vukics 2006–2015. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
// -*- C++ -*-
#ifndef CPPQEDCORE_QUANTUMTRAJECTORY_QUANTUMTRAJECTORY_H_INCLUDED
#define CPPQEDCORE_QUANTUMTRAJECTORY_QUANTUMTRAJECTORY_H_INCLUDED

#include "QuantumSystem.h"
#include "Structure.h"

#include "Trajectory.h"


/// Comprises modules representing trajectory drivers for simulating quantum systems
namespace quantumtrajectory {


/// Forwards to trajectory::initialTimeStep, with the highest frequency of the system taken as structure::QuantumSystem::highestFrequency
template<int RANK>
inline double initialTimeStep(typename structure::QuantumSystem<RANK>::Ptr qs)
{
  return trajectory::initialTimeStep(qs->highestFrequency());
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
class QuantumTrajectory : public BASE
{
protected:
  typedef structure::QuantumSystemWrapper<RANK,true> QuantumSystemWrapper;

  /// Constructor forwarding to `BASE` and QuantumSystemWrapper
  template<typename... ArgumentPack>
  QuantumTrajectory(typename structure::QuantumSystem<RANK>::Ptr qs, bool isNoisy, ArgumentPack&&... argumentPack) 
    : BASE(std::forward<ArgumentPack>(argumentPack)...), t0_(0), qs_(qs, isNoisy) {};

  const QuantumSystemWrapper getQSW() const {return qs_;}

  /// Forwards to `BASE`, but also sets \link getT0 `t0`\endlink
  cpputils::iarchive&  readStateMore_v(cpputils::iarchive &iar) override
    { BASE::readStateMore_v(iar); if (getQSW().getEx()) setT0(); return iar; }

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

  /// Check the dimensions of the stored quantum system against `construct` \tparam CONSTRUCT typically either a quantumdata::StateVector (as in MCWF_Trajectory) or a quantumdata::DensityOperator (as in Master)
  template<typename CONSTRUCT>
  void checkDimension(const CONSTRUCT& construct) const {if (construct!=*qs_.getQS()) throw DimensionalityMismatchException();}
  
private:
  mutable double t0_;
  const QuantumSystemWrapper qs_;

};



} // quantumtrajectory

#endif // CPPQEDCORE_QUANTUMTRAJECTORY_QUANTUMTRAJECTORY_H_INCLUDED
