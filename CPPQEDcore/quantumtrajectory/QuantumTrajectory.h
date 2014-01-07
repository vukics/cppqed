/// \briefFile{Defines quantumtrajectory::initialTimeStep()}
// -*- C++ -*-
#ifndef QUANTUMTRAJECTORY_QUANTUMTRAJECTORY_H_INCLUDED
#define QUANTUMTRAJECTORY_QUANTUMTRAJECTORY_H_INCLUDED

#include "QuantumSystem.h"
#include "Structure.h"

#include "Trajectory.h"


namespace quantumtrajectory {


template<int RANK>
inline double initialTimeStep(typename structure::QuantumSystem<RANK>::Ptr qs)
{
  return trajectory::initialTimeStep(qs->highestFrequency());
}

/// Class hosting common code of MCWF_Trajectory and Master.
/**
 * It contains the structure::QuantumSystemWrapper and a Adaptive::readStateMore_v implementation
 * which sets tIntPic0_. Because this class has to be inserted in the inheritence
 * chain of both MCWF_Trajectory and Master, the base class of QuantumTrajectory is passed as a
 * template parameter BASE and the constructor uses perfect forwarding semantics.
 */
template<int RANK, typename BASE>
class QuantumTrajectory : public BASE
{
protected:
  typedef structure::QuantumSystemWrapper<RANK,true> QuantumSystemWrapper;

  template<typename... ArgumentPack>
  QuantumTrajectory(typename structure::QuantumSystem<RANK>::Ptr qs, bool isNoisy, ArgumentPack&&... argumentPack) 
    : BASE(std::forward<ArgumentPack>(argumentPack)...), tIntPic0_(0), qs_(qs, isNoisy) {};

  const QuantumSystemWrapper getQSW() const {return qs_;}

  cpputils::iarchive&  readStateMore_v(cpputils::iarchive &iar)
    { BASE::readStateMore_v(iar); if (getQSW().getEx()) tIntPic0_=BASE::getTime(); return iar; }

  mutable double tIntPic0_ ; // The time instant of the beginning of the current time step.

private:
  const QuantumSystemWrapper qs_;

};



} // quantumtrajectory

#endif // QUANTUMTRAJECTORY_QUANTUMTRAJECTORY_H_INCLUDED
