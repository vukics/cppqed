// -*- C++ -*-
#ifndef QUANTUMTRAJECTORY_QUANTUMTRAJECTORY_H_INCLUDED
#define QUANTUMTRAJECTORY_QUANTUMTRAJECTORY_H_INCLUDED

#include "QuantumSystem.h"
#include "Trajectory.h"


namespace quantumtrajectory {


template<int RANK>
inline double initialTimeStep(typename structure::QuantumSystem<RANK>::Ptr qs)
{
  return trajectory::initialTimeStep(qs->highestFrequency());
}


} // quantumtrajectory

#endif // QUANTUMTRAJECTORY_QUANTUMTRAJECTORY_H_INCLUDED
