// -*- C++ -*-
#ifndef   QUANTUMTRAJECTORY_IMPL_ENSEMBLEMCWF_TCC_INCLUDED
#define   QUANTUMTRAJECTORY_IMPL_ENSEMBLEMCWF_TCC_INCLUDED

#include "EnsembleMCWF.h"

#include "ParsMCWF_Trajectory.h"

namespace quantumtrajectory {


namespace ensemblemcwf {


template<int RANK>
std::auto_ptr<typename Base<RANK>::StateVectors>
Base<RANK>::stateVectors(const StateVector& psi, size_t nTraj)
{
  StateVectors res;
  for (size_t i=0; i<nTraj; ++i)
    res.push_back(new StateVector(psi));
  // Here trying to use lambda will be much less efficient
  
  return res.release();
}


template<int RANK>
std::auto_ptr<typename Base<RANK>::Trajectories>
Base<RANK>::trajectories(StateVectors& psis, QuantumSystemPtr qs, const ParsMCWF& p, const StateVectorLow& scaleAbs)
{
  Trajectories res;

  p.logLevel=(p.logLevel>0 ? 1 : p.logLevel); // reduced logging for individual trajectories in an Ensemble

  for (typename StateVectors::iterator i=psis.begin(); i!=psis.end(); (++i, ++p.seed))
    res.push_back(new MCWF_Trajectory<RANK>(*i,qs,p,scaleAbs));

  return res.release();
}


template<int RANK>
Base<RANK>::Base(
                 const StateVector& psi,
                 QuantumSystemPtr qs,
                 const ParsMCWF& p,
                 const StateVectorLow& scaleAbs
                 )
  : StateVectorsBase(stateVectors(psi,p.nTraj)),
    Ensemble(trajectories(StateVectorsBase::member,qs,p,scaleAbs),p.logLevel<0),
    rho_(psi.getDimensions(),false), qs_(qs)
{
  std::cerr<<p.logLevel<<std::endl;
}


} // ensemblemcwf

} // quantumtrajectory


#endif // QUANTUMTRAJECTORY_IMPL_ENSEMBLEMCWF_TCC_INCLUDED
