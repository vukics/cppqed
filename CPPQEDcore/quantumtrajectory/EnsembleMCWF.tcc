// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_QUANTUMTRAJECTORY_ENSEMBLEMCWF_TCC_INCLUDED
#define   CPPQEDCORE_QUANTUMTRAJECTORY_ENSEMBLEMCWF_TCC_INCLUDED

#include "EnsembleMCWF.h"

#include "DO_Display.tcc"
#include "MCWF_Trajectory.tcc"
#include "ParsMCWF_Trajectory.h"



template<int RANK>
auto
quantumtrajectory::ensemble::Base<RANK>::stateVectors(const StateVector& psi, size_t nTraj) -> std::auto_ptr<StateVectors>
{
  StateVectors res;
  for (size_t i=0; i<nTraj; ++i)
    res.push_back(new StateVector(psi));
  // Here trying to use lambda will be much less efficient
  
  return res.release();
}


template<int RANK>
auto
quantumtrajectory::ensemble::Base<RANK>::trajectories(StateVectors& psis, QuantumSystemPtr qs, const Pars& p, const StateVectorLow& scaleAbs) -> std::auto_ptr<Trajectories>
{
  Trajectories res;

  p.logLevel=(p.logLevel>0 ? 1 : p.logLevel); // reduced logging for individual trajectories in an Ensemble

  for (auto i=psis.begin(); i!=psis.end(); (++i, ++p.seed))
    res.push_back(new MCWF_Trajectory<RANK>(*i,qs,p,scaleAbs));

  return res.release();
}


template<int RANK>
quantumtrajectory::ensemble::Base<RANK>::Base(
                                              const StateVector& psi,
                                              QuantumSystemPtr qs,
                                              const Pars& p,
                                              const StateVectorLow& scaleAbs
                                              )
  : StateVectorsBase(stateVectors(psi,p.nTraj)),
    Ensemble(trajectories(StateVectorsBase::member,qs,p,scaleAbs),p.logLevel<0),
    rho_(psi.getDimensions(),false), qs_(qs), nBins_(p.nBins), nJumpsPerBin_(p.nJumpsPerBin)
{
}


template<int RANK>
std::ostream&
quantumtrajectory::ensemble::Base<RANK>::logOnEnd_v(std::ostream& os) const
{
  LoggerList loggerList;
  for (auto& i : this->getTrajectories())
    if (const auto traj=dynamic_cast<const MCWF_Trajectory<RANK>*>(&i))
      loggerList.push_back(traj->getLogger());
  
  return displayLog(os,loggerList,nBins_,nJumpsPerBin_);
}



#endif // CPPQEDCORE_QUANTUMTRAJECTORY_ENSEMBLEMCWF_TCC_INCLUDED
