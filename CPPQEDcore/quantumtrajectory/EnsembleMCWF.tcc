// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDCORE_QUANTUMTRAJECTORY_ENSEMBLEMCWF_TCC_INCLUDED
#define   CPPQEDCORE_QUANTUMTRAJECTORY_ENSEMBLEMCWF_TCC_INCLUDED

#include "EnsembleMCWF.h"

#include "DO_Display.tcc"
#include "MCWF_Trajectory.tcc"
#include "ParsMCWF_Trajectory.h"


template<int RANK>
auto
quantumtrajectory::ensemble::Base<RANK>::trajectories(std::shared_ptr<const StateVector> psi, size_t nTraj, QuantumSystemPtr qs, const Pars& p, const StateVectorLow& scaleAbs) -> std::unique_ptr<Trajectories>
{
  Trajectories res;

  p.logLevel=(p.logLevel>0 ? 1 : p.logLevel); // reduced logging for individual trajectories in an Ensemble

  for (size_t i=0; i<nTraj; (++i, ++p.seed) ) 
    res.push_back(new MCWF_Trajectory<RANK>(std::make_shared<StateVector>(*psi),qs,p,scaleAbs));

  return res.release();
}


template<int RANK>
quantumtrajectory::ensemble::Base<RANK>::Base(
                                              std::shared_ptr<const StateVector> psi,
                                              QuantumSystemPtr qs,
                                              const Pars& p,
                                              const StateVectorLow& scaleAbs
                                              )
  : Ensemble(trajectories(psi,p.nTraj,qs,p,scaleAbs),p.logLevel<0),
    qs_(qs), nBins_(p.nBins), nJumpsPerBin_(p.nJumpsPerBin)
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
