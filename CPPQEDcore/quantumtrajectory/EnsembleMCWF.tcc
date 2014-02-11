// -*- C++ -*-
#ifndef   QUANTUMTRAJECTORY_ENSEMBLEMCWF_TCC_INCLUDED
#define   QUANTUMTRAJECTORY_ENSEMBLEMCWF_TCC_INCLUDED

#include "EnsembleMCWF.h"

#include "MCWF_Trajectory.tcc"
#include "ParsMCWF_Trajectory.h"



template<int RANK>
auto
quantumtrajectory::mcwf::ensemble::Base<RANK>::stateVectors(const StateVector& psi, size_t nTraj) -> std::auto_ptr<StateVectors>
{
  StateVectors res;
  for (size_t i=0; i<nTraj; ++i)
    res.push_back(new StateVector(psi));
  // Here trying to use lambda will be much less efficient
  
  return res.release();
}


template<int RANK>
auto
quantumtrajectory::mcwf::ensemble::Base<RANK>::trajectories(StateVectors& psis, QuantumSystemPtr qs, const Pars& p, const StateVectorLow& scaleAbs) -> std::auto_ptr<Trajectories>
{
  Trajectories res;

  p.logLevel=(p.logLevel>0 ? 1 : p.logLevel); // reduced logging for individual trajectories in an Ensemble

  for (auto i=psis.begin(); i!=psis.end(); (++i, ++p.seed))
    res.push_back(new MCWF_Trajectory<RANK>(*i,qs,p,scaleAbs));

  return res.release();
}


template<int RANK>
quantumtrajectory::mcwf::ensemble::Base<RANK>::Base(
                                                    const StateVector& psi,
                                                    QuantumSystemPtr qs,
                                                    const Pars& p,
                                                    const StateVectorLow& scaleAbs
                                                    )
  : StateVectorsBase(stateVectors(psi,p.nTraj)),
    Ensemble(trajectories(StateVectorsBase::member,qs,p,scaleAbs),p.logLevel<0),
    rho_(psi.getDimensions(),false), qs_(qs)
{
}


template<int RANK>
std::ostream&
quantumtrajectory::mcwf::ensemble::Base<RANK>::logOnEnd_v(std::ostream& os) const
{
  LoggerList loggerList;
  for (typename Trajectories::const_iterator i=getTrajectories().begin(); i!=getTrajectories().end(); ++i)
    if (const auto traj=dynamic_cast<const MCWF_Trajectory<RANK>*>(&(*i)))
      loggerList.push_back(traj->getLogger());
  
  return displayLog(os,loggerList);
}



#endif // QUANTUMTRAJECTORY_ENSEMBLEMCWF_TCC_INCLUDED
