// -*- C++ -*-
#ifndef   _ENSEMBLE_OF_MCWF_TRAJECTORIES_IMPL_H
#define   _ENSEMBLE_OF_MCWF_TRAJECTORIES_IMPL_H

#include "EnsembleMCWF.h"

namespace quantumtrajectory {


namespace ensemblemcwf {

namespace details {


template<int RANK>
std::auto_ptr<typename Base<RANK>::StateVectors>
stateVectors(const typename Base<RANK>::StateVector& psi, size_t nTraj)
{
  typedef typename Base<RANK>::StateVector  StateVector ;
  typedef typename Base<RANK>::StateVectors StateVectors;

  StateVectors res;
  for (size_t i=0; i<nTraj; i++)
    res.push_back(new StateVector(psi));
  // Here trying to use lambda will be much less efficient
  
  return res.release();
}
										      

template<int RANK>
std::auto_ptr<typename Base<RANK>::Trajectories>
trajectories(
	     typename Base<RANK>::StateVectors& psis,
	     const typename Base<RANK>::QuantumSystem& sys,
	     const ParsMCWF_Trajectory& p, 
	     const typename Base<RANK>::StateVectorLow& scaleAbs
	     )
{
  typedef typename Base<RANK>::Trajectories Trajectories;
  typedef typename Base<RANK>::StateVectors StateVectors;

  Trajectories res;

  p.logLevel=(p.logLevel>0 ? 1 : 0); // reduced logging for individual trajectories in an Ensemble

  typename StateVectors::iterator i=psis.begin();
  for (size_t j=0; j<p.nTraj; (++i, j++, p.seed++))
    res.push_back(new MCWF_Trajectory<RANK>(*i,sys,p,scaleAbs));

  return res.release();
}


} // details


template<int RANK>
Base<RANK>::Base(
		 const StateVector& psi,
		 const QuantumSystem& sys,
		 const ParsMCWF_Trajectory& p,
		 const StateVectorLow& scaleAbs
		 )
  : StateVectorsBase(details::stateVectors<RANK>(psi,p.nTraj)),
    EnsembleTrajectories(details::trajectories<RANK>(StateVectorsBase::member,sys,p,scaleAbs),p.logLevel<0),
    rho_(psi.getDimensions(),false)
{
}


} // ensemblemcwf

} // quantumtrajectory


#endif // _ENSEMBLE_OF_MCWF_TRAJECTORIES_IMPL_H
