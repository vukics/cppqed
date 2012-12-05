// -*- C++ -*-
#ifndef   QUANTUMTRAJECTORY_IMPL_ENSEMBLEMCWF_TCC_INCLUDED
#define   QUANTUMTRAJECTORY_IMPL_ENSEMBLEMCWF_TCC_INCLUDED

#include "EnsembleMCWF.h"

#include "ParsMCWF_Trajectory.h"

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
	     typename Base<RANK>::QuantumSystemPtr qs,
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
    res.push_back(new MCWF_Trajectory<RANK>(*i,qs,p,scaleAbs));

  return res.release();
}


} // details


template<int RANK>
Base<RANK>::Base(
		 const StateVector& psi,
		 QuantumSystemPtr qs,
		 const ParsMCWF_Trajectory& p,
		 const StateVectorLow& scaleAbs
		 )
  : StateVectorsBase(details::stateVectors<RANK>(psi,p.nTraj)),
    Ensemble(details::trajectories<RANK>(StateVectorsBase::member,qs,p,scaleAbs),p.logLevel<0),
    rho_(psi.getDimensions(),false), qs_(qs)
{
}


} // ensemblemcwf

} // quantumtrajectory


#endif // QUANTUMTRAJECTORY_IMPL_ENSEMBLEMCWF_TCC_INCLUDED
