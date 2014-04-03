// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION_TCC_INCLUDED
#define   CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION_TCC_INCLUDED

#include "Evolution_.h"

#include "EnsembleMCWF.tcc"
#include "Master.tcc"
#include "TimeAveragingMCWF_Trajectory.tcc"

#include "Trajectory.tcc"

#include <boost/make_shared.hpp>


template<typename V, int RANK>
const typename quantumdata::LazyDensityOperator<RANK>::Ptr
evolve(quantumdata::DensityOperator<RANK>& rho,
       typename structure::QuantumSystem<RANK>::Ptr sys,
       const evolution::Pars& pe)
{
  if (pe.evol==evolution::MASTER_FAST) {
    Master<RANK,V,true> traj(rho,sys,pe,pe.negativity);
    trajectory::run(traj,pe);
  }
  else {
    Master<RANK,V> traj(rho,sys,pe,pe.negativity);
    trajectory::run(traj,pe);
  }

  return boost::make_shared<quantumdata::DensityOperator<RANK> >(rho); // deep copy

}


template<typename V, int RANK>
const typename quantumdata::LazyDensityOperator<RANK>::Ptr
evolve(quantumdata::StateVector<RANK>& psi,
       typename structure::QuantumSystem<RANK>::Ptr sys,
       const evolution::Pars& pe)
{
  typedef quantumdata::StateVector    <RANK> SV;
  typedef quantumdata::DensityOperator<RANK> DO;

  if      (pe.evol==evolution::SINGLE) {
    trajectory::run(*makeMCWF(psi,sys,pe),pe);
    return boost::make_shared<SV>(psi); // deep copy
  }
  else if (pe.evol==evolution::ENSEMBLE) {
    EnsembleMCWF<RANK,V> traj(psi,sys,pe,pe.negativity);
    trajectory::run(traj,pe);
    return boost::make_shared<DO>(traj.toBeAveraged()); // deep copy
  }
  else {
    DO rho(psi);
    return evolve<V>(rho,sys,pe);
  }

}


template<int RANK, typename SYS>
const boost::shared_ptr<MCWF_Trajectory<RANK> > evolution::makeMCWF(quantumdata::StateVector<RANK>& psi, const SYS& sys, const evolution::Pars& pe)
{
  if (pe.timeAverage) return boost::make_shared<TimeAveragingMCWF_Trajectory<RANK> >(psi,sys,pe,pe.relaxationTime);
  else                return boost::make_shared<             MCWF_Trajectory<RANK> >(psi,sys,pe                  );
}


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION_TCC_INCLUDED
