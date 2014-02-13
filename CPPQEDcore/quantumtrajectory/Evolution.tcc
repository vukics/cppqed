// -*- C++ -*-
#ifndef   QUANTUMTRAJECTORY_EVOLUTION_TCC_INCLUDED
#define   QUANTUMTRAJECTORY_EVOLUTION_TCC_INCLUDED

#include "Evolution_.h"

#include "DO_Display.tcc"
#include "EnsembleMCWF.h"
#include "Master.h"
#include "TimeAveragingMCWF_Trajectory.tcc"

#include "Trajectory.tcc"

#include <boost/make_shared.hpp>

#include <iostream>
#include <string>



template<typename V, int RANK>
void evolve(quantumdata::StateVector<RANK>& psi,
            typename structure::QuantumSystem<RANK>::Ptr sys,
            const evolution::Pars& pe)
{
  using namespace std; using namespace evolution;

  switch (pe.evol) {


  case SINGLE: {

    trajectory::run(*makeMCWF(psi,sys,pe),pe);

    break;

  }


  case ENSEMBLE: {

    EnsembleMCWF<RANK,V>
      traj(psi,sys,pe,pe.negativity);

    trajectory::run(traj,pe);

    break;

  }


  case MASTER: {

    quantumdata::DensityOperator<RANK> rho(psi);

    Master<RANK,V>
      traj(rho,sys,pe,pe.negativity);

    trajectory::run(traj,pe);
    
    break;

  }


  case MASTER_FAST: {

    quantumdata::DensityOperator<RANK> rho(psi);

    Master<RANK,V,true>
      traj(rho,sys,pe,pe.negativity);

    trajectory::run(traj,pe);

    break;

  }

  }
  
}


// For some reason, this does not compile:
template<int RANK, typename SYS>
const boost::shared_ptr<MCWF_Trajectory<RANK> > evolution::makeMCWF(quantumdata::StateVector<RANK>& psi, const SYS& sys, const evolution::Pars& pe)
{
  if (pe.timeAverage) return boost::make_shared<TimeAveragingMCWF_Trajectory<RANK> >(psi,sys,pe,pe.relaxationTime);
  else                return boost::make_shared<             MCWF_Trajectory<RANK> >(psi,sys,pe                  );
}


#endif // QUANTUMTRAJECTORY_EVOLUTION_TCC_INCLUDED
