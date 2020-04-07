// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION_TCC_INCLUDED
#define   CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION_TCC_INCLUDED

#include "Evolution_.h"

#include "EnsembleMCWF.tcc"
#include "Master.tcc"
#include "TimeAveragingMCWF_Trajectory.tcc"

#include "Trajectory.tcc"


template<typename V, int RANK>
const typename quantumdata::LazyDensityOperator<RANK>::Ptr
evolve(std::shared_ptr<quantumdata::DensityOperator<RANK>> rho,
       typename structure::QuantumSystem<RANK>::Ptr sys,
       const evolution::Pars<>& pe)
{
  Master<RANK,V> traj(rho,sys,pe,pe.negativity);
  trajectory::run(traj,pe);

  return std::make_shared<quantumdata::DensityOperator<RANK> >(*rho); // deep copy

}


template<typename V, int RANK>
const typename quantumdata::LazyDensityOperator<RANK>::Ptr
evolve(std::shared_ptr<quantumdata::StateVector<RANK>> psi,
       typename structure::QuantumSystem<RANK>::Ptr sys,
       const evolution::Pars<>& pe)
{
  typedef quantumdata::StateVector    <RANK> SV;
  typedef quantumdata::DensityOperator<RANK> DO;

  if      (pe.evol==evolution::SINGLE) {
    trajectory::run(*makeMCWF(psi,sys,pe),pe);
    return std::make_shared<SV>(*psi); // deep copy
  }
  else if (pe.evol==evolution::ENSEMBLE) {
    EnsembleMCWF<RANK,V> traj(psi,sys,pe,pe.negativity);
    trajectory::run(traj,pe);
    return traj.averaged();
  }
  else {
    return evolve<V>(std::make_shared<DO>(*psi),sys,pe);
  }

}


template<int RANK>
const std::shared_ptr<MCWF_Trajectory<RANK> > evolution::makeMCWF(std::shared_ptr<quantumdata::StateVector<RANK>> psi, 
                                                                  typename structure::QuantumSystem<RANK>::Ptr sys, 
                                                                  const evolution::Pars<>& pe)
{
  if (pe.timeAverage) return std::make_shared<TimeAveragingMCWF_Trajectory<RANK> >(psi,sys,pe,pe.relaxationTime);
  else                return std::make_shared<             MCWF_Trajectory<RANK> >(psi,sys,pe                  );
}


template<typename Base>
evolution::Pars<Base>::Pars(parameters::ParameterTable& p, const std::string& mod) 
  : ParsRun(p,mod),
    Base(p,mod),
    evol(p.addTitle("Evolution",mod).addMod("evol",mod,"Evolution mode (single, ensemble, master)",SINGLE)),
    negativity(p.addMod("negativity",mod,"Calculates negativity in ensemble & master",false)),
    timeAverage(p.addMod("timeAverage",mod,"Calculates time averages in MCWF trajectory",false)),
    relaxationTime(p.addMod("relaxationTime",mod,"Relaxation time for time averaging",0.))
{}


#endif // CPPQEDCORE_QUANTUMTRAJECTORY_EVOLUTION_TCC_INCLUDED
