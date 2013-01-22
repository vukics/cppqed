// -*- C++ -*-
#ifndef   QUANTUMTRAJECTORY_IMPL_EVOLUTION_TCC_INCLUDED
#define   QUANTUMTRAJECTORY_IMPL_EVOLUTION_TCC_INCLUDED

#include "Evolution_.h"

#include "ParsEvolution.h"

#include "impl/DO_Display.tcc"
#include "EnsembleMCWF.h"
#include "Master.h"
#include "impl/Trajectory.tcc"

#include <iostream>
#include <string>



using namespace quantumtrajectory;


template<typename V, int RANK>
void evolve(quantumdata::StateVector<RANK>& psi,
	    typename structure::QuantumSystem<RANK>::Ptr sys,
	    const ParsEvolution& pe)
{
  using namespace std;

  switch (pe.evol) {


  case EM_SINGLE: {

    MCWF_Trajectory<RANK>
      traj(psi,sys,pe);

    trajectory::evolve(traj,pe);

    break;

  }


  case EM_ENSEMBLE: {

    EnsembleMCWF<RANK,V>
      traj(psi,sys,pe,pe.negativity);

    if      (pe.Dt) runDt(traj,pe.T,pe.Dt,pe.displayInfo);
    else cout<<"Nonzero Dt required!"<<endl;

    break;

  }


  case EM_MASTER: {

    quantumdata::DensityOperator<RANK> rho(psi);

    Master<RANK,V>
      traj(rho,sys,pe,pe.negativity);

    trajectory::evolve(traj,pe);
    
    break;

  }


  case EM_MASTER_FAST: {

    quantumdata::DensityOperator<RANK> rho(psi);

    Master<RANK,V,true>
      traj(rho,sys,pe,pe.negativity);

    trajectory::evolve(traj,pe);

    break;

  }

    /*
  case EM_CONVERGENCE: {
    DensityOperator rho(psi);
    EnsembleMCWF EM(psi,sys,pe);
    Master M(rho,sys,pe);

    M.o()<<"# Comparing the results of"<<endl<<endl;
    EM.Parameters();
    M.o()<<endl<<"# With those of"<<endl<<endl;
    M.Parameters();

    M.o()<<"# Displaying in every "<<pe.Dt<<endl<<endl;

    M.Display();
    while (M()<pe.T) {
      M.evolve (min(pe.Dt,pe.T-M ())); 
      M.Display();
      EM.evolve(min(pe.Dt,pe.T-EM()));
      M.o()<<"# "<<FrobeniusNorm(EM.ToBeAveraged()-rho)<<endl;      
      }
  }
    */
  }
  
}



#endif // QUANTUMTRAJECTORY_IMPL_EVOLUTION_TCC_INCLUDED
