// -*- C++ -*-
#ifndef   _EVOLUTION_MODE_IMPL_H
#define   _EVOLUTION_MODE_IMPL_H

#include "ParsEvolution.h"

#include "EnsembleMCWF.h"
#include "Master.h"

#include<iostream>
#include<string>



using namespace quantumtrajectory;


template<int RANK, typename V>
void evolve(quantumdata::StateVector<RANK>& psi, const structure::QuantumSystem<RANK>& sys,
	    const ParsEvolution& pe, V)
{
  using namespace std;

  switch (pe.evol) {


  case EM_SINGLE: {

    MCWF_Trajectory<RANK>
      traj(psi,sys,pe);
    
    if      (pe.dc) run  (traj,pe.T,pe.dc,pe.displayInfo);
    else if (pe.Dt) runDt(traj,pe.T,pe.Dt,pe.displayInfo);
    else cout<<"Nonzero dc OR Dt required!"<<endl;

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

    if      (pe.dc) run  (traj,pe.T,pe.dc,pe.displayInfo);
    else if (pe.Dt) runDt(traj,pe.T,pe.Dt,pe.displayInfo);
    else cout<<"Nonzero dc OR Dt required!"<<endl;

    break;

  }


  case EM_MASTER_FAST: {

    quantumdata::DensityOperator<RANK> rho(psi);

    Master<RANK,V,true>
      traj(rho,sys,pe,pe.negativity);

    if      (pe.dc) run  (traj,pe.T,pe.dc,pe.displayInfo);
    else if (pe.Dt) runDt(traj,pe.T,pe.Dt,pe.displayInfo);
    else cout<<"Nonzero dc OR Dt required!"<<endl;

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



#endif // _EVOLUTION_MODE_IMPL_H
