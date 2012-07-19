#include "EvolutionHigh.h"

#include "Mode.h"

#include <boost/assign/list_of.hpp>

using namespace std ;
using namespace mode;



int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters
  ParsPumpedLossy pplm(p); 

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  pe.nTraj=1000/*00*/;
  const list<size_t> sections=boost::assign::list_of(100)(200)(500)(1000)/*(2000)(5000)(10000)(20000)(50000)(100000)*/;

  SmartPtr modeMaster(make(pplm,QMP_UIP)), mode(make(pplm,QMP_IP));
  StateVector psi(init(pplm));

  DensityOperator rho(psi);

  double epsRelOld=1e-12;

  swap(epsRelOld,pe.epsRel);
  Master      <1> m(rho,*modeMaster,pe,false);
  m.displayParameters();

  swap(epsRelOld,pe.epsRel);
  EnsembleMCWF<1> e(psi,*mode,pe,false);
  e.displayParameters();

  while (m.getTime()<pe.T) {
    {
      double timestep=min(pe.Dt,pe.T-m.getTime());
      m.evolve(timestep);
      cout<<m.getTime()<<'\t';
      e.evolve(timestep);
    }
    for (list<size_t>::const_iterator i=sections.begin(); i!=sections.end(); ++i) {
      double avr=0;
      for (size_t begin=0; begin<pe.nTraj; begin+=*i)
	avr+=frobeniusNorm(rho-e.toBeAveraged(begin,*i));
      cout<<avr*(double(*i)/pe.nTraj)<<'\t';
    }

  }

  } catch (const ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}


}



// bin/scripts/gcc-4.6/release/PumpedLossyModeConvergence --nTh 10 --cutoff 100 --kappa 5 --eta "(-10,9)" --T 0.5 --Dt 0.0005 --precision 6 --minitFock 20 --seed 1000 --logLevel 1 --logLevel -1
