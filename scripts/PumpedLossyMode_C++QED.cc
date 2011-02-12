#include "EvolutionHigh.h"

#include "Mode.h"
#include "ParsMode.h"

using namespace std;
using namespace mode;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters
  ParsPumpedLossy pplm(p); 

  QM_Picture& qmp=p.add("picture","Quantum mechanical picture",QMP_IP);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  if (pe.evol==EM_MASTER && qmp==QMP_IP) qmp=QMP_UIP;

  SmartPtr mode(maker(pplm,qmp,AveragedQuadratures()));

  StateVector psi(mode::init(pplm));

  evolve(psi,*mode,pe);

  } catch (ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}


}

