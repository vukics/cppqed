#include "Evolution.h"

#include "Mode.h"

int main(int argc, char* argv[])
{
  ParameterTable p;

  evolution::Pars pe(p);       // Driver parameters
  mode::ParsPumpedLossy pm(p); // Mode parameters

  pe.evol=evolution::MASTER;
  pm.cutoff=30;
  // ... other default values may follow

  update(p,argc,argv); // Parsing the command line

  // ****** ****** ****** ****** ****** ******

  mode::Ptr mode(make(pm,QMP_UIP));

  mode::StateVector psi(init(pm));

  evolve(psi,mode,pe);

}