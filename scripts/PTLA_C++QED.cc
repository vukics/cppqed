#include "Evolution.h"
#include "MCWF.h"

#include "ParsQbit.h"
#include "PumpedTwoLevelAtom.h"
#include "Qbit.h"

using namespace std;
using namespace qbit;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters
  ParsPumpedLossy pp2la(p); 

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  PumpedTwoLevelAtomSch atom(pp2la);

  StateVector psi(init(pp2la));

  evolve(psi,atom,pe);

  } catch (const ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}


}

