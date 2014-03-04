#include "Evolution.h"

#include "Qbit.h"

using namespace qbit;
using namespace std;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  
  ParameterTable p;

  evolution::Pars pe(p); // Driver Parameters
  ParsPumpedLossy pplqb(p); 

  // Parameter finalization
  QM_Picture& qmp=updateWithPicture(p,argc,argv);
  
  // ****** ****** ****** ****** ****** ******

  Ptr qbit(make(pplqb,qmp));

  StateVector psi(init(pplqb));
  /*
  StateVector psi(state0()()+state1()());
  psi.renorm();
  */
  evolve(psi,qbit,pe);




}

