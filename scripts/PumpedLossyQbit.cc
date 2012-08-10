#include "Evolution.h"
#include "MCWF.h"

#include "ParsQbit.h"
#include "Qbit.h"

using namespace qbit;
using namespace std;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters
  ParsPumpedLossy pplqb(p); 

  QM_Picture& qmp=p.add("picture","Quantum mechanical picture",QMP_IP);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  if (pe.evol==EM_MASTER && qmp==QMP_IP) qmp=QMP_UIP;

  SmartPtr qbit(make(pplqb,qmp));

  StateVector psi(init(pplqb));
  /*
  StateVector psi(state0()()+state1()());
  psi.renorm();
  */
  evolve(psi,qbit,pe);

  } catch (const ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}


}

