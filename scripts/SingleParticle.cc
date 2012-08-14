#include "Evolution.h"

#include "Particle.h"
#include "ParticleInitialCondition.h"

using namespace std;
using namespace particle;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  ParsEvolution pe   (p); // Driver Parameters
  ParsPumped    ppart(p); 

  QM_Picture& qmp=p.add("picture","Quantum mechanical picture",QMP_IP);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  // if (pe.evol==EM_MASTER && qmp==QMP_IP) qmp=QMP_UIP;

  SmartPtr part(make(ppart,qmp));

  if (!ppart.init.getSig() && !ppart.vClass) {cerr<<"Incorrect initial condition"<<endl; abort();}

  StateVector psi(init(ppart));

  evolve(psi,part,pe);

  } catch (const ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}


}

