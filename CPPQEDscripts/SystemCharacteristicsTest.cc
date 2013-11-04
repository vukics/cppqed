#include "EvolutionComposite.h"

#include "JaynesCummings.h"
#include "NX_CoupledModes.h"


using namespace std;

typedef quantumdata::StateVector<3> StateVector;

typedef composite::result_of::Make<Act<0,2>,Act<1,2> >::type System;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  
  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters

  mode::ParsPumpedLossy pplqb1(p,"1");
  mode::ParsPumpedLossy pplqb2(p,"2");

  mode::ParsPumpedLossy pplm  (p);

  double &u1=p.addTitle("Optomechanics1","").add("u1","N-Q coupling",1.), &u2=p.addTitle("Optomechanics2","").add("u2","N-Q coupling",1.);

  /*
  jaynescummings::Pars  pjc1  (p,"1");
  jaynescummings::Pars  pjc2  (p,"2");
  */
  
  QM_Picture& qmp=p.add("picture","Quantum mechanical picture",QMP_IP);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  if ((pe.evol==EM_MASTER /* || pe.evol==EM_convergence */) && qmp==QMP_IP) qmp=QMP_UIP;

  mode::Ptr mode1(mode::make(pplqb1,qmp)), mode2(mode::make(pplqb2,qmp));
  mode::Ptr mode (mode::make(pplm  ,qmp));

  nxcoupledmodes::Ptr jc1(nxcoupledmodes::make(mode1,mode,u1)), jc2(nxcoupledmodes::make(mode2,mode,u2));

  // BinarySystem system(jc);

  StateVector psi(mode::init(pplqb1)*mode::init(pplqb2)*mode::init(pplm)), psiv(psi);

  System system(composite::make(Act<0,2>(jc1),Act<1,2>(jc2)));
  
  structure::QuantumSystemWrapper<3,true> qs(system);
  
  system->displayParameters(cout);
  
  qs.addContribution(0,psi(),psiv(),0);
  
  qs.displayCharacteristics(cout)<<endl;
  
/*  

evolve<tmptools::Vector<0> >(psi,binary::make(jc),pe);

*/




}
