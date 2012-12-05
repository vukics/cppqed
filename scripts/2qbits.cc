#include "EvolutionComposite.h"
#include "JaynesCummings.h"
#include "QbitModeCorrelations.h"

using namespace std;

typedef quantumdata::StateVector<3> StateVector;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters
  qbit::ParsPumpedLossy pplqb(p);
  qbit::ParsPumpedLossy pplqb2(p,"2");
  mode::ParsPumpedLossy pplm (p); 
  jaynescummings::Pars  pjc  (p);
  jaynescummings::Pars  pjc2  (p,"2");

  QM_Picture& qmp=p.add("picture","Quantum mechanical picture",QMP_IP);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  if ((pe.evol==EM_MASTER /* || pe.evol==EM_convergence */) && qmp==QMP_IP) qmp=QMP_UIP;

  qbit::Ptr qbit(qbit::make(pplqb,qmp));
  qbit::Ptr qbit2(qbit::make(pplqb2,qmp));

  mode::Ptr mode(mode::make(pplm ,qmp));

  StateVector psi(qbit::init(pplqb)*qbit::init(pplqb2)*mode::init(pplm));
  psi.renorm();

  evolve<tmptools::Vector<0> >(psi,composite::make(Act<0,2>(jaynescummings::make<QbitModeCorrelations>(qbit,mode,pjc)),
						   Act<1,2>(jaynescummings::make(qbit2,mode,pjc2))),pe);
 
  } catch (const cpputils::TaggedException& te) {cerr<<"Caught exception with tag: "<<te.getTag()<<endl; exit(1);}


}
