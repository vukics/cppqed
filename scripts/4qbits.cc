#include "EvolutionComposite.h"
#include "JaynesCummings.h"
#include "QbitModeCorrelations.h"

using namespace std;

typedef quantumdata::StateVector<5> StateVector;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters
  qbit::ParsPumpedLossy pplqb(p);
  qbit::ParsPumpedLossy pplqb2(p,"2");
  qbit::ParsPumpedLossy pplqb3(p,"3"); 
  qbit::ParsPumpedLossy pplqb4(p,"4");
  mode::ParsPumpedLossy pplm (p); 
  jaynescummings::Pars  pjc  (p);
  jaynescummings::Pars  pjc2  (p,"2");
  jaynescummings::Pars  pjc3  (p,"3");
  jaynescummings::Pars  pjc4  (p,"4"); 

  QM_Picture& qmp=p.add("picture","Quantum mechanical picture",QMP_IP);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  if ((pe.evol==EM_MASTER /* || pe.evol==EM_convergence */) && qmp==QMP_IP) qmp=QMP_UIP;

  qbit::Ptr qbit(qbit::make(pplqb,qmp));
  qbit::Ptr qbit2(qbit::make(pplqb2,qmp));
  qbit::Ptr qbit3(qbit::make(pplqb3,qmp));
  qbit::Ptr qbit4(qbit::make(pplqb4,qmp));

  mode::Ptr mode(mode::make(pplm ,qmp));

  StateVector psi(qbit::init(pplqb)*qbit::init(pplqb2)*qbit::init(pplqb3)*qbit::init(pplqb4)*mode::init(pplm));
  psi.renorm();

  evolve<tmptools::Vector<0> >(psi,composite::make(Act<0,4>(jaynescummings::make<QbitModeCorrelations>(qbit,mode,pjc)),
						   Act<1,4>(jaynescummings::make(qbit2,mode,pjc2)),
						   Act<2,4>(jaynescummings::make(qbit3,mode,pjc3)),
						   Act<3,4>(jaynescummings::make(qbit4,mode,pjc4))),pe);
 
  } catch (const cpputils::TaggedException& te) {cerr<<"Caught exception with tag: "<<te.getTag()<<endl; exit(1);}


}
