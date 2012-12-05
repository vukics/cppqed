#define BLITZ_ARRAY_LARGEST_RANK 19

#include "EvolutionComposite.h"
#include "JaynesCummings.h"
#include "QbitModeCorrelations.h"

using namespace std;

typedef quantumdata::StateVector<9> StateVector;


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
  qbit::ParsPumpedLossy pplqb5(p,"5");
  qbit::ParsPumpedLossy pplqb6(p,"6");
  qbit::ParsPumpedLossy pplqb7(p,"7");
  qbit::ParsPumpedLossy pplqb8(p,"8"); 
  mode::ParsPumpedLossy pplm (p); 
  jaynescummings::Pars  pjc  (p);
  jaynescummings::Pars  pjc2  (p,"2");
  jaynescummings::Pars  pjc3  (p,"3");
  jaynescummings::Pars  pjc4  (p,"4"); 
  jaynescummings::Pars  pjc5  (p,"5");
  jaynescummings::Pars  pjc6  (p,"6");
  jaynescummings::Pars  pjc7  (p,"7");
  jaynescummings::Pars  pjc8  (p,"8");

  QM_Picture& qmp=p.add("picture","Quantum mechanical picture",QMP_IP);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  if ((pe.evol==EM_MASTER /* || pe.evol==EM_convergence */) && qmp==QMP_IP) qmp=QMP_UIP;

  qbit::Ptr qbit(qbit::make(pplqb,qmp));
  qbit::Ptr qbit2(qbit::make(pplqb2,qmp));
  qbit::Ptr qbit3(qbit::make(pplqb3,qmp));
  qbit::Ptr qbit4(qbit::make(pplqb4,qmp));
  qbit::Ptr qbit5(qbit::make(pplqb5,qmp));
  qbit::Ptr qbit6(qbit::make(pplqb6,qmp));
  qbit::Ptr qbit7(qbit::make(pplqb7,qmp));
  qbit::Ptr qbit8(qbit::make(pplqb8,qmp));

  mode::Ptr mode(mode::make(pplm ,qmp));

  StateVector psi(qbit::init(pplqb)*qbit::init(pplqb2)*qbit::init(pplqb3)*qbit::init(pplqb4)*qbit::init(pplqb5)*qbit::init(pplqb6)*qbit::init(pplqb7)*qbit::init(pplqb8)*mode::init(pplm));
  psi.renorm();

  evolve<tmptools::Vector<0> >(psi,composite::make(Act<0,8>(jaynescummings::make<QbitModeCorrelations>(qbit,mode,pjc)),
						   Act<1,8>(jaynescummings::make(qbit2,mode,pjc2)),
						   Act<2,8>(jaynescummings::make(qbit3,mode,pjc3)),
						   Act<3,8>(jaynescummings::make(qbit4,mode,pjc4)),
						   Act<4,8>(jaynescummings::make(qbit5,mode,pjc5)),
						   Act<5,8>(jaynescummings::make(qbit6,mode,pjc6)),
						   Act<6,8>(jaynescummings::make(qbit7,mode,pjc7)),
						   Act<7,8>(jaynescummings::make(qbit8,mode,pjc8))),pe);
 
  } catch (const cpputils::TaggedException& te) {cerr<<"Caught exception with tag: "<<te.getTag()<<endl; exit(1);}


}
