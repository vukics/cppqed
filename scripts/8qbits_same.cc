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
  mode::ParsPumpedLossy pplm (p); 
  jaynescummings::Pars  pjc  (p);

  QM_Picture& qmp=p.add("picture","Quantum mechanical picture",QMP_IP);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  if ((pe.evol==EM_MASTER /* || pe.evol==EM_convergence */) && qmp==QMP_IP) qmp=QMP_UIP;

  const qbit::Ptr qbit(qbit::make(pplqb,qmp));

  const mode::Ptr mode(mode::make(pplm ,qmp));

  const jaynescummings::Ptr jcCorr(jaynescummings::make<QbitModeCorrelations>(qbit,mode,pjc)), jc(jaynescummings::make(qbit,mode,pjc));
  
  StateVector psi(qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*mode::init(pplm));
  psi.renorm();

  evolve<tmptools::Vector<0> >(psi,composite::make(Act<0,8>(jcCorr),
						   Act<1,8>(jc),
						   Act<2,8>(jc),
						   Act<3,8>(jc),
						   Act<4,8>(jc),
						   Act<5,8>(jc),
						   Act<6,8>(jc),
						   Act<7,8>(jc)),pe);
 
  } catch (const cpputils::TaggedException& te) {cerr<<"Caught exception with tag: "<<te.getTag()<<endl; exit(1);}


}
