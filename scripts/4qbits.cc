// #define BLITZ_ARRAY_LARGEST_RANK 15

#include "EvolutionComposite.h"

#include "AveragingUtils.h"
#include "JaynesCummings.h"
#include "QbitModeCorrelations.h"

using namespace std;

typedef quantumdata::StateVector<5> StateVector;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

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

  const mode::Ptr mode(mode::make<mode::AveragedQuadratures>(pplm ,qmp));

  const jaynescummings::Ptr jcCorr(jaynescummings::make<ReducedDensityOperator<2> >(qbit,mode,pjc,"JaynesCummings0-4",pplm.cutoff,true)), jc(jaynescummings::make(qbit,mode,pjc));
  
  StateVector psi(qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*mode::init(pplm));
  psi.renorm();

  evolve(psi,composite::make(Act<0,4>(jcCorr),
                             Act<1,4>(jc),
                             Act<2,4>(jc),
                             Act<3,4>(jc)),
         pe);
 



}
