#define BLITZ_ARRAY_LARGEST_RANK 15

#include "EvolutionComposite.h"

#include "JaynesCummings.h"
#include "QbitModeCorrelations.h"

using namespace std;

typedef quantumdata::StateVector<5> StateVector;


class JaynesCummings_ReducedDensityOperator : public jaynescummings::Base<true>, public structure::ElementAveraged<2,false>
{
public:
  JaynesCummings_ReducedDensityOperator(qbit::Ptr qbit, mode::Ptr mode, const dcomp& g)
    : jaynescummings::Base<true>(qbit,mode,g),
      structure::ElementAveraged<2,false>
};


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

  const mode::Ptr mode(mode::make(pplm ,qmp,structure::averaged::ReducedDensityOperator("Mode",pplm.cutoff,true)));

  const jaynescummings::Ptr jcCorr(jaynescummings::make<QbitModeCorrelations>(qbit,mode,pjc)), jc(jaynescummings::make(qbit,mode,pjc));
  
  StateVector psi(qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*mode::init(pplm));
  psi.renorm();

  evolve<tmptools::Vector<0> >(psi,composite::make(Act<0,4>(jcCorr),Act<1,4>(jc),Act<2,4>(jc),Act<3,4>(jc)),pe);
 



}
