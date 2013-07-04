// #define BLITZ_ARRAY_LARGEST_RANK 15

#include "EvolutionComposite.h"

#include "JaynesCummings.h"
#include "QbitModeCorrelations.h"

using namespace std;

typedef quantumdata::StateVector<4> StateVector;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters
  qbit::ParsLossy pq0(p,"0"), pq1(p,"1"), pq2(p,"2");
  mode::ParsPumpedLossy pplm (p);
  jaynescummings::Pars  pjc0(p,"0"), pjc1(p,"1"), pjc2(p,"2");

  double& gamma_parallel=p.add("gamma_parallel","gamma_parallel",1.);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  pq2.gamma=(pq1.gamma=pq0.gamma);
  
  // ****** ****** ****** ****** ****** ******

  const LossyQbitWithPhaseNoise q0(pq0,gamma_parallel), q1(pq1,gamma_parallel), q2(pq2,gamma_parallel);

  const mode::Ptr mode(mode::make<mode::AveragedQuadratures>(pplm,QMP_IP));

  const jaynescummings::Ptr jc0(jaynescummings::make(q0,mode,pjc0)), jc1(jaynescummings::make(q1,mode,pjc1)),
                            jc2(jaynescummings::make<QbitModeCorrelations>(q2,mode,pjc2));
                            
  
  StateVector psi(qbit::init(pq0)*qbit::init(pq1)*qbit::init(pq2)*mode::init(pplm)
    // qbit::state0()*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*mode::fock(0,pplm.cutoff)+qbit::state1()*qbit::init(pplqb)*qbit::init(pplqb)*qbit::init(pplqb)*mode::fock(1,pplm.cutoff)
  );
  psi.renorm();

  evolve(psi,composite::make(Act<0,3>(jc0),
                             Act<1,3>(jc1),
                             Act<2,3>(jc2)),
         pe);
 



}
