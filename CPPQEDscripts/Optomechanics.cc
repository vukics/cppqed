#include "EvolutionBinary.h"
#include "Mode.h"

#include "ModeCorrelations.h"
#include "NX_CoupledModes.h"


using namespace std;
using namespace mode;

typedef quantumdata::StateVector<2> StateVector2;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters
  ParsPumpedLossy pA(p,"Opt");
  ParsLossy       pB(p,"Mech"); 

  double& u=p.addTitle("Optomechanics","").add("u","N-Q coupling",1.);

  QM_Picture& qmp=p.add("picture","QM_Picture for particle (IP=UIP or Sch) and field (IP, UIP, or Sch)",QMP_UIP);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  StateVector2 psi(init(pA)*init(pB));
  psi.renorm();

  evolve<0>
    (psi,
     binary::make(nxcoupledmodes::make<ModeCorrelations>(make<DoNotAverage>(pA,qmp),
							 make<DoNotAverage>(pB,qmp),
							 -sqrt(2)*u)),
     pe);


  // Another way which works:
  /*
  evolve<0>
    (psi,
     binary::make(nxcoupledmodes::make<ModeCorrelations>(Mode<DoNotAverage>(pA),
							 Mode<DoNotAverage>(pB),
							 -sqrt(2)*u)),
     pe);
  */


  // A misuse:
  /*
  NX_CoupledModes<ModeCorrelations> nx(Mode<DoNotAverage>(pA),Mode<DoNotAverage>(pB),-sqrt(2)*u);
  // Here, it doesn't matter whether the interaction is created with or without its maker.

  evolve<0>
    (psi,
     binary::make(nx),
     pe);
  */

  // But this is fine again:
  /*
  Mode<DoNotAverage> m1(pA), m2(pB);

  NX_CoupledModes<ModeCorrelations> nx(m1,m2,-sqrt(2)*u);

  evolve<0>
    (psi,
     binary::make(nx),
     pe);
  */




}
