#include "Evolution.h"
#include "Mode.h"

#include "ModeCorrelations.h"
#include "NX_CoupledModes.h"

#include "BinarySystem.h"


using namespace std;
using namespace mode;

typedef quantumdata::StateVector<2> StateVector2;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters
  ParsPumpedLossy pA(p,"Opt");
  ParsLossy       pB(p,"Mech"); 

  double& u=p.addTitle("Optomechanics","").add("u","N-Q coupling",1.);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  StateVector2 psi(init(pA)*init(pB));
  psi.renorm();

  evolve<tmptools::Vector<0> >
    (psi,
     binary::make(NX_CoupledModes<ModeCorrelations>(make<DoNotAverage>(pA,QMP_UIP),
						    make<DoNotAverage>(pB,QMP_UIP),
						    -sqrt(2)*u)),
     pe);

  } catch (const ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}


}
