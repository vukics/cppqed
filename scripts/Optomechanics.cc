#include "EvolutionHigh.h"

#include "ModeCorrelations.h"
#include "NX_CoupledModes.h"

#include "BinarySystem.h"


using namespace std;

typedef quantumdata::StateVector<2> StateVector;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  ParsEvolution pe(p); // Driver Parameters
  mode::ParsPumpedLossy pA(p,"Opt");
  mode::ParsLossy       pB(p,"Mech"); 

  double& u=p.addTitle("Optomechanics","").add("u","N-Q coupling",1.);

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  mode::SmartPtr mA(mode::maker(pA,QMP_UIP,mode::DoNotAverage()));
  mode::SmartPtr mB(mode::maker(pB,QMP_UIP,mode::DoNotAverage()));

  NX_CoupledModes<ModeCorrelations> nx(mA,mB,-sqrt(2)*u);

  StateVector psi(mode::init(pA)*mode::init(pB));
  psi.renorm();

  evolve(psi,BinarySystem(nx),pe,tmptools::Vector<0>());


  } catch (ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}


}
