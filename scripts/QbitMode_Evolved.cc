#include "ParsEvolution.h"
#include "SimulatedHigh.h"
#include "Mode.h"

#include "ParsQbit.h"
#include "Qbit.h"

#include "JaynesCummings.h"
#include "BinarySystem.h"


using namespace std       ;
using namespace cpputils  ;
using namespace trajectory;

typedef TTD_CARRAY(1) Array;


void derivs(double, const Array& b, Array& dbdt, 
	    const mode::ParsPumpedLossy& plm, const qbit::ParsPumpedLossy& pqb, const dcomp& g)
{
  dbdt(0)=dcomp(-plm.kappa,plm.delta)*b(0)+g*b(1)+plm.eta;
  dbdt(1)=dcomp(-pqb.gamma,pqb.delta)*b(1)-g*b(0)-pqb.eta;
}


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  try {

  ParameterTable p;

  ParsEvolution pt(p);

  qbit::ParsPumpedLossy pplqb(p); 
  mode::ParsPumpedLossy pplm (p); 
  jaynescummings::Pars  pjc  (p); 

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  PumpedLossyQbitSch   qbit(pplqb);
  PumpedLossyModeSch<> mode(pplm);

  JaynesCummings<> jc(qbit,mode,pjc);

  BinarySystem<> system(jc);

  double dtinit=.1/static_cast<structure::QuantumSystem<2>*>(&system)->highestFrequency();

  Array alpha(2);
  
  alpha(0)=pplm.minit;
  alpha(1)=pplqb.qbitInit;

  Simulated<Array> S(alpha,bind(derivs,_1,_2,_3,pplm,pplqb,pjc.g),dtinit,Array(),pt);

  evolve(S,pt);

  } catch (const ParsNamedException& pne) {cerr<<"Pars named error: "<<pne.getName()<<endl;}


}
