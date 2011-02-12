#include "Algorithm.h"
#include "SimulatedHigh.h"

#include "ParsMode.h"
#include "Mode.h"

using namespace std       ;
using namespace cpputils  ;
using namespace trajectory;
using namespace mode      ;

typedef TTD_CARRAY(1) Array;

void derivs(double, const Array& b, Array& dbdt, const ParsPumpedLossy& p)
{
  dbdt(0)=dcomp(-p.kappa,p.delta)*b(0)+p.eta;
}


int main(int argc, char* argv[])
{
  ParameterTable p;

  ParsTrajectory pt(p);
  ParsPumpedLossy pplm(p); 

  // Parameter finalization
  update(p,argc,argv,"--");

  PumpedLossyMode<> plm(pplm);

  double dtinit=.1/plm.highestFrequency();
 
  Array alpha(1);
  
  alpha(0)=dcomp(1,-.5);

  Simulated<Array> S(alpha,bind(derivs,_1,_2,_3,pplm),dtinit,Array(),pt);
  
  evolve(S,pt);
  
}
