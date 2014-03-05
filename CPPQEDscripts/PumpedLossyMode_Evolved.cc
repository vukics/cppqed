#include "Evolution_.h"
#include "Simulated.h"

#include "BichromaticMode.h"

using namespace trajectory;
using namespace mode      ;

typedef CArray<1> Array;

void derivs(double t, const Array& b, Array& dbdt, const ParsBichromatic& p)
{
  dbdt(0)=dcomp(-p.kappa,p.delta)*b(0)+p.eta+p.etaOther*exp(DCOMP_I*t*(p.deltaOther+p.delta));
}


int main(int argc, char* argv[])
{
  ParameterTable p;

  evolution::Pars pt(p);
  ParsBichromatic pplm(p); 

  // Parameter finalization
  update(p,argc,argv,"--");

  BichromaticMode<> plm(pplm);

  double dtinit=.1/plm.highestFrequency();
 
  Array alpha(1);
  
  alpha=pplm.minit;

  Simulated<Array> S(alpha,bind(derivs,_1,_2,_3,pplm),dtinit,pt);
  
  run(S,pt);
  
}
