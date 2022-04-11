// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution_.h"
#include "Simulated.h"

#include "BichromaticMode.h"

using namespace trajectory;
using namespace mode      ;

typedef CArray<1> Array;


int main(int argc, char* argv[])
{
  ParameterTable p;

  evolution::Pars<> pt(p);
  ParsBichromatic pplm(p); 

  // Parameter finalization
  update(p,argc,argv,"--");

  BichromaticMode<> plm(pplm);

  double dtinit=.1/plm.highestFrequency();
 
  Array alpha(1);
  
  alpha=pplm.minit;

  run(simulated::makeBoost(alpha,[&](const Array& b, Array& dbdt, double t) {
    dbdt(0)=dcomp(-pplm.kappa,pplm.delta)*b(0)+pplm.eta+pplm.etaOther*exp(DCOMP_I*t*(pplm.deltaOther+pplm.delta));
  },{"alpha"},dtinit,pt),pt);
  
}
