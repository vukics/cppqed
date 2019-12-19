// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution_.h"
#include "Simulated.h"
#include "Mode.h"

#include "Qbit.h"

#include "JaynesCummings.h"
#include "BinarySystem.h"


using namespace std       ;
using namespace cpputils  ;
using namespace trajectory;

typedef CArray<1> Array;


void derivs(double, const Array& b, Array& dbdt, 
            const mode::ParsPumpedLossy& plm, const qbit::ParsPumpedLossy& pqb, const dcomp& g)
{
  dbdt=
    dcomp(-plm.kappa,plm.delta)*b(0)+g*b(1)+plm.eta,
    dcomp(-pqb.gamma,pqb.delta)*b(1)-g*b(0)-pqb.eta;
}


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  
  ParameterTable p;

  evolution::Pars<> pt(p);

  qbit::ParsPumpedLossy pplqb(p); 
  mode::ParsPumpedLossy pplm (p); 
  jaynescummings::Pars  pjc  (p); 

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  PumpedLossyQbitSch   qbit(pplqb);
  PumpedLossyModeSch<> mode(pplm);

  double dtinit=.1/binary::make(jaynescummings::make(qbit,mode,pjc))->highestFrequency();

  Array alpha(2);
  
  alpha=pplm.minit,pplqb.qbitInit;

  Simulated<Array> S(alpha,bind(derivs,_1,_2,_3,pplm,pplqb,pjc.g),dtinit,pt);

  run(S,pt);




}
