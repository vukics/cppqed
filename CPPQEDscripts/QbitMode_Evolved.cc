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

  
  double dtinit=.1/binary::make(jaynescummings::make(std::make_shared<const PumpedLossyQbitSch>(pplqb),std::make_shared<const PumpedLossyModeSch<>>(pplm),pjc))->highestFrequency();

  Array alpha(2);
  
  alpha=pplm.minit,pplqb.qbitInit;

  Simulated<Array> S(alpha,[&](double, const Array& b, Array& dbdt) {
    dbdt=
      dcomp(-pplm.kappa,pplm.delta)*b(0)+pjc.g*b(1)+pplm.eta,
      dcomp(-pplqb.gamma,pplqb.delta)*b(1)-pjc.g*b(0)-pplqb.eta;
  },dtinit,pt);

  run(S,pt);




}
