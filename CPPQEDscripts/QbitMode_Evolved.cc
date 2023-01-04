// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution_.h"
#include "Simulated.h"
#include "Mode.h"

#include "Qbit.h"

#include "JaynesCummings.h"
#include "BinarySystem.h"


using namespace std       ;
using namespace cppqedutils  ;
using namespace trajectory;

typedef CArray<1> Array;



int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  
  ParameterTable p;

  evolution::Pars<> pt(p);

  qbit::ParsDrivenDissipative pplqb(p); 
  mode::ParsDrivenDissipative pplm (p); 
  jaynescummings::Pars  pjc  (p); 

  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  
  double dtinit=.1/binary::make(jaynescummings::make(std::make_shared<const DrivenDissipativeQbitSch>(pplqb),std::make_shared<const DrivenDissipativeModeSch<>>(pplm),pjc))->highestFrequency();

  Array alpha(2);
  
  alpha=pplm.minit,pplqb.qbitInit;

  run(simulated::makeBoost(alpha,[&](const Array& b, Array& dbdt, double) {
    dbdt=
      dcomp(-pplm.kappa,pplm.delta)*b(0)+pjc.g*b(1)+pplm.eta,
      dcomp(-pplqb.gamma,pplqb.delta)*b(1)-pjc.g*b(0)-pplqb.eta;
  },{"alpha","sigma"},dtinit,pt),pt);




}
