// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution.h"
#include "Mode.h"

#include "BichromaticMode.h"
#include "AveragingUtils.tcc"


using namespace std ;
using namespace mode;


typedef averagingUtils::Collecting<1> Collecting;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  evolution::Pars<> pe(p); // Driver Parameters
  ParsBichromatic pplm(p); 

  bool
    &alternative=p.add("alternative","Alternative mode",false),
    &doDisplay=p.add("doDisplay","Display diagonal elements of density operator",false),
    &doOffDiag=p.add("doOffDiag","Display offdiagonal elements of density operator",false);

  // Parameter finalization
  QM_Picture& qmp=updateWithPicture(p,argc,argv);
  
  // ****** ****** ****** ****** ****** ******

  Collecting::Collection collection; collection.push_back(new AveragedQuadratures());
  if (doDisplay) collection.push_back(new ReducedDensityOperator<1>("",pplm.cutoff,doOffDiag));

  Ptr mode(alternative ? Ptr(new PumpedLossyModeIP_NoExact(pplm)) : make<Collecting>(pplm,qmp,collection));

  StateVector psi(mode::init(pplm));

  evolve(psi,mode,pe);




}

