// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution.h"
#include "Mode.h"

#include "BichromaticMode.h"
#include "AveragingUtils.h"


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
    &doStream=p.add("doStream","Stream diagonal elements of density operator",false),
    &doOffDiag=p.add("doOffDiag","Stream offdiagonal elements of density operator",false);

  // Parameter finalization
  QM_Picture& qmp=updateWithPicture(p,argc,argv);
  
  // ****** ****** ****** ****** ****** ******

  Collecting::Collection collection; collection.push_back(new AveragedQuadratures());
  if (doStream) collection.push_back(new ReducedDensityOperator<1>("",pplm.cutoff,doOffDiag));

  evolve(mode::init(pplm),(alternative ? Ptr(new PumpedLossyModeIP_NoExact(pplm)) : make<Collecting>(pplm,qmp,collection)),pe);

}

