// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution.h"
#include "Mode.h"

#include "HomodynedMode.h"
#include "AveragingUtils.tcc"


using namespace std ;
using namespace mode;


typedef averagingUtils::Collecting<1> Collecting;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  evolution::Pars<> pe(p); // Driver Parameters
  ParsHomodyned<ParsPumpedLossy> pplm(p);

  bool
    &doDisplay=p.add("doDisplay","Display diagonal elements of density operator",false),
    &doOffDiag=p.add("doOffDiag","Display offdiagonal elements of density operator",false);

  // Parameter finalization
  update(p,argc,argv);
  
  // ****** ****** ****** ****** ****** ******

  Collecting::Collection collection; collection.push_back(new AveragedQuadratures());
  if (doDisplay) collection.push_back(new ReducedDensityOperator<1>("",pplm.cutoff,doOffDiag));

  Ptr mode(make<Collecting>(pplm,collection));

  evolve(mode::init(pplm),mode,pe);

}

