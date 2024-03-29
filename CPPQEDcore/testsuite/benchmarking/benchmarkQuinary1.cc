// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution.h"

#include "benchmarking.h"

#include "Mode_.h"

using namespace mode;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  ParsPumpedLossy pplm(p);

  pplm.eta=dcomp(4.,-2.);
  pplm.cutoff=16;

  // ****** ****** ****** ****** ****** ******

  const PumpedLossyModeSch<false> modeH(pplm);
  const ModeBase mode0(pplm.cutoff);

  const Composite<composite::result_of::make_vector<Act<0,1,2,3,4> >::type>
    sys(makeComposite(Act<0,1,2,3,4>(structure::Interaction<5>(structure::Interaction<5>::Frees(&mode0,&mode0,&modeH,&mode0,&mode0)))));
  
  benchmark(sys,modeH,Vector<2>());

}


