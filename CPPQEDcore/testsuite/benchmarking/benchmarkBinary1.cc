// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution.h"

#include "benchmarking.h"

#include "Mode_.h"

using namespace mode;

const unsigned nRepeat=1000;

int main(int, char**)
{
  // ****** Parameters of the Problem

  ParameterTable p;

  ParsPumpedLossy pplm(p);

  pplm.eta=dcomp(4.,-2.);
  pplm.cutoff=1000;

  // ****** ****** ****** ****** ****** ******


  const PumpedLossyModeSch<false> modeH(pplm);
  const ModeBase mode0(pplm.cutoff);

  const Composite<composite::result_of::make_vector<Act<0,1> >::type> sys(makeComposite(Act<0,1>(structure::Interaction<2>(structure::Interaction<2>::Frees(&modeH,&mode0)))));
  
  benchmark(sys,modeH,Vector<0>());

}


