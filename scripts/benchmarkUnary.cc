#include "EvolutionHigh.h"

#include "benchmarking.h"

#include "Mode.h"

using namespace mode;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  ParsPumpedLossy pplm(p);

  pplm.eta=dcomp(4.,-2.);
  pplm.cutoff=1000000;

  // ****** ****** ****** ****** ****** ******

  const PumpedLossyModeSch<false> mode(pplm);

  benchmark(mode);

}


