// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution.h"
#include "Mode.h"

#include <chrono>

using namespace std ;
using namespace mode;


int main(int argc, char* argv[])
{
  ParameterTable p;

  evolution::Pars<> pe(p); // Driver Parameters
  ParsPumpedLossy pplm(p); 

  size_t &nRep=p.add("nRep","numberOfBenchmarkRepetitions",size_t(10));

  // Parameter finalization
  QM_Picture& qmp=updateWithPicture(p,argc,argv);
  
  // ****** ****** ****** ****** ****** ******

  for (size_t i=0; i<nRep; ++i) {
    auto psi{mode::init(pplm)};
    auto m{make(pplm,qmp)};
    evolve(psi,m,pe);
  }

}

