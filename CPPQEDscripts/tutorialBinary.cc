// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "EvolutionBinary.h"

#include "JaynesCummings.h"

int main(int argc, char* argv[])
{
  ParameterTable p;

  evolution::Pars<> pe(p); 

  qbit::ParsPumpedLossy pq(p);
  mode::ParsPumpedLossy pm(p);
  jaynescummings::Pars pjc(p);

  auto qmp=updateWithPicture(p,argc,argv,"--"); // Parsing the command line

  const auto psi{std::make_shared<quantumdata::StateVector<2>>(init(pq)*init(pm))};
  
  evolve(psi,
         binary::make(jaynescummings::make(make(pq,qmp),
                                           make(pm,qmp),
                                           pjc)),
         pe);

}
