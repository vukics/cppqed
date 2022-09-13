// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "EvolutionBinary.h"

#include "JaynesCummings.h"

int main(int argc, char* argv[])
{
  ParameterTable p;

  evolution::Pars<> pe(p); 

  qbit::ParsDrivenDissipative pq(p);
  mode::ParsDrivenDissipative pm(p);
  jaynescummings::Pars pjc(p);

  auto qmp=updateWithPicture(p,argc,argv,"--"); // Parsing the command line

  evolve(init(pq)*init(pm),
         binary::make(jaynescummings::make(make(pq,qmp),
                                           make(pm,qmp),
                                           pjc)),
         pe);

}
