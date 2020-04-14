// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "EvolutionComposite.h"
#include "Mode.h"

using Interaction=structure::Interaction<4>;
using Frees=Interaction::Frees;


int main(int argc, char* argv[])
{
  ParameterTable p;

  evolution::Pars<> pe(p);

  mode::ParsPumpedLossy
    pm0(p,"0"),
    pm1(p,"1"),
    pm2(p,"2"),
    pm3(p,"3");

  update(p,argc,argv,"--");

  auto psi{std::make_shared<quantumdata::StateVector<4>>(init(pm0)*init(pm1)*init(pm2)*init(pm3))};

  evolve(psi,
         composite::make(
                         _<0,1,2,3>(std::make_shared<Interaction>(Frees(std::make_shared<PumpedLossyMode<>>(pm0),
                                                                        std::make_shared<PumpedLossyModeAlternative<false>>(pm1),
                                                                        std::make_shared<PumpedLossyMode<>>(pm2),
                                                                        std::make_shared<PumpedLossyModeAlternative<false>>(pm3))))
                        ),
         pe);


}
