// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution.h"

#include "DrivenTwoLevelAtom.h"
#include "Qbit.h"

using namespace std;
using namespace qbit;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  evolution::Pars<> pe(p); // Driver Parameters
  ParsDrivenDissipativePhaseNoise pp2la(p);
  
  auto& fullImpl=p.addTitle("Script specific").add("fullImpl","Implement with Qbit (instead of DrivenTwoLevelAtomSch)",false);

  // Parameter finalization
  QM_Picture& qmp=updateWithPicture(p,argc,argv);
  
  // ****** ****** ****** ****** ****** ******

  evolve(init(pp2la),
         fullImpl ? std::static_pointer_cast<const structure::Free>(make(pp2la,qmp)) : 
                    std::static_pointer_cast<const structure::Free>(std::make_shared<const DrivenTwoLevelAtomSch>(pp2la)),
         pe);

}

