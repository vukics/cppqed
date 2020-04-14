// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "EvolutionBinary.h"
#include "Mode.h"

#include "ModeCorrelations.h"
#include "CoupledModes.h"


using namespace std;
using namespace mode;

typedef quantumdata::StateVector<2> StateVector2;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  evolution::Pars<> pe(p); // Driver Parameters
  ParsPumpedLossy pA(p,"Opt");
  ParsLossy       pB(p,"Mech"); 

  double& u=p.addTitle("Optomechanics","").add("u","N-Q coupling",1.);

  // Parameter finalization
  QM_Picture& qmp=updateWithPicture(p,argc,argv);
  
  // ****** ****** ****** ****** ****** ******

  const auto psi{std::make_shared<StateVector2>(init(pA)*init(pB))};
  psi->renorm();

  evolve
    (psi,
     binary::make(coupledmodes::make<coupledmodes::CM_NX,ModeCorrelations>(make<DoNotAverage>(pA,qmp),
                                                                           make<DoNotAverage>(pB,qmp),
                                                                           -sqrt(2)*u)),
     pe);


  // Another way which works:
  /*
  evolve<0>
    (psi,
     binary::make(coupledmodes::make<ModeCorrelations>(Mode<DoNotAverage>(pA),
                                                         Mode<DoNotAverage>(pB),
                                                         -sqrt(2)*u)),
     pe);
  */


  // A misuse:
  /*
  NX_CoupledModes<ModeCorrelations> nx(Mode<DoNotAverage>(pA),Mode<DoNotAverage>(pB),-sqrt(2)*u);
  // Here, it doesn't matter whether the interaction is created with or without its maker.

  evolve<0>
    (psi,
     binary::make(nx),
     pe);
  */

  // But this is fine again:
  /*
  Mode<DoNotAverage> m1(pA), m2(pB);

  NX_CoupledModes<ModeCorrelations> nx(m1,m2,-sqrt(2)*u);

  evolve<0>
    (psi,
     binary::make(nx),
     pe);
  */




}
