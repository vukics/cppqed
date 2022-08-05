// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "EvolutionBinary.h"

#include "JaynesCummings.h"


using namespace std;

typedef quantumdata::StateVector<2> StateVector;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem
  
  ParameterTable p;

  evolution::Pars<> pe(p); // Driver Parameters
  qbit::ParsDrivenDissipativePhaseNoise pplqb(p); 
  mode::ParsDrivenDissipative pplm (p); 
  jaynescummings::Pars  pjc  (p); 

  // Parameter finalization
  QM_Picture& qmp=updateWithPicture(p,argc,argv);
  
  // ****** ****** ****** ****** ****** ******

  /*
  mode::StateVector psiMode(pplm.cutoff);
  psiMode=(mode::fock(0,pplm.cutoff,mathutils::PI/2.)+mode::fock(1,pplm.cutoff)+mode::fock(2,pplm.cutoff,-mathutils::PI/2.))/sqrt(3.);
  */
  // entangled state: 
  // psi=qbit::state0()*mode::fock(1,pplm.cutoff)+qbit::state1()*mode::fock(0,pplm.cutoff);
  // psi.renorm();

  evolve(qbit::init(pplqb)*mode::init(pplm),
         binary::make(jaynescummings::make(qbit::make(pplqb,qmp),
                                           mode::make(pplm ,qmp/*,structure::averaged::ReducedDensityOperator("Mode",pplm.cutoff,true)*/),
                                           pjc)),
         pe);


}
