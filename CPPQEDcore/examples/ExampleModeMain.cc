// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "EvolutionBinary.h"

#include "ExampleInteraction.h"


using std::make_shared;

int main(int argc, char* argv[])
{
  { // basic
  {
    basic::DrivenDissipativeModeIP mode(1,1,1i,1,10);
  }

  basic::InteractionX_X(make_shared<basic::DrivenDissipativeMode>(1,1,1i,1,10),make_shared<basic::DrivenDissipativeMode>(1,1,1i,1,10),2.);
  }
  
  { // hierarchical
  {
    hierarchical::DrivenDissipativeMode mode(1,1,1i,1,10);
  }

  auto m0{make_shared<hierarchical::DrivenDissipativeModeIP>(1,1,1i,1,10)};
  auto m1{make_shared<hierarchical::DrivenDissipativeMode  >(1,1,1i,1,20)};
  
  hierarchical::InteractionX_X i(m0,m1,2.);
  
  ParameterTable p;

  evolution::Pars<> pe(p); // Driver Parameters

  update(p,argc,argv,"--");

  typedef quantumdata::StateVector<2> StateVector;
  StateVector psi(StateVector::Dimensions(10,20)); psi(0,0)=1;
  
  evolve<0>(psi,binary::make(make_shared<hierarchical::InteractionX_X_Correlations>(m0,m1,2.)),pe);

  }
  
}
