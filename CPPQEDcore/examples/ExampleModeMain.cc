// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ExampleInteraction.h"

#include "BinarySystem.h"
#include "BlitzArrayTraits.h"
#include "EnsembleMCWF.tcc"
#include "Evolution.tcc"
#include "Master.tcc"

#include "EvolvedGSL.tcc"
#include "Pars.h"


using namespace parameters;

int main(int argc, char* argv[])
{
  { // basic
  {
    basic::PumpedLossyModeIP mode(1,1,DCOMP_I,1,10);
  }

  basic::InteractionX_X(std::make_shared<basic::PumpedLossyMode>(1,1,DCOMP_I,1,10),std::make_shared<basic::PumpedLossyMode>(1,1,DCOMP_I,1,10),2.);
  }
  
  { // hierarchical
  {
    hierarchical::PumpedLossyMode mode(1,1,DCOMP_I,1,10);
  }

  auto m0{std::make_shared<hierarchical::PumpedLossyModeIP>(1,1,DCOMP_I,1,10)};
  auto m1{std::make_shared<hierarchical::PumpedLossyMode  >(1,1,DCOMP_I,1,20)};
  
  hierarchical::InteractionX_X i(m0,m1,2.);
  
  ParameterTable p;

  evolution::Pars<> pe(p); // Driver Parameters

  update(p,argc,argv,"--");

  typedef quantumdata::StateVector<2> StateVector;
  StateVector psi(StateVector::Dimensions(10,20)); psi(0,0)=1;
  
  evolve<tmptools::Vector<0> >(std::make_shared<StateVector>(psi),binary::make(std::make_shared<hierarchical::InteractionX_X_Correlations>(m0,m1,2.)),pe);

  }
  
}
