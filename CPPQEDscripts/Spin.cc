// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution.h"
#include "Spin.h"

using namespace std ;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  evolution::Pars<> pe(p); // Driver Parameters

  spin::Pars ps(p);
  
  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  const auto spin{std::make_shared<LossySpinSch>(ps)};
  
  const auto psi{std::make_shared<structure::freesystem::StateVector>(spin->getDimensions())};

  (*psi)(psi->getArray().ubound(0))=1;
  
  evolve(psi,spin,pe);

}

