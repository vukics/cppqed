// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Evolution.h"
#include "Spin.h"

using namespace std ;


int main(int argc, char* argv[])
{
  // ****** Parameters of the Problem

  ParameterTable p;

  evolution::Pars pe(p); // Driver Parameters

  spin::Pars ps(p);
  
  // Parameter finalization
  update(p,argc,argv,"--");
  
  // ****** ****** ****** ****** ****** ******

  LossySpin spin(ps);
  
  structure::freesystem::StateVector psi(spin.getDimensions());

  psi(psi.getArray().ubound(0))=1;
  
  evolve(psi,spin,pe);

}

