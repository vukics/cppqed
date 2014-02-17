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

