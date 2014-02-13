#include "ExampleInteraction.h"

#include "BinarySystem.h"
#include "BlitzArrayTraits.h"
#include "EnsembleMCWF.tcc"
#include "Evolution.tcc"
#include "Master.tcc"

#include "EvolvedGSL.tcc"
#include "Pars.tcc"


using namespace parameters;

int main(int argc, char* argv[])
{
  { // basic
  {
    basic::PumpedLossyModeIP mode(1,1,DCOMP_I,1,10);
  }

  basic::PumpedLossyMode m0(1,1,DCOMP_I,1,10), 
                         m1(1,1,DCOMP_I,1,10);
  
  basic::InteractionX_X(m0,m1,2.);
  }
  
  { // hierarchical
  {
    hierarchical::PumpedLossyMode mode(1,1,DCOMP_I,1,10);
  }

  hierarchical::PumpedLossyModeIP m0(1,1,DCOMP_I,1,10);
  hierarchical::PumpedLossyMode   m1(1,1,DCOMP_I,1,20);
  
  hierarchical::InteractionX_X i(m0,m1,2.);
  hierarchical::InteractionX_X_Correlations ii(m0,m1,2.);
  
  ParameterTable p;

  evolution::Pars pe(p); // Driver Parameters

  update(p,argc,argv,"--");

  typedef quantumdata::StateVector<2> StateVector;
  StateVector psi(StateVector::Dimensions(10,20)); psi(0,0)=1;
  
  evolve<tmptools::Vector<0> >(psi,binary::make(ii),pe);

  }
  
}
