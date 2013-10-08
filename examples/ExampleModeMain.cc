#include "ExampleInteraction.h"

int main()
{
  {
  {
    basic::PumpedLossyModeIP mode(1,1,DCOMP_I,1,10);
  }

  basic::PumpedLossyMode m0(1,1,DCOMP_I,1,10), 
                         m1(1,1,DCOMP_I,1,10);
  
  basic::InteractionX_X(m0,m1,2.);
  }
  
  {
  {
    hierarchical::PumpedLossyMode mode(1,1,DCOMP_I,1,10);
  }

  hierarchical::PumpedLossyModeIP m0(1,1,DCOMP_I,1,10);
  hierarchical::PumpedLossyMode   m1(1,1,DCOMP_I,1,10);
  
  hierarchical::InteractionX_X(m0,m1,2.);
  }
  
}
