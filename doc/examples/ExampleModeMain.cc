#include "ExampleInteraction.h"

int main()
{
  {
    PumpedLossyMode mode(1,1,DCOMP_I,1,10);
  }

  PumpedLossyModeIP m0(1,1,DCOMP_I,1,10);
  PumpedLossyModeIP m1(1,1,DCOMP_I,1,10);
  
  InteractionX_X(m0,m1,2.);

}
