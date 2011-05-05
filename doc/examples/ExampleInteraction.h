// -*- C++ -*-
#include "ExampleMode.h"
#include "Interaction.h"


class InteractionX_X
  : public Interaction<2>, public TridiagonalHamiltonian<2,true>
{
public:
  InteractionX_X(const PumpedLossyModeIP&, const PumpedLossyModeIP&, double g);

};
