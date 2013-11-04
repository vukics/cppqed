// -*- C++ -*-
#include "ExampleMode.h"

#include "Interaction.h"
#include "TridiagonalHamiltonian.tcc"


namespace basic {

class InteractionX_X
  : public Interaction<2>, public TridiagonalHamiltonian<2,false>
{
public:
  InteractionX_X(const PumpedLossyMode&, const PumpedLossyMode&, double g);

};

} // basic


namespace hierarchical {

class InteractionX_X
  : public Interaction<2>, public TridiagonalHamiltonian<2,true>
{
public:
  InteractionX_X(const ModeBase&, const ModeBase&, double g);

};

} // hierarchical