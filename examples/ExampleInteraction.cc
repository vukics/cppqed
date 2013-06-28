#include "ExampleInteraction.h"

#include "SmartPtr.h"
using cpputils::sharedPointerize;


InteractionX_X::InteractionX_X(const PumpedLossyModeIP& m0, const PumpedLossyModeIP& m1, double g)
  : Interaction<2>(Frees(sharedPointerize(m0),sharedPointerize(m1)),FREQS("g",g,sqrt(m0.getDimension()*m1.getDimension()))),
    TridiagonalHamiltonian<2,true>(g*
                                   (aop(m0)+aop(m0).dagger())*
                                   (aop(m1)+aop(m1).dagger())
                                   )
{}
