#include "ExampleInteraction.h"

#include "SmartPtr.h"
using cpputils::sharedPointerize;


basic::InteractionX_X::InteractionX_X(const PumpedLossyMode& m0, const PumpedLossyMode& m1, double g)
  : Interaction<2>(Frees(sharedPointerize(m0),sharedPointerize(m1)),FREQS("g",g,sqrt(m0.getDimension()*m1.getDimension()))),
    TridiagonalHamiltonian<2,false>(g*
                                    (aop(m0.getDimension())+aop(m0.getDimension()).dagger())*
                                    (aop(m1.getDimension())+aop(m1.getDimension()).dagger())
                                    )
{}


hierarchical::InteractionX_X::InteractionX_X(const ModeBase& m0, const ModeBase& m1, double g)
  : Interaction<2>(Frees(sharedPointerize(m0),sharedPointerize(m1)),FREQS("g",g,sqrt(m0.getDimension()*m1.getDimension()))),
    TridiagonalHamiltonian<2,true>(g*
                                   (aop(m0)+aop(m0).dagger())*
                                   (aop(m1)+aop(m1).dagger())
                                   )
{}
