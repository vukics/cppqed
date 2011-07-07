#include "ExampleInteraction.h"

#include <boost/assign/list_of.hpp>
using boost::assign::tuple_list_of;

InteractionX_X::InteractionX_X(const PumpedLossyModeIP& m0, const PumpedLossyModeIP& m1, double g)
  : Interaction<2>(Frees(&m0,&m1),tuple_list_of("g",g,sqrt(m0.getDimension()*m1.getDimension()))),
    TridiagonalHamiltonian<2,true>(g*
				   (aop(m0)+aop(m0).dagger())*
				   (aop(m1)+aop(m1).dagger())
				   )
{}
