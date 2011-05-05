#include "ExampleInteraction.h"

#include <boost/assign/list_of.hpp>
using boost::assign::tuple_list_of;

InteractionX_X::InteractionX_X(const PumpedLossyModeIP& m0, const PumpedLossyModeIP& m1, double g)
  : Interaction<2>(Frees(&m0,&m1),tuple_list_of("g",g,sqrt(m0.getDimension()*m1.getDimension()))),
    TridiagonalHamiltonian<2,true>(g*
				   (aop(m0.getDimension())+aop(m0.getDimension()).dagger())*
				   (aop(m1.getDimension())+aop(m1.getDimension()).dagger()),
				   freqs(m0.getDelta(),m0.getDimension())*freqs(m1.getDelta(),m1.getDimension())
				   )
{}
