#include "NX_CoupledModes.h"

#include "impl/Tridiagonal.tcc"

#include<boost/assign/list_of.hpp>

using namespace boost::assign;

using namespace mode;

NX_CoupledModes::NX_CoupledModes(mode::SmartPtr m1, mode::SmartPtr m2, double u)
  : structure::Interaction<2>(Frees(m1,m2),tuple_list_of("u",u,m1->getDimension()*sqrt(m2->getDimension()))),
    structure::TridiagonalHamiltonian<2,true>(u*nop(m1)*xop(m2)/DCOMP_I)
{
}
