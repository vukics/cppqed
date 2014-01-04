#include "NX_CoupledModes.h"

#include "Tridiagonal.tcc"


using namespace mode;


nxcoupledmodes::Base<false>::Base(mode::Ptr m1, mode::Ptr m2, double u)
  : structure::Interaction<2>(Frees(m1,m2),RF{"u",u,m1->getDimension()*sqrt(m2->getDimension())})
{
  getParsStream()<<"# N-X coupling between modes\n";
}


nxcoupledmodes::Base<true >::Base(mode::Ptr m1, mode::Ptr m2, double u)
  : Base<false>(m1,m2,u), quantumoperator::TridiagonalHamiltonian<2,true>(u*nop(m1)*xop(m2)/DCOMP_I)
{
}
