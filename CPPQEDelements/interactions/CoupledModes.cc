// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "CoupledModes.h"

#include "Tridiagonal.tcc"


using namespace mode;
using namespace quantumoperator;


namespace coupledmodes {

Base<false>::Base(mode::Ptr m1, mode::Ptr m2, dcomp u)
  : structure::Interaction<2>({m1,m2},{},{CF{"u",u,m1->getDimension()*sqrt(m2->getDimension())}})
{
  if (!abs(u)) getParsStream()<<"Dummy coupling between modes\n";
}


template<>
Base<true ,CM_NX>::Base(mode::Ptr m1, mode::Ptr m2, dcomp u)
  : Base<false>(m1,m2,u), TridiagonalHamiltonian<2,true>(real(u)*nop(m1)*xop(m2)/1i)
{
  getParsStream()<<"N-X coupling between modes\n";
}


template<>
Base<true ,CM_XX>::Base(mode::Ptr m1, mode::Ptr m2, dcomp u)
  : Base<false>(m1,m2,u), TridiagonalHamiltonian<2,true>(real(u)*xop(m1)*xop(m2)/1i)
{
  getParsStream()<<"X-X coupling between modes\n";
}


template<>
Base<true ,CM_XX_RWA>::Base(mode::Ptr m1, mode::Ptr m2, dcomp u)
  : Base<false>(m1,m2,u), TridiagonalHamiltonian<2,true>(tridiagMinusHC(conj(u)*aop(m1)*aop(m2).dagger()))
{
  getParsStream()<<"X-X coupling in RWA between modes\n";
}

}
