// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "GeneralDicke.h"

#include "Tridiagonal.tcc"


using namespace structure;
using namespace spin;
using namespace mode;



generaldicke::Base::Base(mode::Ptr mode, spin::Ptr spin, dcomp u, dcomp y)
  : Interaction<2>({mode,spin},{},{CF{"u",u,mode->getDimension()*spin->getDimension()},CF{"y",y,sqrt(mode->getDimension()*spin->getDimension())}}),
    TridiagonalHamiltonian<2,true>((u*nop(mode)*(sz(spin)+spin->getTwoS()/2.*quantumoperator::identity(spin->getDimension()))
                                   +
                                    y*(aop(mode).dagger()+aop(mode))*sx(spin))/DCOMP_I)
{
  getParsStream()<<"General Dicke interaction\n";
  // getParsStream()<<sx(spin->getDimension()-1)<<std::endl;
}
