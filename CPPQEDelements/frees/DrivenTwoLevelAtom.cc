// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "DrivenTwoLevelAtom.h"

#include "ParsQbit.h"

#include "LazyDensityOperator.h"

using namespace linalg;
using namespace std;
using namespace structure;


DrivenTwoLevelAtom::DrivenTwoLevelAtom(const qbit::ParsDrivenDissipative& pp2la)
  : Free(2,{},{CF{"(gamma,deltaA)",{pp2la.gamma,pp2la.delta},1.},CF{"eta",pp2la.eta,1.}}), 
    structure::ElementLiouvillian<1,1>("Driven Two-Level Atom","atomic decay"),
    za_(-pp2la.gamma,pp2la.delta), eta_(pp2la.eta)
{
  getParsStream()<<"Driven Two-Level Atom\n";
}


double DrivenTwoLevelAtom::rate(NoTime, const qbit::LazyDensityOperator& matrix) const
{
  return -2.*real(za_)*matrix(1);
}


void DrivenTwoLevelAtom::doActWithJ(NoTime, qbit::StateVectorLow& psi) const
{
  psi(0)=sqrt(-2.*real(za_))*psi(1);
  psi(1)=0;
}


DrivenTwoLevelAtomSch::DrivenTwoLevelAtomSch(const qbit::ParsDrivenDissipative& pp2la)
  : DrivenTwoLevelAtom(pp2la), hamiltonianOverI_(hamiltonianOverI(getZa(),pp2la.eta)) 
{
  getParsStream()<<"Schroedinger picture.\n"; 
}


const CMatrix DrivenTwoLevelAtomSch::hamiltonianOverI(dcomp za, dcomp eta)
{
  CMatrix res(2,2);
  res=
    0,    conj(eta),
    -eta, za       ;
  return res;
}
