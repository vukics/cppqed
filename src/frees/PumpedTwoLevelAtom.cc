#include "PumpedTwoLevelAtom.h"

#include "ParsQbit.h"

#include "LazyDensityOperator.h"

using namespace linalg;
using namespace std;


PumpedTwoLevelAtom::PumpedTwoLevelAtom(const qbit::ParsPumpedLossy& pp2la)
  : Free(2,RealFreqs(),FREQS("(gamma,deltaA)",-conj(dcomp(-pp2la.gamma,pp2la.delta)),1.)("eta",pp2la.eta,1.)), 
    Base("Pumped Two-Level Atom","atomic decay"),
    za_(-pp2la.gamma,pp2la.delta), eta_(pp2la.eta)
{
  getParsStream()<<"# Pumped Two-Level Atom\n";
}


double PumpedTwoLevelAtom::probability(const LazyDensityOperator& matrix) const
{
  return -2.*real(za_)*matrix(1);
}


void PumpedTwoLevelAtom::doActWithJ(StateVectorLow& psi) const
{
  psi(0)=sqrt(-2.*real(za_))*psi(1);
  psi(1)=0;
}


PumpedTwoLevelAtomSch::PumpedTwoLevelAtomSch(const qbit::ParsPumpedLossy& pp2la)
  : Base(pp2la), hamiltonianOverI_(hamiltonianOverI(getZa(),pp2la.eta)) 
{
  getParsStream()<<"# Schroedinger picture.\n"; 
}


const CMatrix PumpedTwoLevelAtomSch::hamiltonianOverI(const dcomp& za, const dcomp& eta)
{
  CMatrix res(2,2);
  res=
    0,    conj(eta),
    -eta, za       ;
  return res;
}
