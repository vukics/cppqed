#include "PumpedTwoLevelAtom.h"

#include "ParsQbit.h"

#include "LazyDensityOperator.h"

#include <boost/assign/list_of.hpp>


using namespace linalg;
using namespace std;

using namespace boost;
using namespace assign;


PumpedTwoLevelAtom::PumpedTwoLevelAtom(const qbit::ParsPumpedLossy& pp2la)
  : Free(2,RealFreqs(),tuple_list_of("(gamma,deltaA)",-conj(dcomp(-pp2la.gamma,pp2la.delta)),1.)("eta",pp2la.eta,1.)), 
    Base("Pumped Two-Level Atom","atomic decay"),
    Averaged("Pumped Two-Level Atom",list_of("rho00")("rho11")("real(rho01)")("imag(rho01)")),
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


const PumpedTwoLevelAtom::Averages PumpedTwoLevelAtom::average(const LazyDensityOperator& M) const
// Calculates the complete density matrix
{
  Averages avr(4);
  avr(0)=M(0);
  avr(1)=M(1);
  avr(2)=real(M(1,0));
  avr(3)=imag(M(1,0));
  return avr;
}


void PumpedTwoLevelAtom::process(Averages&) const
{
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
