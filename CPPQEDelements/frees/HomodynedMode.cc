#include "HomodynedMode.h"


mode::HomodynedBase::HomodynedBase(const ParsLossy& p,
                                   const dcomp& homodyneAmplitude, const dcomp& eta)
  : Hamiltonian<false>(dcomp(finiteTemperatureHamiltonianDecay(p,boost::mpl::true_()),-p.delta),eta,p.cutoff),
    structure::ElementLiouvillean<1,2>("HomodynedMode",{"homodyned loss","homodyned absorption"}),
    homodyneAmplitude_(homodyneAmplitude), kappa_(p.kappa), nTh_(p.nTh)
{
  getH_OverIs().front()-=sqrt(p.kappa/2.)*(sqrt(p.nTh)+sqrt(p.nTh+1))*tridiagPlusHC(conj(homodyneAmplitude_)*aop(p.cutoff));
}


void mode::HomodynedBase::doActWithJ(NoTime, StateVectorLow& psi, LindbladNo<0>) const
{
  double fact=sqrt(2.*kappa_*(nTh_+1));
  int ubound=psi.ubound(0);
  for (int n=0; n<ubound; ++n)
    psi(n)=homodyneAmplitude_*psi(n)+fact*sqrt(n+1)*psi(n+1);
  psi(ubound)*=homodyneAmplitude_;
}

void mode::HomodynedBase::doActWithJ(NoTime, StateVectorLow& psi, LindbladNo<1>) const
{
  double fact=sqrt(2.*kappa_*nTh_);
  for (int n=psi.ubound(0); n>0; --n)
    psi(n)=conj(homodyneAmplitude_)*psi(n)+fact*sqrt(n)*psi(n-1);
  psi(0)*=conj(homodyneAmplitude_);
}
