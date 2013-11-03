#include "ParticleCavity_.h"

#include "MathExtensions.h"

#include "TridiagonalHamiltonian.tcc"


using namespace std;
using namespace mathutils;

using particle:: mfNKX;
using particle::cosNKX;

using namespace mode;


namespace particlecavity {


const Tridiagonal dispersive(mode::Ptr mode, particle::Ptr particle, double uNot, const ModeFunction& mf)
{
  return uNot*nop(mode)*(mf.get<0>()==MFT_SIN ? -1 : 1)*cosNKX(particle,mf.get<1>()<<1)/(2.*DCOMP_I);
}


const Tridiagonal interferic(mode::Ptr mode, particle::Ptr particle, double uNotTimesEtaeff, double uNot, const ModeFunction& mf)
{
  if (uNotTimesEtaeff>=0) return sign(uNot)*sqrt(uNotTimesEtaeff)*tridiagPlusHC_overI(aop(mode).dagger()*mfNKX(particle,mf));
  else                    throw UnotEtaeffSignDiscrepancy();
}


const Tridiagonals fillT(mode::Ptr mode, particle::Ptr particle, double uNot, double etaeff, const ModeFunction& mf)
{
  Tridiagonals res;

  if (uNot && !isComplex(mf.get<0>())) res.push_back(dispersive(mode,particle,uNot,mf));

  if (double factor=uNot*etaeff) res.push_back(interferic(mode,particle,factor,uNot,mf));

  return res;

}


Base::Base(mode::Ptr mode, particle::Ptr particle, double uNot, double etaeff)
  : structure::Interaction<2>(Frees(mode,particle),
			      FREQS("Unot",uNot,mode->getDimension())("etaeff",etaeff,sqrt(mode->getDimension()))
			      )
{
  getParsStream()<<"# Particle-Cavity Interaction\n";
}


InterferenceBase::InterferenceBase(mode::Ptr mode, particle::Ptr particle, double u, size_t kCav, ModeFunctionType modeCav)
  : MF_Base(modeCav,kCav),
    structure::Interaction<2>(Frees(mode,particle),
			      FREQS("u",u,mode->getDimension())),
    TridiagonalHamiltonian(interferic(mode,particle,sqr(u),u,MF_Base::member))
{
  getParsStream()<<"# Interference term with "<<getMF()<<endl;
}



POC_Base::POC_Base(mode::Ptr mode, particle::PtrPumped particle, double uNot)
  : particlecavity::Base(mode,particle,uNot,particle->getV_Class()),
    TridiagonalHamiltonian(interferic(mode,particle,uNot*particle->getV_Class(),uNot,particle->getMF()))
{
  getParsStream()<<"# Particle moving orthogonal to cavity\n"; /* photons/(particle number)^2="
								  <<uNot*particle->getV_Class()/sqrAbs(mode->getComplexFreqs(""))<<std::endl; */
}


PAC_Base::PAC_Base(mode::Ptr mode, particle::Ptr particle, double uNot, size_t kCav, ModeFunctionType modeCav, double etaeff)
  : MF_Base(modeCav,kCav),
    particlecavity::Base(mode,particle,uNot,etaeff),
    TridiagonalHamiltonian(fillT(mode,particle,uNot,etaeff,MF_Base::member))
{
  getParsStream()<<"# Particle moving along cavity with "<<getMF()<<endl;
}


} // particlecavity

