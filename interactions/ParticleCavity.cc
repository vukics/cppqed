#include "ParticleCavity.h"

#include "MathExtensions.h"
#include "ParsParticleCavity.h"

#include "Mode.h"
#include "impl/TridiagonalHamiltonian.tcc"

#include<boost/assign/list_of.hpp>


using namespace std;
using namespace boost::assign;
using namespace mathutils;

using particle:: mfNKX;
using particle::cosNKX;

using namespace mode;


namespace particlecavity {


const Tridiagonal dispersive(const ModeBase* mode, const ParticleBase* particle, double uNot, const ModeFunction& mf)
{
  return uNot*nop(mode)*(mf.get<0>()==MFT_SIN ? -1 : 1)*cosNKX(particle,mf.get<1>()<<1)/(2.*DCOMP_I);
}


const Tridiagonal interferic(const ModeBase* mode, const ParticleBase* particle, double uNotTimesEtaeff, double uNot, const ModeFunction& mf)
{
  if (uNotTimesEtaeff>=0) return sign(uNot)*sqrt(uNotTimesEtaeff)*tridiagPlusHC_overI(aop(mode).dagger()*mfNKX(particle,mf));
  else                    throw UnotEtaeffSignDiscrepancy();
}


const Tridiagonals fillT(const ModeBase* mode, const ParticleBase* particle, double uNot, double etaeff, const ModeFunction& mf)
{
  Tridiagonals res;

  if (uNot && !isComplex(mf.get<0>())) res.push_back(dispersive(mode,particle,uNot,mf));

  if (double factor=uNot*etaeff) res.push_back(interferic(mode,particle,factor,uNot,mf));

  return res;

}


Base::Base(const ModeBase* mode, const ParticleBase* particle, double uNot, double etaeff)
  : structure::Interaction<2>(Frees(mode,particle),
			      tuple_list_of("Unot",uNot,mode->getDimension())("etaeff",etaeff,sqrt(mode->getDimension()))
			      )
{
  getParsStream()<<"# Particle-Cavity Interaction\n";
}


InterferenceBase::InterferenceBase(const ModeBase* mode, const ParticleBase* particle, double u, size_t kCav, ModeFunctionType modeCav)
  : MF_Base(modeCav,kCav),
    structure::Interaction<2>(Frees(mode,particle),
			      tuple_list_of("u",u,mode->getDimension())),
    TridiagonalHamiltonian(interferic(mode,particle,sqr(u),u,MF_Base::member))
{
  getParsStream()<<"# Interference term with "<<getMF()<<endl;
}



POC_Base::POC_Base(const ModeBase* mode, const PumpedParticleBase* particle, double uNot)
  : particlecavity::Base(mode,particle,uNot,particle->getV_Class()),
    TridiagonalHamiltonian(interferic(mode,particle,uNot*particle->getV_Class(),uNot,particle->getMF()))
{
  getParsStream()<<"# Particle moving orthogonal to cavity\n"; /* photons/(particle number)^2="
								  <<uNot*particle->getV_Class()/sqrAbs(mode->getComplexFreqs(""))<<std::endl; */
}


PAC_Base::PAC_Base(const ModeBase* mode, const ParticleBase* particle, double uNot, size_t kCav, ModeFunctionType modeCav, double etaeff)
  : MF_Base(modeCav,kCav),
    particlecavity::Base(mode,particle,uNot,etaeff),
    TridiagonalHamiltonian(fillT(mode,particle,uNot,etaeff,MF_Base::member))
{
  getParsStream()<<"# Particle moving along cavity with "<<getMF()<<endl;
}


} // particlecavity

