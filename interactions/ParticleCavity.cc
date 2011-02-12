#include "ParticleCavity.h"

#include "ParsParticleCavity.h"

#include "Mode.h"

#include "ModeFunctionType.h"

#include<boost/assign/list_of.hpp>
#include<boost/assign/list_inserter.hpp>

#include<boost/tuple/tuple_io.hpp>


using namespace std;
using namespace boost::assign;
using namespace mathutils;

using particle:: mfNKX;
using particle::cosNKX;
using particle::ModeFunction;

using namespace mode;


namespace particlecavity {


const Tridiagonal dispersive(size_t modeDim, size_t particleDim, double uNot, const ModeFunction& mf)
{
  return uNot*nop(modeDim)*(mf.get<0>()==MFT_SIN ? -1 : 1)*cosNKX(particleDim,mf.get<1>()<<1)/(2.*DCOMP_I);
}


const Tridiagonal interferic(size_t modeDim, size_t particleDim, double uNotTimesEtaeff, double uNot, const ModeFunction& mf)
{
  if (uNotTimesEtaeff>=0) return sign(uNot)*sqrt(uNotTimesEtaeff)*tridiagPlusHC_overI(aop(modeDim).dagger()*mfNKX(particleDim,mf));
  else                    throw UnotEtaeffSignDiscrepancy();
}


const TridiagonalIPs fill(const ModeBase* mode, const ParticleBase* particle, double uNot, double etaeff, const ModeFunction& mf)
{
  TridiagonalIPs res;

  if (uNot && !isComplex(mf.get<0>())) push_back(res)(dispersive(mode->getDimension(),particle->getDimension(),uNot,mf),
						      freqs(mode)*freqs(particle,mf.get<1>()<<1));

  if (double factor=uNot*etaeff) push_back(res)(interferic(mode->getDimension(),particle->getDimension(),factor,uNot,mf),
						freqs(mode)*freqs(particle,mf.get<1>()));

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
    TridiagonalHamiltonian(interferic(mode->getDimension(),particle->getDimension(),sqr(u),u,MF_Base::member),
			   freqs(mode)*freqs(particle,MF_Base::member.get<1>()))
{
  getParsStream()<<"# Interference term with "<<getMF()<<endl;
}



POC_Base::POC_Base(const ModeBase* mode, const PumpedParticleBase* particle, double uNot)
  : particlecavity::Base(mode,particle,uNot,particle->getV_Class()),
    TridiagonalHamiltonian(interferic(mode->getDimension(),particle->getDimension(),uNot*particle->getV_Class(),uNot,particle->getMF()),
			   freqs(mode)*freqs(particle,particle->getMF().get<1>()))
{
  getParsStream()<<"# Particle moving orthogonal to cavity\n"; /* photons/(particle number)^2="
								  <<uNot*particle->getV_Class()/sqrAbs(mode->getComplexFreqs(""))<<std::endl; */
}


PAC_Base::PAC_Base(const ModeBase* mode, const ParticleBase* particle, double uNot, size_t kCav, ModeFunctionType modeCav, double etaeff)
  : MF_Base(modeCav,kCav),
    particlecavity::Base(mode,particle,uNot,etaeff),
    TridiagonalHamiltonian(fill(mode,particle,uNot,etaeff,MF_Base::member))
{
  getParsStream()<<"# Particle moving along cavity with "<<getMF()<<endl;
}


} // particlecavity

