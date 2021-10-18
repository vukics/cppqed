// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ParticleCavity_.h"

#include "TridiagonalHamiltonian.h"

#include "ContainerIO.h"
#include "MathExtensions.h"
#include "SliceIterator.tcc"

using namespace std;
using namespace cppqedutils;

using particle:: mfNKX;
using particle::cosNKX;
using particle::mfComposition;


namespace {

double etaEff(double uNot, double vClass)
{
  if (uNot*vClass>=0) return sign(uNot)*sqrt(uNot*vClass);
  else                throw std::domain_error("Unot VClass sign discrepancy");
}

const particlecavity::Tridiagonals dispersive(mode::Ptr mode, particle::Ptr particle, double uNot, const ModeFunction& mf)
{
  particlecavity::Tridiagonals res;
  if (uNot && !isComplex(get<0>(mf)))
    res.push_back(uNot*mode::nop(mode)*(get<0>(mf)==MFT_SIN ? -1 : 1)*cosNKX(particle,get<1>(mf)<<1)/(2.*DCOMP_I));
  return res;
}

particlecavity::Tridiagonal interfericOneModeFN(mode::Ptr mode, particle::Ptr particle, double etaeff, const ModeFunction& mf)
{
  return etaeff*tridiagPlusHC_overI(mode::aop(mode)*mfNKX(particle,mf));
}

const particlecavity::Tridiagonal interfericTwoModeFN(mode::Ptr mode, particle::Ptr particle, double etaeff, const ModeFunction& mfCav, const ModeFunction& mfPart)
{
  return etaeff*tridiagPlusHC_overI(mode::aop(mode).dagger()*mfComposition(particle, mfCav, mfPart));
}

const particlecavity::Tridiagonals fillT(mode::Ptr mode, particle::Ptr particle, double uNot, double etaeff, const ModeFunction& mf)
{
  particlecavity::Tridiagonals res = dispersive(mode,particle,uNot,mf);

  if (etaeff) res.push_back(interfericOneModeFN(mode,particle,etaeff,mf));

  return res;
}

const particlecavity::Tridiagonals fillTTwoModeFN(mode::Ptr mode, particle::PtrPumped particle, double uNot, double etaeff, const ModeFunction& mf, bool isSpecial)
{
  particlecavity::Tridiagonals res = dispersive(mode,particle,uNot,mf);
  const ModeFunction mfPart = particle->getMF();

  if (etaeff && isSpecial) res.push_back(interfericTwoModeFN(mode,particle,etaeff,mf,mfPart));

  return res;
}

} // anonymous namespace


particlecavity::Base::Base(mode::Ptr mode, particle::Ptr particle, double uNot, double etaeff)
  : structure::Interaction<2>({mode,particle},{RF{"Unot",uNot,mode->getDimension()},RF{"etaeff",etaeff,sqrt(mode->getDimension())}})
{
  getParsStream()<<"Particle-Cavity Interaction\n";
}

ParticleOrthogonalToCavity::ParticleOrthogonalToCavity(mode::Ptr mode, particle::PtrPumped particle, double uNot)
  : ParticleOrthogonalToCavity(mode, particle, uNot, etaEff(uNot, particle->getV_Class())) {}


ParticleOrthogonalToCavity::ParticleOrthogonalToCavity(mode::Ptr mode, particle::PtrPumped particle, double uNot, double etaeff)
  : particlecavity::Base(mode,particle,uNot,etaeff),
    TridiagonalHamiltonian(interfericOneModeFN(mode,particle,etaeff,particle->getMF()))
{
  getParsStream()<<"Particle moving orthogonal to cavity\n"; /* photons/(particle number)^2="
                                                                  <<uNot*particle->getV_Class()/sqrAbs(mode->getComplexFreqs(""))<<std::endl; */
}


ParticleAlongCavity::ParticleAlongCavity(mode::Ptr mode, particle::Ptr particle, const particlecavity::ParsAlong& p, double vClass, const ThePrivateOne&)
  : ParticleAlongCavity(mode, particle, p.uNot, p.kCav, p.modeCav, etaEff(p.uNot,vClass)) {}

ParticleAlongCavity::ParticleAlongCavity(mode::Ptr mode, particle::PtrPumped particle, const particlecavity::ParsAlong& p, const ThePrivateOne&)
  : ParticleAlongCavity(mode, particle, p.uNot, p.kCav, p.modeCav, etaEff(p.uNot,particle->getV_Class())) {}

ParticleAlongCavity::ParticleAlongCavity(mode::Ptr mode, particle::Ptr particle, double uNot, size_t kCav, ModeFunctionType modeCav, double etaeff)
  : MF_Base(modeCav,kCav),
    particlecavity::Base(mode,particle,uNot,etaeff),
    isSpecialH_(true),
    tridiagonalH_(fillT(mode,particle,uNot,etaeff,MF_Base::member)),
    firstH_(quantumoperator::identity(0)*quantumoperator::identity(0)),firstHT_(firstH_),    // not used
    secondH_(quantumoperator::identity(0)),secondHT_(secondH_)  // not used
{
  getParsStream()<<"Particle moving along cavity with "<<getMF()<<endl;
}

ParticleAlongCavity::ParticleAlongCavity(mode::Ptr mode, particle::PtrPumped particle, double uNot, size_t kCav, ModeFunctionType modeCav, double etaeff)
  : MF_Base(modeCav,kCav),
    particlecavity::Base(mode,particle,uNot,etaeff),
    isSpecialH_(kCav==abs(get<1>(particle->getMF()))),
    tridiagonalH_(fillTTwoModeFN(mode,particle,uNot,etaeff,MF_Base::member,isSpecialH_)),
    firstH_(etaeff*mode::aop(mode)*mfNKX(particle,MF_Base::member)/DCOMP_I), firstHT_(-firstH_.dagger()),
    secondH_(mfNKX(particle,particle->getMF()).dagger()), secondHT_(secondH_.dagger())
{
  getParsStream()<<"Particle with "<< particle->getMF() <<" pump moving along cavity with "<<getMF()<<endl;
}

void ParticleAlongCavity::addContribution_v(double t, const quantumdata::StateVectorLow<2>& psi, quantumdata::StateVectorLow<2>& dpsidt, double t0) const
{
  tridiagonalH_.addContribution(t,psi,dpsidt,t0);

  if (isSpecialH_) return;

  using cppqedutils::sliceiterator::fullRange;

  typedef tmptools::Vector<1> V1;

  {
    double dt=t-t0;
    firstH_ .propagate(dt); firstHT_. propagate(dt);
    secondH_.propagate(dt); secondHT_.propagate(dt);
  }

  quantumdata::StateVectorLow<2> dpsidtTemp{psi.shape()}; // NEEDS_WORK check whether putting this into class scope saves time (only one dynamic allocation)
  {
    dpsidtTemp=0;
    apply(psi,dpsidtTemp,firstH_);
    boost::for_each(fullRange<V1>(dpsidtTemp),fullRange<V1>(dpsidt),[&](const auto& p1, auto& p2){secondH_ .apply(p1,p2);});
  }
  {
    dpsidtTemp=0;
    apply(psi,dpsidtTemp,firstHT_);
    boost::for_each(fullRange<V1>(dpsidtTemp),fullRange<V1>(dpsidt),[&](const auto& p1, auto& p2){secondHT_.apply(p1,p2);});
  }
}
