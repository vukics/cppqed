// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Particle_.h"

#include "ParsParticle.h"
#include "ParticleInitialCondition.h"

#include "LazyDensityOperatorFFT.h"
#include "TridiagonalHamiltonian.h"

#include <boost/math/special_functions/hermite.hpp>

#include "TMP_Tools.h"

#include <boost/assign/list_of.hpp>
#include <utility>

using namespace std;
using namespace cppqedutils;
using namespace fft;



////////
//
// Exact
//
////////


namespace {

const particle::Exact::Diagonal fExpFill(const particle::Spatial& space, double omrec)
{
  return particle::Exact::Diagonal(-1i*omrec*blitz::sqr(space.getK()));
}

} 


particle::Exact::Exact(const Spatial& space, double omrec)
  : FreeExact(space.getDimension()), details::Storage(space,omrec), factorExponents_(fExpFill(space,omrec))
{
}

void particle::Exact::updateU(structure::OneTime t) const
{
  getDiagonal()=exp(factorExponents_*t);
}
 

//////////////
//
// Hamiltonian
//
//////////////


namespace {

const particle::Tridiagonal expINKX(size_t dim, ptrdiff_t nK)
{
  using particle::Tridiagonal;

  typedef Tridiagonal::Diagonal Diag;

  ptrdiff_t D=dim-(nK>0 ? nK : -nK);
  if (D<0) return Tridiagonal();

  Diag temp(D); temp=1;
  if (nK>0) return Tridiagonal(Diag(),nK,temp);
  else if (nK==0) return Tridiagonal(temp);
  return Tridiagonal(Diag(),-nK,Diag(),temp);
}


const particle::Tridiagonal cosNKX(size_t dim, ptrdiff_t nK)
{
  return (expINKX(dim,nK)+expINKX(dim,-nK))/2.;
}


const particle::Tridiagonal hOverI(size_t dim, double vClass, const ModeFunction& mf)
{
  return (vClass && !isComplex(get<0>(mf))) ? vClass*(get<0>(mf)==MFT_SIN ? -1 : 1)*cosNKX(dim,get<1>(mf)<<1)/(2.*1i) : particle::Tridiagonal();
}


const particle::Tridiagonal::Diagonal mainDiagonal(const particle::Spatial& space, double omrec)
{
  return particle::Tridiagonal::Diagonal(1i*omrec*blitz::sqr(space.getK()));
}


} 


namespace particle { // open it for the partial specializations

template<>
Hamiltonian<true >::Hamiltonian(const Spatial& space, double omrec, double vClass, const ModeFunction& mf)
  : Base(furnishWithFreqs(hOverI(space.getDimension(),vClass,mf),mainDiagonal(space,omrec))), Exact(space,omrec)
{
}


template<>
Hamiltonian<false>::Hamiltonian(const Spatial& space, double omrec, double vClass, const ModeFunction& mf)
  : Base(
         Tridiagonal(mainDiagonal(space,-omrec))
         +
         hOverI(space.getDimension(),vClass,mf)
         )
{
}


template<> template<>
Hamiltonian<false>::Hamiltonian(const Spatial& space, double omrec)
  : Base(
         Tridiagonal(mainDiagonal(space,-omrec))
         )
{
}

} // particle

///////////
//
// Averaged
//
///////////


particle::Averaged::Averaged(const Spatial& space)
  : Base("Particle",{"<P>","VAR(P)","<X>","DEV(X)"}),
    space_(space)
{
}


auto particle::Averaged::average_v(NoTime, const LazyDensityOperator& matrix) const -> const Averages
{
  int dim=space_.getDimension();

  auto averages(initializedAverages());
  
  const auto matrixX(quantumdata::ffTransform<tmptools::Vector<0> >(matrix,DIR_KX));
  
  for (int i=0; i<dim; i++) {

    double diag=matrix(i);
    averages(0)+=    space_.k(i) *diag;
    averages(1)+=sqr(space_.k(i))*diag;
  }

  for (int i=0; i<dim; i++) {
    double diag=matrixX->operator()(i);
    averages(2)+=    space_.x(i) *diag; 
    averages(3)+=sqr(space_.x(i))*diag;      
  }

  return averages;
}



void particle::Averaged::process_v(Averages& averages) const
{
  averages(1)-=sqr(averages(0));
  averages(3)=sqrt(averages(3)-sqr(averages(2)));
}



////////////////
//
// Highest level
//
////////////////

ParticleBase::ParticleBase(size_t fin, 
                           const RealFreqs& realFreqs, const ComplexFreqs& complexFreqs)
  : Free(1<<fin,realFreqs,complexFreqs), Averaged(particle::Spatial(fin))
{
  getParsStream()<<"Particle\n";
  getSpace().header(getParsStream());
}


DrivenParticleBase::DrivenParticleBase(size_t fin, double vClass, const ModeFunction& mf,
                                       const RealFreqs& realFreqs, const ComplexFreqs& complexFreqs)
  : ParticleBase(fin,
                 boost::assign::list_of(*realFreqs.begin()).range(next(realFreqs.begin()),realFreqs.end())({"vClass",vClass,1.}),
                 complexFreqs),
    vClass_(vClass), mf_(mf)
{
  getParsStream()<<"Pump "<<mf<<endl;
}



Particle::Particle(const particle::Pars& p)
  : ParticleBase(p.fin,{RF{"omrec",p.omrec,1<<p.fin}}),
    Exact(getSpace(),p.omrec)
{
}


ParticleSch::ParticleSch(const particle::Pars& p)
  : ParticleBase(p.fin,{RF{"omrec",p.omrec,sqr(1<<p.fin)}}),
    Hamiltonian<false>(getSpace(),p.omrec)
{
  getParsStream()<<"Schroedinger picture.\n";
}



DrivenParticle::DrivenParticle(const particle::ParsDriven& p)
  : DrivenParticleBase(p.fin,p.vClass,ModeFunction(p.modePart,p.kPart),{RF{"omrec",p.omrec,1<<p.fin}}),
    Hamiltonian<true>(getSpace(),p.omrec,p.vClass,getMF())
{
}


DrivenParticleSch::DrivenParticleSch(const particle::ParsDriven& p)
  : DrivenParticleBase(p.fin,p.vClass,ModeFunction(p.modePart,p.kPart),{RF{"omrec" ,p.omrec,sqr(1<<p.fin)}}),
    Hamiltonian<false>(getSpace(),p.omrec,p.vClass,getMF())
{
  getParsStream()<<"Schroedinger picture.\n";
}


//////////
//
// Helpers
//
//////////


namespace {

using particle::Spatial;

const Spatial::Array fill(size_t fin, double d, double m)
{
  Spatial::Array res(1<<fin);
  return res=d*blitz::tensor::i-m;
}

} 


particle::Spatial::Spatial(size_t fin, double deltaK)
  : fin_(fin),
    xMax_(PI/deltaK), deltaX_(2*xMax_/(1<<fin)), kMax_(PI/deltaX_), deltaK_(deltaK),
    x_(fill(fin,deltaX_,xMax_)), k_(fill(fin,deltaK_,kMax_))
{
}


void particle::Spatial::header(std::ostream& os) const
{
  os<<"Spatial Degree of Freedom finesse="<<fin_<<" xMax="<<xMax_<<" deltaX="<<deltaX_<<" kMax="<<kMax_<<" deltaK="<<deltaK_<<std::endl;
}



auto particle::expINKX(particle::Ptr particle, ptrdiff_t nK) -> const Tridiagonal
{
  size_t dim=particle->getDimension();
  Tridiagonal res(::expINKX(dim,nK));
  if (const auto exact=dynamic_cast<const particle::Exact*>(particle.get())) res.furnishWithFreqs(mainDiagonal(get<0>(*exact),get<1>(*exact)));
  return res;
}


auto particle::mfNKX(particle::Ptr particle, const ModeFunction& modeFunction) -> const Tridiagonal
{
  ModeFunctionType mf(get<0>(modeFunction));
  ptrdiff_t        nK(get<1>(modeFunction));
  switch (mf) {
  case MFT_SIN  : return sinNKX (particle,nK);
  case MFT_COS  : return cosNKX (particle,nK);
  case MFT_PLUS : return expINKX(particle,nK);
  case MFT_MINUS:                            ;
  }
  return expINKX(particle,-nK);
}

auto particle::mfComposition(particle::Ptr particle, const ModeFunction& modeFunction1, const ModeFunction& modeFunction2) -> const Tridiagonal
{
  ModeFunctionType mf1(get<0>(modeFunction1));
  ModeFunctionType mf2(get<0>(modeFunction2));
  ptrdiff_t        nK1(get<1>(modeFunction1));
  ptrdiff_t        nK2(get<1>(modeFunction2));
  ptrdiff_t        nK(abs(nK1));
  ptrdiff_t        sign1 = mf1==MFT_SIN && nK1<0 ? -1 : 1;
  ptrdiff_t        sign2 = mf2==MFT_SIN && nK2<0 ? -1 : 1;
  ptrdiff_t        sign  = sign1*sign2;

  if (abs(nK1) != abs(nK2)) throw std::runtime_error("Not a Tridiagonal");

  if (mf1==MFT_PLUS  && nK1<0) { mf1=MFT_MINUS; nK1*=-1; }
  if (mf1==MFT_MINUS && nK1<0) { mf1=MFT_PLUS;  nK1*=-1; }
  if (mf2==MFT_PLUS  && nK2<0) { mf2=MFT_MINUS; nK2*=-1; }
  if (mf2==MFT_MINUS && nK2<0) { mf2=MFT_PLUS;  nK2*=-1; }

  Tridiagonal id(quantumoperator::identity(particle->getSpace().getDimension()));

  typedef pair<ModeFunctionType,ModeFunctionType> MFPair;

  MFPair mfs(mf1,mf2);

  if (mfs == MFPair(MFT_SIN,MFT_SIN)){
    return sign * (id-cosNKX(particle,2*nK)) / 2.;
  }
  if (mfs == MFPair(MFT_COS,MFT_COS)){
    return (id+cosNKX(particle,2*nK)) / 2.;
  }
  if (mfs == MFPair(MFT_SIN,MFT_COS) || mfs == MFPair(MFT_COS,MFT_SIN)){
    return sign * sinNKX(particle,2*nK) / 2.;
  }
  if (mfs == MFPair(MFT_PLUS,MFT_PLUS) || mfs == MFPair(MFT_MINUS,MFT_MINUS)){
    return id;
  }
  if (mfs == MFPair(MFT_PLUS,MFT_MINUS)){
    return expINKX(particle,-2*nK);
  }
  if (mfs == MFPair(MFT_MINUS,MFT_PLUS)){
    return expINKX(particle,2*nK);
  }
  if (mfs == MFPair(MFT_PLUS,MFT_SIN) || mfs == MFPair(MFT_SIN,MFT_MINUS)) {
    return sign * (expINKX(particle, -2*nK)-id) * 1i/2.;
  }
  if (mfs == MFPair(MFT_PLUS,MFT_COS) || mfs == MFPair(MFT_COS,MFT_MINUS)) {
    return (id + expINKX(particle,-2*nK)) / 2.;
  }
  if (mfs == MFPair(MFT_MINUS,MFT_COS) || mfs == MFPair(MFT_COS,MFT_PLUS)) {
    return (id + expINKX(particle, 2*nK)) / 2.;
  }
  if (!(mfs == MFPair(MFT_MINUS,MFT_SIN) || mfs == MFPair(MFT_SIN,MFT_PLUS))) {
    // this should never be reached
    throw std::logic_error("In mfComposition");
  }
  return -sign * (expINKX(particle,2*nK)-id) * 1i/2.;
}

/*
const Tridiagonal mfNKX_AbsSqr(ModeFunctionType mf, size_t dim, ptrdiff_t K)
{
  return (mfcomplex(mf) ? Identity(dim) : (mf==MFT_Sin ? -1 : 1)*cosNKX(dim,2*K))/2.;
  // The other 1/2 should be taken into account separately!
}
*/

auto particle::wavePacket(const InitialCondition& init, const Spatial& space, bool kFlag) -> StateVector
{
  double 
    offset1=PI*init.getX0(),
    offset2=   init.getK0();

  if (init.isInK()) {swap(offset1,offset2); offset2=-offset2;}

  const Spatial::Array array(init.isInK() ? space.getK() : space.getX());

  StateVectorLow psiLow(exp(-blitz::sqr(array-offset1)/(4*sqr(init.getSig()))+1i*array*offset2));

  if      ( kFlag && !init.isInK()) quantumdata::ffTransform(psiLow,DIR_XK);
  else if (!kFlag &&  init.isInK()) quantumdata::ffTransform(psiLow,DIR_KX);

  StateVector res(psiLow,quantumdata::byReference); res.renorm();

  return res;

}


auto particle::wavePacket(const Pars& p, bool kFlag) -> StateVector
{
  return wavePacket(p.init,Spatial(p.fin),kFlag);
}


namespace {

const particle::InitialCondition coherent(const particle::ParsDriven& p)
{
  return particle::InitialCondition(p.init.getX0(),p.init.getK0(),pow(p.omrec/fabs(p.vClass),.25)/sqrt(2),false);
}

}

auto particle::wavePacket(const ParsDriven& p, bool kFlag) -> StateVector
{
  if (p.init.getSig()) return wavePacket(static_cast<const Pars&>(p),kFlag);
  else                 return wavePacket(coherent(p),
                                         Spatial(p.fin),
                                         kFlag);
}


auto particle::hoState(int n, const InitialCondition& init, const Spatial& space, bool kFlag) -> StateVector
{
  if (n<0) n=0;

  double kx0(PI*init.getX0());

  size_t dim=space.getDimension();

  StateVectorLow psiLow(dim);

  for (size_t j=0; j<dim; j++) {
    double temp=(space.getX()(j)-kx0)/init.getSig();
    psiLow(j)=exp(-sqr(temp)/2.)*boost::math::hermite(n,temp);
  }

  if (kFlag) quantumdata::ffTransform(psiLow,DIR_XK);
  
  StateVector res(psiLow,quantumdata::byReference); res.renorm();

  return res;

}


auto particle::hoState(const Pars      & p, bool kFlag) -> StateVector
{
  return hoState(p.hoInitn,p.init,Spatial(p.fin),kFlag);
}


auto particle::hoState(const ParsDriven& p, bool kFlag) -> StateVector
{
  if (p.init.getSig()) return hoState(static_cast<const Pars&>(p),kFlag);
  else                 return hoState(p.hoInitn,
                                      coherent(p),
                                      Spatial(p.fin),
                                      kFlag);
}


auto particle::init(const Pars& p) -> StateVector
{
  if (const auto pp=dynamic_cast<const ParsDriven*>(&p))
    return p.hoInitn<0 ? wavePacket(*pp) : hoState(*pp);
  else
    return p.hoInitn<0 ? wavePacket( p ) : hoState( p );
}




particle::Ptr particle::make(const Pars& p, QM_Picture qmp)
{
  return qmp==QMP_SCH ? Ptr(std::make_shared<ParticleSch>(p)) : Ptr(std::make_shared<Particle>(p));
}


particle::PtrDriven particle::makeDriven(const ParsDriven& p, QM_Picture qmp)
{
  return qmp==QMP_SCH ? PtrDriven(std::make_shared<DrivenParticleSch>(p)) : PtrDriven(std::make_shared<DrivenParticle>(p));
}
