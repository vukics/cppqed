// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "Particle_.h"

#include "ParsParticle.h"
#include "ParticleInitialCondition.h"

#include "LazyDensityOperatorFFT.tcc"
#include "TridiagonalHamiltonian.tcc"

#include <boost/math/special_functions/hermite.hpp>

#include "TMP_Tools.h"

#include <boost/assign/list_of.hpp>
#include <boost/bind.hpp>

using namespace std;
using namespace mathutils;
using namespace cpputils;
using namespace fft;

using boost::make_shared;


////////
//
// Exact
//
////////


namespace {

const particle::Exact::Diagonal fExpFill(const particle::Spatial& space, double omrec)
{
  return particle::Exact::Diagonal(-DCOMP_I*omrec*blitz::sqr(space.getK()));
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
  return (vClass && !isComplex(mf.get<0>())) ? vClass*(mf.get<0>()==MFT_SIN ? -1 : 1)*cosNKX(dim,mf.get<1>()<<1)/(2.*DCOMP_I) : particle::Tridiagonal();
}


const particle::Tridiagonal::Diagonal mainDiagonal(const particle::Spatial& space, double omrec)
{
  return particle::Tridiagonal::Diagonal(DCOMP_I*omrec*blitz::sqr(space.getK()));
}


} 


template<>
particle::Hamiltonian<true >::Hamiltonian(const Spatial& space, double omrec, double vClass, const ModeFunction& mf)
  : Base(furnishWithFreqs(hOverI(space.getDimension(),vClass,mf),mainDiagonal(space,omrec))), Exact(space,omrec)
{
}


template<>
particle::Hamiltonian<false>::Hamiltonian(const Spatial& space, double omrec, double vClass, const ModeFunction& mf)
  : Base(
         Tridiagonal(mainDiagonal(space,-omrec))
         +
         hOverI(space.getDimension(),vClass,mf)
         )
{
}


template<>
particle::Hamiltonian<false>::Hamiltonian(const Spatial& space, double omrec, boost::mpl::bool_<false>)
  : Base(
         Tridiagonal(mainDiagonal(space,-omrec))
         )
{
}


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
  
  const LazyDensityOperator::Ptr matrixX(quantumdata::ffTransform<tmptools::Vector<0> >(matrix,DIR_KX));
  
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
  getParsStream()<<"# Particle\n";
  getSpace().header(getParsStream());
}


PumpedParticleBase::PumpedParticleBase(size_t fin, double vClass, const ModeFunction& mf,
                                       const RealFreqs& realFreqs, const ComplexFreqs& complexFreqs)
  : ParticleBase(fin,
                 boost::assign::list_of(*realFreqs.begin()).range(next(realFreqs.begin()),realFreqs.end())(make_tuple("vClass",vClass,1.)),
                 complexFreqs),
    vClass_(vClass), mf_(mf)
{
  getParsStream()<<"# Pump "<<mf<<endl;
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
  getParsStream()<<"# Schroedinger picture.\n";
}



PumpedParticle::PumpedParticle(const particle::ParsPumped& p)
  : PumpedParticleBase(p.fin,p.vClass,ModeFunction(p.modePart,p.kPart),{RF{"omrec",p.omrec,1<<p.fin}}),
    Hamiltonian<true>(getSpace(),p.omrec,p.vClass,getMF())
{
}


PumpedParticleSch::PumpedParticleSch(const particle::ParsPumped& p)
  : PumpedParticleBase(p.fin,p.vClass,ModeFunction(p.modePart,p.kPart),{RF{"omrec" ,p.omrec,sqr(1<<p.fin)}}),
    Hamiltonian<false>(getSpace(),p.omrec,p.vClass,getMF())
{
  getParsStream()<<"# Schroedinger picture.\n";
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
  os<<"# Spatial Degree of Freedom finesse="<<fin_<<" xMax="<<xMax_<<" deltaX="<<deltaX_<<" kMax="<<kMax_<<" deltaK="<<deltaK_<<std::endl;
}



auto particle::expINKX(particle::Ptr particle, ptrdiff_t nK) -> const Tridiagonal
{
  size_t dim=particle->getDimension();
  Tridiagonal res(::expINKX(dim,nK));
  if (const auto exact=dynamic_cast<const particle::Exact*>(particle.get())) res.furnishWithFreqs(mainDiagonal(exact->get<0>(),exact->get<1>()));
  return res;
}


auto particle::mfNKX(particle::Ptr particle, const ModeFunction& modeFunction) -> const Tridiagonal
{
  ModeFunctionType mf(modeFunction.get<0>());
  ptrdiff_t        nK(modeFunction.get<1>());
  switch (mf) {
  case MFT_SIN  : return sinNKX (particle,nK);
  case MFT_COS  : return cosNKX (particle,nK);
  case MFT_PLUS : return expINKX(particle,nK);
  case MFT_MINUS:                            ;
  }
  return expINKX(particle,-nK);
}

/*
const Tridiagonal mfNKX_AbsSqr(ModeFunctionType mf, size_t dim, ptrdiff_t K)
{
  return (mfcomplex(mf) ? Identity(dim) : (mf==MFT_Sin ? -1 : 1)*cosNKX(dim,2*K))/2.;
  // The other 1/2 should be taken into account separately!
}
*/

auto particle::wavePacket(const InitialCondition& init, const Spatial& space, bool kFlag) -> const StateVector
{
  double 
    offset1=PI*init.getX0(),
    offset2=   init.getK0();

  if (init.isInK()) {swap(offset1,offset2); offset2=-offset2;}

  const Spatial::Array array(init.isInK() ? space.getK() : space.getX());

  StateVectorLow psiLow(exp(-blitz::sqr(array-offset1)/(4*sqr(init.getSig()))+DCOMP_I*array*offset2));

  if      ( kFlag && !init.isInK()) quantumdata::ffTransform(psiLow,DIR_XK);
  else if (!kFlag &&  init.isInK()) quantumdata::ffTransform(psiLow,DIR_KX);

  StateVector res(psiLow,quantumdata::byReference); res.renorm();

  return res;

}


auto particle::wavePacket(const Pars& p, bool kFlag) -> const StateVector
{
  return wavePacket(p.init,Spatial(p.fin),kFlag);
}


namespace {

const particle::InitialCondition coherent(const particle::ParsPumped& p)
{
  return particle::InitialCondition(p.init.getX0(),p.init.getK0(),pow(p.omrec/fabs(p.vClass),.25)/sqrt(2),false);
}

}

auto particle::wavePacket(const ParsPumped& p, bool kFlag) -> const StateVector
{
  if (p.init.getSig()) return wavePacket(static_cast<const Pars&>(p),kFlag);
  else                 return wavePacket(coherent(p),
                                         Spatial(p.fin),
                                         kFlag);
}


auto particle::hoState(int n, const InitialCondition& init, const Spatial& space, bool kFlag) -> const StateVector
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


auto particle::hoState(const Pars      & p, bool kFlag) -> const StateVector
{
  return hoState(p.hoInitn,p.init,Spatial(p.fin),kFlag);
}


auto particle::hoState(const ParsPumped& p, bool kFlag) -> const StateVector
{
  if (p.init.getSig()) return hoState(static_cast<const Pars&>(p),kFlag);
  else                 return hoState(p.hoInitn,
                                      coherent(p),
                                      Spatial(p.fin),
                                      kFlag);
}


auto particle::init(const Pars& p) -> const StateVector
{
  if (const auto pp=dynamic_cast<const ParsPumped*>(&p))
    return p.hoInitn<0 ? wavePacket(*pp) : hoState(*pp);
  else
    return p.hoInitn<0 ? wavePacket( p ) : hoState( p );
}




particle::Ptr particle::make(const Pars& p, QM_Picture qmp)
{
  return qmp==QMP_SCH ? Ptr(make_shared<ParticleSch>(p)) : Ptr(make_shared<Particle>(p));
}


particle::PtrPumped particle::makePumped(const ParsPumped& p, QM_Picture qmp)
{
  return qmp==QMP_SCH ? PtrPumped(make_shared<PumpedParticleSch>(p)) : PtrPumped(make_shared<PumpedParticle>(p));
}


particle::Ptr particle::make(const ParsPumped& p, QM_Picture qmp)
{
  if (p.vClass)
    return makePumped(p,qmp);
  else
    return make(static_cast<const Pars&>(p),qmp);

}
