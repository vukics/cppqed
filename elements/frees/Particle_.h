// -*- C++ -*-
#ifndef   _PARTICLE__H
#define   _PARTICLE__H

#include "Particle_Fwd.h"

#include "ParsParticle.h"

#include "QM_PictureFwd.h"
#include "StateVectorFwd.h"

#include "ModeFunction.h"
#include "ElementAveraged.h"
#include "Free.h"
#include "FreeExact.h"
#include "TridiagonalHamiltonian.h"

#include "FFTFwd.h"

#include <boost/shared_ptr.hpp>

namespace particle {

using namespace structure::free;


typedef boost::shared_ptr<const ParticleBase> SmartPtr;

typedef boost::shared_ptr<const PumpedParticleBase> SmartPtrPumped;


const Tridiagonal expINKX(const ParticleBase*, ptrdiff_t);

inline const Tridiagonal sinNKX(const ParticleBase* particle, ptrdiff_t nK) {return (expINKX(particle,nK)-expINKX(particle,-nK))/(2.*DCOMP_I);}
inline const Tridiagonal cosNKX(const ParticleBase* particle, ptrdiff_t nK) {return (expINKX(particle,nK)+expINKX(particle,-nK))/ 2.         ;}

const Tridiagonal mfNKX       (const ParticleBase*, const ModeFunction&);
const Tridiagonal mfNKX_AbsSqr(const ParticleBase*, const ModeFunction&);


const StateVector wavePacket(const InitialCondition&, const Spatial&, bool kFlag=true);
const StateVector wavePacket(const Pars      &,                       bool kFlag=true);
const StateVector wavePacket(const ParsPumped&,                       bool kFlag=true);

const StateVector hoState(size_t n, const InitialCondition&, const Spatial&, bool kFlag=true/*, bool exactRenorm=true*/);
const StateVector hoState(const Pars      &,                                 bool Kflag=true);
const StateVector hoState(const ParsPumped&,                                 bool Kflag=true);

const StateVector init(const Pars&);


void ffTransform(StateVectorLow&, fft::Direction);

SmartPtr make(const Pars      &, QM_Picture);
SmartPtr make(const ParsPumped&, QM_Picture);



namespace details {

typedef boost::tuple<const Spatial&, double> Storage;

class Empty {};

} // details


class Exact : public structure::FreeExact, private details::Storage
{
public:
  Exact(const Spatial&, double omrec);

  using details::Storage::get;

private:
  void updateU(double) const;

  bool isUnitary() const {return true;}

  const Factors factorExponents_;

};


template<bool IS_TD>
class Hamiltonian 
  : public structure::TridiagonalHamiltonian<1,IS_TD>,
    public mpl::if_c<IS_TD,Exact,details::Empty>::type
{
public:
  typedef structure::TridiagonalHamiltonian<1,IS_TD> Base;

  Hamiltonian(const Spatial&, double omrec, double vClass, const ModeFunction&);
  Hamiltonian(const Spatial&, double omrec, mpl::bool_<IS_TD> =mpl::false_());

};



// Spatial is a tool to facilitate state vector representations in both X and K space, where the two are canonically conjugate operators, so that [X,K]=i, and hence the two representations are linked with ffTransform.

class Spatial
{
public:
  typedef TTD_DARRAY(1) Array;

  explicit Spatial(size_t, double deltaK=1);

  size_t getFinesse  () const {return fin_;}
  size_t getDimension() const {return x_.size();}

  // x and k space coordinates, respectively
  double x(size_t i) const {return x_(i);}
  double k(size_t i) const {return k_(i);} 

  const Array& getX() const {return x_;}
  const Array& getK() const {return k_;}

  void header(std::ostream&) const;

private:
  const size_t fin_; // = log2 (number of dimensions)

  const double xMax_, deltaX_, kMax_, deltaK_;

  const Array x_, k_; 

};


class Averaged
  : public structure::ElementAveraged<1>
{
public:
  typedef structure::ElementAveraged<1> Base;

  Averaged(const Spatial&);

  const Spatial& getSpace() const {return space_;}

private:
  const Averages average(const LazyDensityOperator&) const;
  void           process(Averages&)                  const;

  const Spatial space_;

};



} // particle



////////////////
//
// Highest level
//
////////////////


class ParticleBase
  : public structure::Free, public particle::Averaged
{
protected:
  explicit ParticleBase(size_t fin, 
			const RealFreqs& =RealFreqs(), const ComplexFreqs& =ComplexFreqs());

};


class PumpedParticleBase
  : public ParticleBase
{
public:
  double getV_Class() const {return vClass_;}

  const ModeFunction& getMF() const {return mf_;}

protected:
  PumpedParticleBase(size_t fin, double vClass, const ModeFunction&,
		     const RealFreqs& =RealFreqs(), const ComplexFreqs& =ComplexFreqs());

private:
  const double       vClass_;
  const ModeFunction mf_    ;

};


class Particle
  : public ParticleBase, public particle::Exact
{
public:
  explicit Particle(const particle::Pars&);

};


class ParticleSch
  : public ParticleBase, public particle::Hamiltonian<false>
{
public:
  explicit ParticleSch(const particle::Pars&);

};


class PumpedParticle
  : public PumpedParticleBase, public particle::Hamiltonian<true>
{
public:
  explicit PumpedParticle(const particle::ParsPumped&);
};


class PumpedParticleSch
  : public PumpedParticleBase, public particle::Hamiltonian<false>
{
public:
  explicit PumpedParticleSch(const particle::ParsPumped&);
};


#endif // _PARTICLE__H
