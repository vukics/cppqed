// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines the particle-bundle}
#ifndef   CPPQEDELEMENTS_FREES_PARTICLE__H_INCLUDED
#define   CPPQEDELEMENTS_FREES_PARTICLE__H_INCLUDED

#include "Particle_Fwd.h"

#include "ParsParticle.h"

#include "QM_PictureFwd.h"
#include "StateVectorFwd.h"

#include "ModeFunction.h"
#include "ElementAveraged.h"
#include "Exception.h"
#include "Free.h"
#include "FreeExact.h"
#include "TridiagonalHamiltonian.h"

#include <tuple>

namespace particle {

using namespace structure::freesystem; using structure::NoTime;


typedef std::shared_ptr<const ParticleBase> Ptr;

typedef std::shared_ptr<const PumpedParticleBase> PtrPumped;

struct NotATridiagonal : public cpputils::Exception {};

const Tridiagonal expINKX(particle::Ptr, ptrdiff_t);

inline const Tridiagonal sinNKX(particle::Ptr particle, ptrdiff_t nK) {return (expINKX(particle,nK)-expINKX(particle,-nK))/(2.*DCOMP_I);}
inline const Tridiagonal cosNKX(particle::Ptr particle, ptrdiff_t nK) {return (expINKX(particle,nK)+expINKX(particle,-nK))/ 2.         ;}

const Tridiagonal mfNKX       (particle::Ptr, const ModeFunction&);
const Tridiagonal mfNKX_AbsSqr(particle::Ptr, const ModeFunction&);

/// Returns the composition of two mode functions for the special case that this composition is a Tridiagonal.
/**
 * Given two mode functions \f$f(nk_1 x)\f$ and \f$g(nk_2 x)\f$, this functions returns the composition
 * \f[f^\ast(nk_1 x)g(nk_2 x)\f]
 * if this result can be represented as a Tridiagonal, i.e. if \f$|nk_1|=|nk_2|\f$. Otherwise the exception NotATridiagonal is thrown.
 */
const Tridiagonal mfComposition(particle::Ptr particle,                 ///< Determines the dimension of the result.
                                const ModeFunction& modeFunction1,      ///< \f$f(nkx)\f$ above
                                const ModeFunction& modeFunction2);     ///< \f$g(nkx)\f$ above


StateVector wavePacket(const InitialCondition&, const Spatial&, bool kFlag=true);
StateVector wavePacket(const Pars      &,                       bool kFlag=true);
StateVector wavePacket(const ParsPumped&,                       bool kFlag=true);

StateVector hoState(int n, const InitialCondition&, const Spatial&, bool kFlag=true/*, bool exactRenorm=true*/);
StateVector hoState(const Pars      &,                                 bool Kflag=true);
StateVector hoState(const ParsPumped&,                                 bool Kflag=true);

StateVector init(const Pars&);

Ptr make(const Pars&, QM_Picture);

PtrPumped makePumped(const ParsPumped&, QM_Picture);


namespace details {

typedef std::tuple<const Spatial&, double> Storage;

} // details


class Exact : public structure::FreeExact<false>, public details::Storage
{
public:
  Exact(const Spatial&, double omrec);

private:
  void updateU(structure::OneTime) const;

  bool applicableInMaster_v() const {return true;}

  const Diagonal factorExponents_;

};


namespace details { struct EmptyBase {}; }


template<bool IS_TIME_DEPENDENT>
class Hamiltonian 
  : public quantumoperator::TridiagonalHamiltonian<1,IS_TIME_DEPENDENT>,
    public mpl::if_c<IS_TIME_DEPENDENT,Exact,details::EmptyBase>::type
{
public:
  typedef quantumoperator::TridiagonalHamiltonian<1,IS_TIME_DEPENDENT> Base;

  Hamiltonian(const Spatial&, double omrec, double vClass, const ModeFunction&);
  Hamiltonian(const Spatial&, double omrec, mpl::bool_<IS_TIME_DEPENDENT> =mpl::false_());

};


// Spatial is a tool to facilitate state vector representations in both X and K space, where the two are canonically conjugate operators, so that [X,K]=i, and hence the two representations are linked with ffTransform.

class Spatial
{
public:
  typedef DArray<1> Array;

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
  const Averages average_v(NoTime, const LazyDensityOperator&) const;
  void           process_v(        Averages&)                  const;

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
  explicit ParticleBase(size_t fin, const RealFreqs& =emptyRF, const ComplexFreqs& =emptyCF);

};


class PumpedParticleBase
  : public ParticleBase
{
public:
  double getV_Class() const {return vClass_;}

  const ModeFunction& getMF() const {return mf_;}

protected:
  PumpedParticleBase(size_t fin, double vClass, const ModeFunction&,
                     const RealFreqs& =emptyRF, const ComplexFreqs& =emptyCF);

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


#endif // CPPQEDELEMENTS_FREES_PARTICLE__H_INCLUDED
