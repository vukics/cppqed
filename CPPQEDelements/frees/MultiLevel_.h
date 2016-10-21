/// \briefFile{Defines free elements of the \ref multilevelbundle "MultiLevel bundle" }
// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
#ifndef   CPPQEDELEMENTS_FREES_MULTILEVEL__H_INCLUDED
#define   CPPQEDELEMENTS_FREES_MULTILEVEL__H_INCLUDED

#include "MultiLevel_Fwd.h"

#include "ParsMultiLevel.h"

#include "AveragingUtils.tcc"

#include "ElementAveraged.h"
#include "ElementLiouvillean.h"
#include "Free.h"
#include "FreeExact.h"
#include "Hamiltonian.h"

#include <boost/fusion/mpl/size.hpp>

#include <boost/shared_ptr.hpp>


/// Contains helpers for the \ref multilevelbundle "MultiLevel bundle"
namespace multilevel {


const std::string keyTitle="MultiLevel";


using namespace structure::freesystem;

/// Type for storing complex level energies (the \f$z_i\f$s \ref multilevelelements "here") \tparam NL number of levels
template<int NL> using Levels = blitz::TinyVector<dcomp,NL>;

/// Type for storing level energies (the \f$\delta_i\f$s \ref multilevelelements "here") \tparam NL number of levels
template<int NL> using RealLevels = blitz::TinyVector<double,NL>;



////////
//
// Exact
//
////////

struct MultiLevelExactNotImplementedException : public cpputils::Exception {};

/** \todo MultiLevel Exact not fully implemented, cf. the late 'shift' function. Maybe it's easier with something similar to Tridiagonal::freqs_ in Sigma */
template<int NL>
class Exact : public structure::FreeExact<false>
{
public:
  typedef Levels<NL> L;

  Exact(const L& zIs) : FreeExact<false>(NL), zIs_(zIs) {}

  const L& get_zIs() const {return zIs_;}

private:
  void updateU(double) const;

  bool applicableInMaster_v() const;

  const L zIs_;

};



//////////////
//
// Hamiltonian
//
//////////////


template<typename T>
class Storage
{
public:
  Storage(const T& value) : value_(value) {}

  const T& get() const {return value_;}
  void set(const T& value) {value_=value;}

private:
  T value_;

};


template<>
class Storage<double>
{
public:
  Storage(double value) : value_(value) {}

  double get() const {return value_;}
  void set(double value) {value_=value;}

private:
  double value_;

};



template<typename T>
inline
std::ostream& operator<<(std::ostream& os, const Storage<T>& s)
{
  return os<<s.get();
}


template<typename T>
inline
std::istream& operator>>(std::istream& is,       Storage<T>& s)
{
  T temp; is>>temp;
  s.set(temp);
  return is;
}


/// Class representing an elementary pump term (an \f$\eta_{ij}\f$ \ref multilevelactualHamiltonian "here") with a compile-time pair \f$i,j\f$ and a runtime complex value
template<int I, int J>
class Pump : public Storage<dcomp>, public tmptools::pair_c<I,J>
{
public:
  typedef Storage<dcomp> Base;

  Pump(const dcomp& value) : Base(value) {}
  Pump() : Base(dcomp()) {}

};


template<int NL, typename VP>
class HamiltonianIP 
  : public structure::Hamiltonian<1>,
    public Exact<NL>
{
public:
  static const int NPT=mpl::size<VP>::value; // number of pumped transitions

  typedef typename Exact<NL>::L L;

  HamiltonianIP(const L& zSchs, const L& zIs, const VP& etas)
    : Exact<NL>(zIs), zSchs_(zSchs), etas_(etas) {}

private:
  void addContribution_v(double, const StateVectorLow&, StateVectorLow&, double) const;


  const L zSchs_;

  const VP etas_;

};


template<int NL, typename VP>
class HamiltonianSch 
  : public structure::HamiltonianTimeDependenceDispatched<1,structure::NO_TIME>
{
public:
  static const int NPT=mpl::size<VP>::value; // number of pumped transitions

  typedef Levels<NL> L;

  HamiltonianSch(const L& zSchs, const VP& etas) : zSchs_(zSchs), etas_(etas) {}

  const L& get_zSchs() const {return zSchs_;}

private:
  void addContribution_v(structure::NoTime, const StateVectorLow&, StateVectorLow&) const;


  const L zSchs_;

  const VP etas_;

};


//////////////
//
// Liouvillean
//
//////////////


/// Class representing an elementary decay term (a \f$\gamma_{ij}\f$ \ref multilevelelements "here") with a compile-time pair \f$i,j\f$ and a runtime real value
template<int I, int J>
class Decay : public Storage<double>, public tmptools::pair_c<I,J>
{
public:
  typedef Storage<double> Base;

  Decay(double value) : Base(value) {}
  Decay() : Base(0) {}

};


template<int NL, typename VL>
class Liouvillean : public structure::ElementLiouvilleanStrategies<1,mpl::size<VL>::value+NL>
// Note that, at some point, the Fusion sequence VL needs to be converted into a runtime sequence (JumpStrategies & JumpRateStrategies)
{
public:
  static const int NLT=mpl::size<VL>::value; // number of lossy transitions

  typedef structure::ElementLiouvilleanStrategies<1,NLT+NL> Base;
  
  typedef typename Base::JumpStrategies     JumpStrategies    ;
  typedef typename Base::JumpRateStrategies JumpRateStrategies;

  static_assert( blitzplusplus::TinyVectorLengthTraits<JumpRateStrategies>::value==NLT+NL , "Jump number inconsistent." );

  typedef typename Base::KeyLabels KeyLabels;

  Liouvillean(const VL& gammas, double gamma_parallel) : Base(Liouvillean::fillJS(),Liouvillean::fillJRS(),keyTitle,fillKeyLabels()), gammas_(gammas), gamma_parallel_(gamma_parallel) {}

  const JumpStrategies     fillJS () const;
  const JumpRateStrategies fillJRS() const;

  static const KeyLabels fillKeyLabels();

private:
  template<int>
  void jumpStrategy(StateVectorLow&) const;

  template<int>
  double jumpRateStrategy(const LazyDensityOperator&) const;

  // NEED_TO_UNDERSTAND can member TEMPLATES be passed as template parameters? This would be needed to fuse fillJS and fillJRS into a template together with the helper classes below

  void   flipStrategy(StateVectorLow& psi, size_t i) const {psi*=sqrt(2.*gamma_parallel_); psi(i)*=-1.;}
  double flipRateStrategy(const LazyDensityOperator&) const {return -1;} // Being a member is somewhat superfluous here

  class  JS_helper;
  class JRS_helper;

  class KeyHelper;

  const VL gammas_;

  const double gamma_parallel_;

};


} // multilevel


////////////////
//
// Highest level
//
////////////////


template<int NL>
class MultiLevelBase 
  : public structure::Free
{
public:
  typedef boost::shared_ptr<const MultiLevelBase> Ptr;

  using structure::Free::getParsStream;

  MultiLevelBase(const RealFreqs& realFreqs=RealFreqs(), const ComplexFreqs& complexFreqs=ComplexFreqs())
    : structure::Free(NL,realFreqs,complexFreqs)
  {
    getParsStream()<<"# "<<multilevel::keyTitle<<std::endl;
  }

  virtual ~MultiLevelBase() {}

};


/// Implements a free multi-level system with driving and loss \see \ref multilevelbundle
/**
 * \tparam NL number of levels
 * \tparam VP fusion vector of multilevel::Pump types representing all the pump terms (the sequence of \f$\eta_{ij}\f$s \ref multilevelactualHamiltonian "here") in the system
 * \tparam VL fusion vector of multilevel::Decay types representing all the decays (the sequence of \f$\gamma_{ij}\f$s \ref multilevelelements "here") of the system
 *
 * \todo some static sanity checks of `VP` and `VL` in view of `NL` should be done
 */
template<int NL, typename VP, typename VL, typename AveragingType>
class PumpedLossyMultiLevelSch 
  : public multilevel::HamiltonianSch<NL,VP>,
    public multilevel::Liouvillean<NL,VL>,
    public MultiLevelBase<NL>,
    public AveragingType
  // The ordering becomes important here
{
public:
  typedef multilevel::    Levels<NL>     Levels;
  typedef multilevel::RealLevels<NL> RealLevels;

  typedef multilevel::HamiltonianSch<NL,VP> Hamiltonian;
  typedef multilevel::Liouvillean   <NL,VL> Liouvillean;
  typedef MultiLevelBase<NL> Base;

  typedef typename Base::Ptr Ptr;

  template<typename... AveragingConstructorParameters>
  PumpedLossyMultiLevelSch(const RealLevels&, const VP&, const VL&, double gamma_parallel, AveragingConstructorParameters&&...);

};



namespace multilevel {

#define RETURN_type typename MultiLevelBase<NL>::Ptr

/// Maker function for PumpedLossyMultiLevelSch
template<typename AveragingType, int NL, typename VP, typename VL, typename... AveragingConstructorParameters>
inline
RETURN_type
makePumpedLossySch(const RealLevels<NL>& deltas,
                   const VP& etas, const VL& gammas, double gamma_parallel,
                   AveragingConstructorParameters&&... a)
{
  return boost::make_shared<PumpedLossyMultiLevelSch<NL,VP,VL,AveragingType> >(deltas,etas,gammas,gamma_parallel,a...);
}

/// \overload
template<typename AveragingType, int NL, typename VP, typename VL, typename... AveragingConstructorParameters>
inline
RETURN_type
makePumpedLossySch(const multilevel::ParsPumpedLossy<NL,VP,VL>& p, AveragingConstructorParameters&&... a)
{
  return makePumpedLossySch<AveragingType>(p.deltas,p.etas,p.gammas,p.gamma_parallel,a...);
}

/// \overload
template<int NL, typename VP, typename VL>
inline
RETURN_type
makePumpedLossySch(const RealLevels<NL>& deltas, const VP& etas, const VL& gammas, double gamma_parallel, const std::string& keyTitle="PumpedLossyMultiLevelSch", bool offDiagonals=false)
{
  return makePumpedLossySch<ReducedDensityOperator<1> >(deltas,etas,gammas,gamma_parallel,keyTitle,NL,offDiagonals);
}

/// \overload
template<int NL, typename VP, typename VL>
inline
RETURN_type
makePumpedLossySch(const multilevel::ParsPumpedLossy<NL,VP,VL>& p, const std::string& keyTitle="PumpedLossyMultiLevelSch", bool offDiagonals=false)
{
  return makePumpedLossySch<ReducedDensityOperator<1> >(p,keyTitle,NL,offDiagonals);
}


#undef  RETURN_type

} // multilevel


#endif // CPPQEDELEMENTS_FREES_MULTILEVEL__H_INCLUDED
