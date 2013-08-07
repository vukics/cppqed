// -*- C++ -*-
#ifndef   ELEMENTS_FREES_MULTILEVEL__H_INCLUDED
#define   ELEMENTS_FREES_MULTILEVEL__H_INCLUDED

#include "MultiLevel_Fwd.h"

#include "ParsMultiLevel.h"

#include "impl/AveragingUtils.tcc"

#include "ElementAveraged.h"
#include "ElementLiouvillean.h"
#include "Free.h"
#include "FreeExact.h"
#include "Hamiltonian.h"

#include <boost/fusion/mpl/size.hpp>

#include <boost/shared_ptr.hpp>


namespace multilevel {


const std::string keyTitle="MultiLevel";


using namespace structure::free;


template<int NL>
struct LevelsMF : mpl::identity<blitz::TinyVector<dcomp,NL> > {};


template<int NL>
struct RealLevelsMF : mpl::identity<blitz::TinyVector<double,NL> > {};



////////
//
// Exact
//
////////

struct MultiLevelExactNotImplementedException : public cpputils::Exception {};

template<int NL>
class Exact : public structure::FreeExact<false>
{
public:
  typedef typename LevelsMF<NL>::type Levels;

  Exact(const Levels& zIs) : FreeExact<false>(NL), zIs_(zIs) {}

  const Levels& get_zIs() const {return zIs_;}

private:
  void updateU(double) const;

  bool isUnitary_v() const;

  const Levels zIs_;

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



template<int N1, int N2>
class Pump : public Storage<dcomp>, public tmptools::pair_c<N1,N2>
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

  typedef typename Exact<NL>::Levels Levels;

  HamiltonianIP(const Levels& zSchs, const Levels& zIs, const VP& etas)
    : Exact<NL>(zIs), zSchs_(zSchs), etas_(etas) {}

private:
  void addContribution_v(double, const StateVectorLow&, StateVectorLow&, double) const;


  const Levels zSchs_;

  const VP etas_;

};


template<int NL, typename VP>
class HamiltonianSch 
  : public structure::Hamiltonian<1,structure::NO_TIME>
{
public:
  static const int NPT=mpl::size<VP>::value; // number of pumped transitions

  typedef typename LevelsMF<NL>::type Levels;

  HamiltonianSch(const Levels& zSchs, const VP& etas) : zSchs_(zSchs), etas_(etas) {}

  const Levels& get_zSchs() const {return zSchs_;}

private:
  void addContribution_v(const StateVectorLow&, StateVectorLow&) const;


  const Levels zSchs_;

  const VP etas_;

};


//////////////
//
// Liouvillean
//
//////////////


template<int N1, int N2>
class Decay : public Storage<double>, public tmptools::pair_c<N1,N2>
{
public:
  typedef Storage<double> Base;

  Decay(double value) : Base(value) {}
  Decay() : Base(0) {}

};


template<int NL, typename VL>
class Liouvillean : public structure::ElementLiouvillean<1,mpl::size<VL>::value>
// Note that, at some point, the Fusion sequence VL needs to be converted into a runtime sequence (JumpStrategies & JumpRateStrategies)
{
public:
  static const int NLT=mpl::size<VL>::value; // number of lossy transitions

  typedef structure::ElementLiouvillean<1,NLT> Base;
  
  typedef typename Base::JumpStrategies     JumpStrategies    ;
  typedef typename Base::JumpRateStrategies JumpRateStrategies;

  BOOST_STATIC_ASSERT( blitzplusplus::TinyVectorLengthTraits<JumpRateStrategies>::value==NLT );

  typedef typename Base::KeyLabels KeyLabels;

  Liouvillean(const VL& gammas) : Base(Liouvillean::fillJS(),Liouvillean::fillJRS(),keyTitle,fillKeyLabels()), gammas_(gammas) {}

  const JumpStrategies     fillJS () const;
  const JumpRateStrategies fillJRS() const;

  static const KeyLabels fillKeyLabels();

private:
  template<int>
  void jumpStrategy(StateVectorLow&) const;

  template<int>
  double jumpRateStrategy(const LazyDensityOperator&) const;

  // NEED_TO_UNDERSTAND can member TEMPLATES be passed as template parameters? This would be needed to fuse fillJS and fillJRS into a template together with the helper classes below
  
  class  JS_helper;
  class JRS_helper;

  class KeyHelper;

  const VL gammas_;

};



///////////
//
// Averaged
//
///////////


class EmptyBase {};


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


template<int NL, typename VP, typename VL, typename AveragingType>
class PumpedLossyMultiLevelSch 
  : public multilevel::HamiltonianSch<NL,VP>,
    public multilevel::Liouvillean<NL,VL>,
    public MultiLevelBase<NL>,
    public AveragingType
  // The ordering becomes important here
{
public:
  // NEEDS_WORK some static sanity checks of VP and VL in view of NL should be done

  typedef typename multilevel::    LevelsMF<NL>::type     Levels;
  typedef typename multilevel::RealLevelsMF<NL>::type RealLevels;

  typedef multilevel::HamiltonianSch<NL,VP> Hamiltonian;
  typedef multilevel::Liouvillean<NL,VL> Liouvillean;
  typedef MultiLevelBase<NL> Base;

  typedef typename Base::Ptr Ptr;

  using Hamiltonian::get_zSchs;
  using Base::getParsStream;

  template<typename... AveragingConstructorParameters>
  PumpedLossyMultiLevelSch(const RealLevels&, const VP&, const VL&, const AveragingConstructorParameters&...);

};



namespace multilevel {

#define RETURN_type typename MultiLevelBase<NL>::Ptr


template<typename AveragingType, int NL, typename VP, typename VL, typename... AveragingConstructorParameters>
inline
RETURN_type
makePumpedLossySch(const blitz::TinyVector<double,NL>& deltas,
                   // note that, understandably, if we write
                   // const typename RealLevelsMF<NL>::type&
                   // here, the compiler cannot deduce NL anymore
                   const VP& etas, const VL& gammas,
                   const AveragingConstructorParameters&... a)
{
  return boost::make_shared<PumpedLossyMultiLevelSch<NL,VP,VL,AveragingType> >(deltas,etas,gammas,a...);
}


template<typename AveragingType, int NL, typename VP, typename VL, typename... AveragingConstructorParameters>
inline
RETURN_type
makePumpedLossySch(const multilevel::ParsPumpedLossy<NL,VP,VL>& p, const AveragingConstructorParameters&... a)
{
  return makePumpedLossySch<AveragingType>(p.deltas,p.etas,p.gammas,a...);
}


template<int NL, typename VP, typename VL>
inline
RETURN_type
makePumpedLossySch(const blitz::TinyVector<double,NL>& deltas, const VP& etas, const VL& gammas, const std::string& keyTitle="PumpedLossyMultiLevelSch", bool offDiagonals=false)
{
  return makePumpedLossySch<ReducedDensityOperator<1> >(deltas,etas,gammas,keyTitle,NL,offDiagonals);
}


template<int NL, typename VP, typename VL>
inline
RETURN_type
makePumpedLossySch(const multilevel::ParsPumpedLossy<NL,VP,VL>& p, const std::string& keyTitle="PumpedLossyMultiLevelSch", bool offDiagonals=false)
{
  return makePumpedLossySch<ReducedDensityOperator<1> >(p,keyTitle,NL,offDiagonals);
}


#undef  RETURN_type

} // multilevel


#endif // ELEMENTS_FREES_MULTILEVEL__H_INCLUDED
