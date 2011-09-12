// -*- C++ -*-
#ifndef   _COMPOSITE_SYSTEM_INCLUDED
#define   _COMPOSITE_SYSTEM_INCLUDED

#include "CompositeFwd.h"

#include "Act.h"

#include "BlitzArraySliceIterator.h"
// This is included at this point mainly to pull in necessary TMP tools

#include "details/TMP_helpers.h"

#include <boost/fusion/container/generation/make_vector.hpp>

#include <boost/preprocessor/iteration/iterate.hpp>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/arithmetic/mul.hpp>


namespace composite {

using boost::fusion::make_vector;

namespace result_of {

// using boost::fusion::result_of::make_list;
using boost::fusion::result_of::make_vector;

} // result_of

} // composite

template<typename VA>
struct MaxRankMF : composite::MaxMF<typename composite::MaxMF<VA,composite::SeqLess<mpl::_1,mpl::_2> >::type>::type {};


template<typename VA>
// VA should model a fusion sequence of Acts
class Composite 
// The base_from_member idiom appears because the Frees has to be calculated and stored somehow first
  : private boost::base_from_member<const blitz::TinyVector<structure::SubSystemFree,MaxRankMF<VA>::type::value+1> >,
    public structure::QuantumSystem<MaxRankMF<VA>::type::value+1>,
    public structure::Exact        <MaxRankMF<VA>::type::value+1>, 
    public structure::Hamiltonian  <MaxRankMF<VA>::type::value+1>,
    public structure::Liouvillean  <MaxRankMF<VA>::type::value+1>,
    public structure::Averaged     <MaxRankMF<VA>::type::value+1>
{
public:
  // The calculated RANK

  static const int RANK=MaxRankMF<VA>::type::value+1;

  // Public types

  typedef structure::QuantumSystem<RANK> QS_Base;
  typedef structure::Exact        <RANK> Ex_Base;
  typedef structure::Hamiltonian  <RANK> Ha_Base;
  typedef structure::Liouvillean  <RANK> Li_Base;
  typedef structure::Averaged     <RANK> Av_Base;

  typedef blitz::TinyVector<structure::SubSystemFree,RANK> Frees;

  typedef boost::base_from_member<const Frees> FreesBase;

  typedef quantumdata::Types<RANK> Types;

  typedef typename Types::    StateVectorLow     StateVectorLow;
  typedef typename Types::DensityOperatorLow DensityOperatorLow;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

  typedef typename tmptools::OrdinalMF<RANK>::type Ordinals;

  typedef typename QS_Base::Dimensions    Dimensions   ;
  typedef typename Li_Base::Probabilities Probabilities;
  typedef typename Av_Base::Averages      Averages     ;

  // Base class names

  using QS_Base::getDimensions; using QS_Base::getTotalDimension; using Li_Base::defaultArray;

  // Compile-time sanity check

  BOOST_MPL_ASSERT_MSG( ( composite::CheckMeta<RANK,VA>::type::value == true ), COMPOSITE_not_CONSISTENT, (mpl::void_) );

private:
  // Constructor helpers
		      
  static const Frees fillFrees(const VA&);

  static const Dimensions fillDimensions(const Frees&);

  template<typename A> // A must be an Act type (which is also a compile-time vector)
  class FurnishedAct : public A
  {
  public:
    typedef FurnishedAct<A> type; // so that it can act as a metafunction
    typedef blitzplusplus::SlicesData<RANK,A> SlicesData;

    FurnishedAct(const A& a, const StateVectorLow& psiProbe) : A(a), sd_(psiProbe) {}

    const SlicesData& getSlicesData() const {return sd_;}

  private:
    const SlicesData sd_;
  };

  typedef typename boost::fusion::result_of::as_vector<typename mpl::transform<VA,FurnishedAct<mpl::_1> >::type>::type VFA;

  static const VFA furnish(const VA&, const Dimensions&);

public:
  // Constructor

  explicit Composite(const VA& acts) 
    : FreesBase(fillFrees(acts)), QS_Base(fillDimensions(FreesBase::member)), frees_(FreesBase::member), furnishedActs_(furnish(acts,getDimensions())) {}


private:
  // Implementing QS_Base

  double highestFrequency (             ) const;

  void   displayParameters(std::ostream&) const; class DisplayParameters;

  // Implementing Ex_Base

  bool isUnitary() const; class IsUnitary;

  void actWithU(double, StateVectorLow&) const; class ActWithU;

  // Implementing Ha_Base

  void addContribution(double, const StateVectorLow&, StateVectorLow&, double) const; class Hamiltonian;

  // Implementing Li_Base

  size_t              nJumps       ()                                   const; class NJumps;
  const Probabilities probabilities(double, const LazyDensityOperator&) const; class Probas;
  void                actWithJ     (double, StateVectorLow&, size_t)    const; class ActWithJ;

  // Implementing Av_Base

  void   displayKey(std::ostream&, size_t&) const; class DisplayKey;
  size_t nAvr      ()                       const; class NAvr;

  const Averages average(double, const LazyDensityOperator&)  const; class Average;
  void           process(Averages&)                           const; class Process;
  void           display(const Averages&, std::ostream&, int) const; class Display;

  // overall helpers

  template<typename H>
  void worker(const H&) const;

  // Storage
  // Note that Frees are stored by value in FreesBase

  const Frees& frees_;
  // const VA      acts_;

  const VFA furnishedActs_;

};


// The following provides a much more convenient interface:


#define BOOST_PP_ITERATION_LIMITS (1,10)
#define BOOST_PP_FILENAME_1 "details/CompositeMakerImplementationsSpecialization.h"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS


#include "impl/Composite.tcc"

#endif // _COMPOSITE_SYSTEM_INCLUDED
