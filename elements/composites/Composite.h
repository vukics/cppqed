// -*- C++ -*-
#ifndef   ELEMENTS_COMPOSITES_COMPOSITE_H_INCLUDED
#define   ELEMENTS_COMPOSITES_COMPOSITE_H_INCLUDED

#include "CompositeFwd.h"

#include "Act.h"

#include "BlitzArraySliceIterator.h"
// This is included at this point mainly to pull in necessary TMP tools

#include "SubSystemFwd.h"

#include "details/TMP_helpers.h"

#include <boost/fusion/container/generation/make_list.hpp>



namespace composite {

using boost::fusion::make_list;

using ::size_t;

namespace result_of {

using boost::fusion::result_of::make_list;

} // result_of


template<typename VA>
struct MaxRank : MaxMF<typename MaxMF<VA,composite::SeqLess<mpl::_1,mpl::_2> >::type>::type::type {};



template<typename VA>
// VA should model a fusion sequence of Acts
class Base
// The base_from_member idiom appears because the Frees has to be calculated and stored somehow first
  : public structure::QuantumSystem<composite::MaxRank<VA>::value+1>,
    public structure::Averaged     <composite::MaxRank<VA>::value+1>
{
public:
  // The calculated RANK
  static const int RANK=composite::MaxRank<VA>::value+1;

  // Public types

  typedef structure::QuantumSystem<RANK> QS_Base;
  typedef structure::Averaged     <RANK> Av_Base;

  typedef blitz::TinyVector<composite::SubSystemFree,RANK> Frees;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

  typedef tmptools::Ordinals<RANK> Ordinals;

  typedef typename QS_Base::Dimensions Dimensions;
  typedef typename Av_Base::Averages   Averages  ;

  // Base class names

  using QS_Base::getDimensions; using QS_Base::getTotalDimension;

private:
  // Constructor helper
  static const Dimensions fillDimensions(const Frees&);

public:
  // Constructor
  explicit Base(const Frees& frees, const VA& acts) 
    : QS_Base(fillDimensions(frees)), frees_(frees), acts_(acts) {}

private:
  // Implementing QS_Base

  double  highestFrequency_v(             ) const;
  void   displayParameters_v(std::ostream&) const; class DisplayParameters;

  // Implementing Av_Base

  void   displayKey_v(std::ostream&, size_t&) const; class DisplayKey;
  size_t       nAvr_v()                       const; class NAvr;

  const Averages average_v(double, const LazyDensityOperator&)  const; class Average;
  void           process_v(Averages&)                           const; class Process;
  void           display_v(const Averages&, std::ostream&, int) const; class Display;

  const Frees& frees_;
  const VA   &  acts_;

};




template<typename VA>
class Exact
  : public structure::Exact<MaxRank<VA>::value+1>
{
private:
  static const int RANK=MaxRank<VA>::value+1;

  typedef blitz::TinyVector<SubSystemFree,RANK> Frees;

  typedef typename quantumdata::Types<RANK>::StateVectorLow StateVectorLow;

  typedef tmptools::Ordinals<RANK> Ordinals;

protected:
  Exact(const Frees& frees, const VA& acts) : frees_(frees), acts_(acts) {}

private:
  bool isUnitary_v(                       ) const; class IsUnitary;
  void  actWithU_v(double, StateVectorLow&) const; class ActWithU ;

  const Frees& frees_;
  const VA   &  acts_;

};


template<typename VA>
class Hamiltonian
  : public structure::Hamiltonian<MaxRank<VA>::value+1>
{
private:
  static const int RANK=MaxRank<VA>::value+1;

  typedef blitz::TinyVector<SubSystemFree,RANK> Frees;

  typedef typename quantumdata::Types<RANK>::StateVectorLow StateVectorLow;

  typedef tmptools::Ordinals<RANK> Ordinals;

protected:
  Hamiltonian(const Frees& frees, const VA& acts) : frees_(frees), acts_(acts) {}

private:
  void addContribution_v(double, const StateVectorLow&, StateVectorLow&, double) const; class AddContribution;

  const Frees& frees_;
  const VA   &  acts_;

};


template<typename VA>
class Liouvillean
  : public structure::Liouvillean<MaxRank<VA>::value+1>
{
private:
  static const int RANK=MaxRank<VA>::value+1;

  typedef blitz::TinyVector<SubSystemFree,RANK> Frees;

  typedef typename quantumdata::Types<RANK>::StateVectorLow StateVectorLow;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

  typedef tmptools::Ordinals<RANK> Ordinals;

  typedef typename structure::Liouvillean<RANK>::Probabilities Probabilities;

protected:
  Liouvillean(const Frees& frees, const VA& acts) : frees_(frees), acts_(acts) {}

private:
  size_t                     nJumps_v()                                   const; class NJumps;
  const Probabilities probabilities_v(double, const LazyDensityOperator&) const; class Probas;
  void                     actWithJ_v(double, StateVectorLow&, size_t)    const; class ActWithJ;

  void displayKey_v(std::ostream&, size_t&) const; // {}

  const Frees& frees_;
  const VA   &  acts_;

};


template<typename VA>
class EmptyBase
{
public:
  typedef blitz::TinyVector<SubSystemFree,MaxRank<VA>::value+1> Frees;

  EmptyBase(const Frees&, const VA&) {}
  
};


} // composite



#define BASE_CLASS(Aux,Class) mpl::if_c<IS_##Aux,composite::Class<VA>,composite::EmptyBase<VA> >::type

template<typename VA, bool IS_EX=true, bool IS_HA=true, bool IS_LI=true>
// VA should model a fusion sequence of Acts
class Composite 
// The base_from_member idiom appears because the Frees has to be calculated and stored somehow first
  : private boost::base_from_member<const blitz::TinyVector<composite::SubSystemFree,composite::MaxRank<VA>::value+1> >,
    public composite::Base<VA>,
    public BASE_CLASS(EX,Exact),
    public BASE_CLASS(HA,Hamiltonian),
    public BASE_CLASS(LI,Liouvillean)
{
public:
  typedef typename BASE_CLASS(EX,Exact)             ExactBase;
  typedef typename BASE_CLASS(HA,Hamiltonian) HamiltonianBase;
  typedef typename BASE_CLASS(LI,Liouvillean) LiouvilleanBase;
  
  // The calculated RANK
  static const int RANK=composite::MaxRank<VA>::value+1;

  // Public types

  typedef blitz::TinyVector<composite::SubSystemFree,RANK> Frees;

  typedef boost::base_from_member<const Frees> FreesBase;

  // Compile-time sanity check
  BOOST_MPL_ASSERT_MSG( ( composite::CheckMeta<RANK,VA>::type::value == true ), COMPOSITE_not_CONSISTENT, (mpl::void_) );

private:
  // Constructor helper
  static const Frees fillFrees(const VA&);

public:
  // Constructor
  explicit Composite(const VA& acts)
    : FreesBase(fillFrees(acts)), composite::Base<VA>(FreesBase::member,acts),
      ExactBase      (FreesBase::member,acts),
      HamiltonianBase(FreesBase::member,acts),
      LiouvilleanBase(FreesBase::member,acts),
      acts_(acts) {}

private:
  // Storage
  // Note that Frees are stored by value in FreesBase, but the Acts need to be stored by value here:
  const VA acts_;

};

#undef BASE_CLASS


// The following provides a much more convenient interface:

namespace composite {

namespace result_of {


namespace mpl=boost::mpl;

typedef Act<> DefaultArgument;


template<BOOST_PP_ENUM_BINARY_PARAMS(FUSION_MAX_VECTOR_SIZE,typename A,=DefaultArgument BOOST_PP_INTERCEPT)> 
struct Make : boost::mpl::identity<Composite<typename make_list<BOOST_PP_ENUM_PARAMS(FUSION_MAX_VECTOR_SIZE,A)>::type> >
{};


} // result_of

} // composite


#define DEFAULT_print(z, n, data) DefaultArgument

#define BOOST_PP_ITERATION_LIMITS (1,BOOST_PP_SUB(FUSION_MAX_VECTOR_SIZE,1) )
#define BOOST_PP_FILENAME_1 "../composites/details/CompositeMakerImplementationsSpecialization.h"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS

#undef DEFAULT_print


#endif // ELEMENTS_COMPOSITES_COMPOSITE_H_INCLUDED
