// -*- C++ -*-
#if !BOOST_PP_IS_ITERATING

#ifndef   COMPOSITES_COMPOSITE_H_INCLUDED
#define   COMPOSITES_COMPOSITE_H_INCLUDED

#include "CompositeFwd.h"

#include "Act.h"

#include "BlitzArraySliceIterator.h"
// This is included at this point mainly to pull in necessary TMP tools

#include "SubSystem.h"

#include "details_TMP_helpers.h"

#include <boost/fusion/container/generation/make_list.hpp>

#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_trailing.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>


namespace composite {

using boost::fusion::make_list;

using ::size_t;

namespace result_of {

using boost::fusion::result_of::make_list;

} // result_of


template<typename VA>
struct MaxRank : MaxMF<typename MaxMF<VA,SeqLess<mpl::_1,mpl::_2> >::type>::type::type {};



template<int N_RANK>
// Factoring out code that depends only on RANK:
class RankedBase : public structure::QuantumSystem<N_RANK>
{
public:
  static const int RANK=N_RANK;

  typedef boost::shared_ptr<const RankedBase<RANK> > Ptr;

  typedef blitz::TinyVector<SubSystemFree,RANK> Frees;

  typedef structure::QuantumSystem<RANK> QS_Base;

  typedef typename QS_Base::Dimensions Dimensions;

  typedef tmptools::Ordinals<RANK> Ordinals;

private:
  // Constructor helper
  static const Dimensions fillDimensions(const Frees&);

protected:
  explicit RankedBase(const Frees& frees)
    : QS_Base(fillDimensions(frees)), frees_(frees) {}
  
  const Frees& getFrees() const {return frees_;}

private:
  const Frees frees_;
 
};


template<typename VA>
// VA should model a fusion sequence of Acts
class Base
  : public RankedBase<MaxRank<VA>::value+1>,
    public structure::Averaged<MaxRank<VA>::value+1>
{
public:
  typedef boost::shared_ptr<const Base<VA> > Ptr;
  
  // The calculated RANK
  static const int RANK=MaxRank<VA>::value+1;

  // Public types
  typedef VA Acts;

  typedef RankedBase<RANK>          RBase  ;
  typedef structure::Averaged<RANK> Av_Base;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

  typedef typename Av_Base::Averages Averages  ;

  typedef typename RBase::   Frees    Frees;
  typedef typename RBase::Ordinals Ordinals;
  
  template<structure::LiouvilleanAveragedTag>
  static std::ostream& displayKeyLA(std::ostream& , size_t&, const Frees&, const VA& acts);

  template<structure::LiouvilleanAveragedTag>
  class DisplayKey;
  
  template<structure::LiouvilleanAveragedTag>
  static size_t nAvrLA(const Frees& frees, const VA& acts);

  template<structure::LiouvilleanAveragedTag>
  class NAvr;
 
  template<structure::LiouvilleanAveragedTag>
  static const Averages averageLA(double t, const LazyDensityOperator& ldo, const Frees& frees, const VA& acts, size_t numberAvr);

  template<structure::LiouvilleanAveragedTag>
  class Average;

protected:
  // Constructor
  explicit Base(const Frees& frees, const VA& acts)
    : RBase(frees), frees_(RBase::getFrees()), acts_(acts) {}
    
  const VA& getActs() const {return acts_;}

private:
  // Implementing QuantumSystem interface

  double         highestFrequency_v(             ) const;
  std::ostream& displayParameters_v(std::ostream&) const; class DisplayParameters;

  // Implementing Av_Base

  std::ostream& displayKey_v(std::ostream& os, size_t& i)              const {return displayKeyLA<structure::LA_Av>(os,i, frees_,acts_         );}
  size_t              nAvr_v()                                         const {return       nAvrLA<structure::LA_Av>(      frees_,acts_         );}
  const Averages   average_v(double t, const LazyDensityOperator& ldo) const {return    averageLA<structure::LA_Av>(t,ldo,frees_,acts_,nAvr_v());}
  
  void          process_v(Averages&)                           const; class Process;
  std::ostream& display_v(const Averages&, std::ostream&, int) const; class Display;

  const Frees& frees_;
  const VA      acts_;

};


// Constructor helper
template<typename VA>
const typename Base<VA>::Frees fillFrees(const VA& acts);



template<typename VA>
const typename Base<VA>::Ptr doMake(const VA&);



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
  bool applicableInMaster_v(                     ) const; class IsUnitary;
  void actWithU_v(double, StateVectorLow&, double) const; class ActWithU ;

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

  typedef typename structure::Liouvillean<RANK>::Rates Rates;

protected:
  Liouvillean(const Frees& frees, const VA& acts) : frees_(frees), acts_(acts) {}

private:
  std::ostream& displayKey_v(std::ostream& os, size_t& i)              const {return Base<VA>::template displayKeyLA<structure::LA_Li>(os,i, frees_,acts_         );}
  size_t              nAvr_v()                                         const {return Base<VA>::template       nAvrLA<structure::LA_Li>(      frees_,acts_         );}
  const Rates      average_v(double t, const LazyDensityOperator& ldo) const {return Base<VA>::template    averageLA<structure::LA_Li>(t,ldo,frees_,acts_,nAvr_v());}

  void actWithJ_v(double, StateVectorLow&, size_t) const; class ActWithJ;

  const Frees& frees_;
  const VA   &  acts_;

};


template<typename>
class EmptyBase
{
public:
  template<typename VA>
  EmptyBase(const blitz::TinyVector<SubSystemFree,MaxRank<VA>::value+1>&, const VA&) {}
  
};


} // composite



#define BASE_class(Aux,Class) mpl::if_c<IS_##Aux,composite::Class<VA>,composite::EmptyBase<composite::Class<VA> > >::type

template<typename VA, bool IS_EX=true, bool IS_HA=true, bool IS_LI=true>
// VA should model a fusion sequence of Acts
class Composite
  : public composite::Base<VA>,
    public BASE_class(EX,Exact),
    public BASE_class(HA,Hamiltonian),
    public BASE_class(LI,Liouvillean)
{
public:
  typedef composite::Base<VA> Base;
  
  typedef typename BASE_class(EX,Exact)             ExactBase;
  typedef typename BASE_class(HA,Hamiltonian) HamiltonianBase;
  typedef typename BASE_class(LI,Liouvillean) LiouvilleanBase;
  
  typedef typename composite::Base<VA>::Frees Frees;
  
  // The calculated RANK
  static const int RANK=composite::MaxRank<VA>::value+1;

  // Compile-time sanity check
  BOOST_MPL_ASSERT_MSG( ( composite::CheckMeta<RANK,VA>::type::value == true ), COMPOSITE_not_CONSISTENT, (composite::CheckMeta<RANK,VA>) );

private:
  using Base::getFrees; using Base::getActs ;
  
public:
  // Constructor
  explicit Composite(const VA& acts)
    : Base(composite::fillFrees(acts),acts),
      ExactBase      (getFrees(),getActs()),
      HamiltonianBase(getFrees(),getActs()),
      LiouvilleanBase(getFrees(),getActs()) {}
      
// private:
  Composite(const Frees& frees, const VA& acts)
    : Base(frees ,acts),
      ExactBase      (getFrees(),getActs()),
      HamiltonianBase(getFrees(),getActs()),
      LiouvilleanBase(getFrees(),getActs()) {}
  
  friend const typename composite::Base<VA>::Ptr composite::doMake<VA>(const VA&);

  // Note that Frees and Acts are stored by value in Base

};

#undef BASE_class


// The following provides a much more convenient interface:

namespace composite {

namespace result_of {


namespace mpl=boost::mpl;

typedef _<> DefaultArgument;


template<bool IS_EX, bool IS_HA, bool IS_LI, BOOST_PP_ENUM_BINARY_PARAMS(FUSION_MAX_VECTOR_SIZE,typename A,=DefaultArgument BOOST_PP_INTERCEPT)> 
struct MakeConcrete : boost::mpl::identity<Composite<typename make_list<BOOST_PP_ENUM_PARAMS(FUSION_MAX_VECTOR_SIZE,A)>::type, IS_EX, IS_HA, IS_LI> >
{};


template<BOOST_PP_ENUM_BINARY_PARAMS(FUSION_MAX_VECTOR_SIZE,typename A,=DefaultArgument BOOST_PP_INTERCEPT)> 
struct Make : boost::mpl::identity<typename Base<typename make_list<BOOST_PP_ENUM_PARAMS(FUSION_MAX_VECTOR_SIZE,A)>::type>::Ptr>
{};



} // result_of

} // composite


#define DEFAULT_print(z, n, data) DefaultArgument

#define BOOST_PP_ITERATION_LIMITS (1,BOOST_PP_SUB(FUSION_MAX_VECTOR_SIZE,1) )
#define BOOST_PP_FILENAME_1 "Composite.h"

#include BOOST_PP_ITERATE()

#undef BOOST_PP_FILENAME_1
#undef BOOST_PP_ITERATION_LIMITS

#undef DEFAULT_print

#endif // COMPOSITES_COMPOSITE_H_INCLUDED


#else  // BOOST_PP_IS_ITERATING

#define ITER BOOST_PP_ITERATION()

namespace composite {


namespace result_of {


template<bool IS_EX, bool IS_HA, bool IS_LI, BOOST_PP_ENUM_PARAMS(ITER,typename A)>
struct MakeConcrete<IS_EX, IS_HA, IS_LI, BOOST_PP_ENUM_PARAMS(ITER,A) BOOST_PP_ENUM_TRAILING(BOOST_PP_SUB(FUSION_MAX_VECTOR_SIZE,ITER),DEFAULT_print,~) >
  : boost::mpl::identity<Composite<typename make_list<BOOST_PP_ENUM_PARAMS(ITER,A) >::type, IS_EX, IS_HA, IS_LI> >
{};


template<BOOST_PP_ENUM_PARAMS(ITER,typename A)>
struct Make<BOOST_PP_ENUM_PARAMS(ITER,A) BOOST_PP_ENUM_TRAILING(BOOST_PP_SUB(FUSION_MAX_VECTOR_SIZE,ITER),DEFAULT_print,~) >
  : boost::mpl::identity<typename Base<typename make_list<BOOST_PP_ENUM_PARAMS(ITER,A) >::type>::Ptr>
{};


} // result_of


#define RETURN_type typename result_of::Make<BOOST_PP_ENUM_PARAMS(ITER,A) >::type

template<BOOST_PP_ENUM_PARAMS(ITER,typename A)> 
const RETURN_type
make(BOOST_PP_ENUM_BINARY_PARAMS(ITER,const A,& act) )
{
  return doMake(make_list(BOOST_PP_ENUM_PARAMS(ITER,act)));
}

#undef  RETURN_type


} // composite


#undef  ITER


#endif // BOOST_PP_IS_ITERATING
