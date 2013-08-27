// -*- C++ -*-
/// \briefFile{Defines the hierarchical partial specializations of structure::ElementLiouvillean}
#ifndef STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED
#define STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED

#include "ElementLiouvilleanFwd.h"

#include "Liouvillean.h"
#include "ElementLiouvilleanAveragedCommon.h"

#include <boost/function.hpp>

#include <boost/mpl/inherit_linearly.hpp>
#include <boost/mpl/inherit.hpp>


namespace structure {


namespace details {


template<int RANK, typename JUMP_ORDINAL_ICW, bool IS_TIME_DEPENDENT>
class Lindblad;

/*
template<int RANK, typename JUMP_ORDINAL_ICW>
class Lindblad<RANK,JUMP_ORDINAL_ICW,true>
{
public:
  static const int JUMP_ORDINAL=JUMP_ORDINAL_ICW::value;
  
protected:
  // Type-erasure & non-virtual interface idiom in one:
  
  void typeErasedActWithJ(double t, typename quantumdata::Types<RANK>::StateVectorLow& psi) const {doActWithJ(t,psi,mpl::int_<JUMP_ORDINAL>());}
  
  double typeErasedRate(double t, const quantumdata::LazyDensityOperator<RANK>& matrix) const {return rate(t,matrix,mpl::int_<JUMP_ORDINAL>());}
  
private:
  virtual void doActWithJ(double t, typename quantumdata::Types<RANK>::StateVectorLow&, mpl::int_<JUMP_ORDINAL>) const = 0;
  
  virtual double rate(double t, const quantumdata::LazyDensityOperator<RANK>&, mpl::int_<JUMP_ORDINAL>) const = 0;
  
};
*/

template<int RANK, typename JUMP_ORDINAL_ICW>
class Lindblad<RANK,JUMP_ORDINAL_ICW,false>
{
public:
  static const int JUMP_ORDINAL=JUMP_ORDINAL_ICW::value;
  
protected:
  void typeErasedActWithJ(typename quantumdata::Types<RANK>::StateVectorLow& psi) const {doActWithJ(psi,mpl::int_<JUMP_ORDINAL>());}
  
  double typeErasedRate(const quantumdata::LazyDensityOperator<RANK>& matrix) const {return rate(matrix,mpl::int_<JUMP_ORDINAL>());}
  
private:
  virtual void doActWithJ(typename quantumdata::Types<RANK>::StateVectorLow&, mpl::int_<JUMP_ORDINAL>) const = 0;
  
  virtual double rate(const quantumdata::LazyDensityOperator<RANK>&, mpl::int_<JUMP_ORDINAL>) const = 0;
  
};

} // details


template<int RANK, int NOJ>
class ElementLiouvillean<RANK,NOJ,false> : public ElementLiouvilleanAveragedCommon<Liouvillean<RANK,false> >,
                                           private mpl::inherit_linearly<tmptools::Ordinals<NOJ>,mpl::inherit<mpl::_2,details::Lindblad<RANK,mpl::_1,false> > >::type
{
private:
  typedef ElementLiouvilleanAveragedCommon<Liouvillean<RANK,false> > Base;
  
public:
  typedef typename Base::StateVectorLow StateVectorLow;

  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  typedef typename Base::Rates Rates;
  
  template<int JUMP_ORDINAL>
  class JumpNo : mpl::int_<JUMP_ORDINAL> {BOOST_STATIC_ASSERT(JUMP_ORDINAL<NOJ);}; // or some enable_if-based solution
  
  typedef boost::function<void  (      StateVectorLow&     )> JumpStrategy    ;
  typedef boost::function<double(const LazyDensityOperator&)> JumpRateStrategy;

  typedef blitz::TinyVector<JumpStrategy    ,NOJ> JumpStrategies    ;
  typedef blitz::TinyVector<JumpRateStrategy,NOJ> JumpRateStrategies;

protected:
  template<typename... KeyLabelsPack>
  ElementLiouvillean(const std::string& keyTitle, KeyLabelsPack&&... keyLabelsPack) : Base(keyTitle,keyLabelsPack...), jumps_(fillJumpStrategies()), jumpRates_(fillJumpRateStrategies()) {}
  
private:
  const Rates average_v(const LazyDensityOperator&) const;

  void actWithJ_v(StateVectorLow& psi, size_t jumpNo) const;
  
  const JumpStrategies     fillJumpStrategies    () const;
  const JumpRateStrategies fillJumpRateStrategies() const;
  
  const JumpStrategies     jumps_    ;
  const JumpRateStrategies jumpRates_;

};



template<int RANK>
class ElementLiouvillean<RANK,1,false> : public ElementLiouvilleanAveragedCommon<Liouvillean<RANK,false> >
{
private:
  typedef ElementLiouvilleanAveragedCommon<Liouvillean<RANK,false> > Base;
  
public:
  typedef typename Base::StateVectorLow StateVectorLow;

  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  typedef typename Base::Rates Rates;
  
  struct ElementLiouvilleanException : cpputils::TaggedException
  {
    ElementLiouvilleanException(const std::string& tag) : cpputils::TaggedException(tag) {}
  };
  
protected:
  ElementLiouvillean(const std::string& keyTitle, const std::string& keyLabel) : Base(keyTitle,1,keyLabel) {}
  
private:
  const Rates average_v(const LazyDensityOperator& matrix) const {Rates rates(1); rates(0)=rate(matrix); return rates;}

  void actWithJ_v(StateVectorLow& psi, size_t jumpNo) const {if (jumpNo) throw ElementLiouvilleanException(Base::getTitle()); doActWithJ(psi);}

  virtual void doActWithJ(StateVectorLow&) const = 0;
  
  virtual double rate(const LazyDensityOperator&) const = 0;

};


/*
template<int RANK, int NOJ>
const LiouvilleanAveragedCommon::DArray1D ElementLiouvillean<RANK,NOJ,ISTD_true_false>::average_v(COND_ARG_T(ISTD) const LazyDensityOperator& matrix) const
{
  Rates rates(NOJ);
  // Note that this cannot be anything like static because of the by-reference semantics of blitz::Array

  boost::transform(jumpRates_,rates.begin(),
                   bind(&JumpRateStrategy::operator(),_1,BOOST_PP_EXPR_IIF(ISTD,t) BOOST_PP_COMMA_IF(ISTD) boost::cref(matrix)));
  return rates;
}
*/

} // structure


#endif // STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED


