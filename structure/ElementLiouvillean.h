// -*- C++ -*-
/// \briefFile{Defines the hierarchical partial specializations of structure::Liouvillean}
#ifndef STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED
#define STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED

#include "ElementLiouvilleanFwd.h"

#include "Liouvillean.h"
#include "ElementLiouvilleanAveragedCommon.h"

#include "Range.h"

#include <boost/function.hpp>
#include <boost/bind.hpp>


namespace structure {


template<int RANK, int NOJ, bool IS_TIME_DEPENDENT>
class ElementLiouvillean : public ElementLiouvilleanAveragedCommon<Liouvillean<RANK,IS_TIME_DEPENDENT> >
{
private:
  typedef Liouvillean<RANK,IS_TIME_DEPENDENT> Base;
  
public:
  typedef typename Base::StateVectorLow StateVectorLow;

  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  typedef typename Base::Rates Rates;

protected:
  template<typename... KeyLabelsPack>
  ElementLiouvillean(const JumpStrategies& jumps, const JumpRateStrategies& jumpRates, const std::string& keyTitle, KeyLabelsPack&&... keyLabelsPack) 
    : jumps_(jumps), jumpRates_(jumpRates), keyPrinter_(keyTitle,keyLabelsPack) {}

private:
  size_t nAvr_v() const {return NOJ;}

  const Rates average_v(double t, const LazyDensityOperator&) const;
  const Rates average_v(          const LazyDensityOperator&) const;

  void actWithJ_v(double t, StateVectorLow& psi, size_t jumpNo) const;
  void actWithJ_v(          StateVectorLow& psi, size_t jumpNo) const;

};


template<int RANK, bool IS_TIME_DEPENDENT>
class ElementLiouvillean<RANK,1,IS_TIME_DEPENDENT> : public ElementLiouvilleanAveragedCommon<Liouvillean<RANK,IS_TIME_DEPENDENT> >
{
private:
  typedef ElementLiouvilleanAveragedCommon<Liouvillean<RANK,IS_TIME_DEPENDENT> >
  
public:
  typedef typename Base::StateVectorLow StateVectorLow;

  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  typedef typename Base::Rates Rates;
  
protected:
  ElementLiouvillean(const std::string& keyTitle, const std::string& keyLabel) : Base(keyTitle,1,keyLabel) {}
  
private:
  const Rates average_v(double t, const LazyDensityOperator& matrix) const {Rates rates(1); rates(0)=rate(t,matrix); return rates;}

  void actWithJ_v(double t, StateVectorLow& psi, size_t jumpNo) const {if (jumpNo) throw typename Base::ElementLiouvilleanException(); doActWithJ(t,psi);}

  virtual void doActWithJ(double t, StateVectorLow&) const = 0;
  
  virtual double rate(double t, const LazyDensityOperator&) const = 0;

};


template<int RANK>
class ElementLiouvillean<RANK,1,false> : public ElementLiouvilleanOneJump<RANK,false>
{
private:
  typedef ElementLiouvilleanOneJump<RANK,false> Base;
  
public:
  typedef typename Base::StateVectorLow StateVectorLow;

  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  typedef typename Base::Rates Rates;
  
protected:
  ElementLiouvillean(const std::string& keyTitle, const std::string& keyLabel) : Base(keyTitle,keyLabel) {}
  
private:
  const Rates average_v(const LazyDensityOperator& matrix) const {Rates rates(1); rates(0)=rate(matrix); return rates;}

  void actWithJ_v(StateVectorLow& psi, size_t jumpNo) const {if (jumpNo) throw typename Base::ElementLiouvilleanException(); doActWithJ(psi);}

  virtual void doActWithJ(StateVectorLow&) const = 0;
  
  virtual double rate(const LazyDensityOperator&) const = 0;

};

template<int RANK, int NOJ, bool IS_TIME_DEPENDENT>
void ElementLiouvillean<RANK,NOJ,IS_TIME_DEPENDENT>::actWithJ_v(double t, StateVectorLow& psi, size_t jumpNo) const {typename boost::enable_if_c< IS_TIME_DEPENDENT>::type();/*jumps_(jumpNo)(t,psi);*/}

template<int RANK, int NOJ, bool IS_TIME_DEPENDENT>
void ElementLiouvillean<RANK,NOJ,IS_TIME_DEPENDENT>::actWithJ_v(          StateVectorLow& psi, size_t jumpNo) const {typename boost::enable_if_c<!IS_TIME_DEPENDENT>::type();/*jumps_(jumpNo)(  psi);*/}

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


