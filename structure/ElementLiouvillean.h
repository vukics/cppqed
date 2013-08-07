// -*- C++ -*-
#ifndef ELEMENTLIOUVILLEAN_ISTD_REENTERING

#ifndef STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED
#define STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED

#include "ElementLiouvilleanFwd.h"

#include "Liouvillean.h"

#include "KeyPrinter.h"
#include "Range.h"

#ifndef   NDEBUG
#include "Exception.h"
#endif // NDEBUG

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <boost/preprocessor/control/expr_iif.hpp>
#include <boost/preprocessor/control/iif.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>

namespace structure {


#define ELEMENTLIOUVILLEAN_ISTD_REENTERING 0
#include "ElementLiouvillean.h"

#define ELEMENTLIOUVILLEAN_ISTD_REENTERING 1
#include "ElementLiouvillean.h"


} // structure



#endif // STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED


#else  // ELEMENTLIOUVILLEAN_ISTD_REENTERING


#define ISTD ELEMENTLIOUVILLEAN_ISTD_REENTERING

#define COND_ARG(IS_TIME_DEPENDENT) BOOST_PP_EXPR_IIF(ISTD,double) BOOST_PP_COMMA_IF(ISTD)
#define COND_ARG_T(IS_TIME_DEPENDENT) BOOST_PP_EXPR_IIF(ISTD,double t) BOOST_PP_COMMA_IF(ISTD)

#define ISTD_true_false BOOST_PP_IIF(ISTD,true,false)


template<int RANK, int NOJ>
class ElementLiouvillean<RANK,NOJ,ISTD_true_false> : public Liouvillean<RANK,ISTD_true_false>
{
public:
  typedef Liouvillean<RANK,ISTD_true_false> Base;

  typedef typename Base::     StateVectorLow      StateVectorLow;

  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  typedef typename LiouvilleanAveragedCommon::DArray1D Rates;

  typedef boost::function<void  (COND_ARG(ISTD)       StateVectorLow&     )> JumpStrategy;
  typedef boost::function<double(COND_ARG(ISTD) const LazyDensityOperator&)> JumpRateStrategy;

  typedef blitz::TinyVector<JumpStrategy           ,NOJ> JumpStrategies;
  typedef blitz::TinyVector<JumpRateStrategy,NOJ> JumpRateStrategies;

  typedef cpputils::KeyPrinter::KeyLabels KeyLabels;

protected:
  ElementLiouvillean(const JumpStrategies& jumps, const JumpRateStrategies& jumpRates, const std::string& keyTitle, const KeyLabels& keyLabels) 
    : jumps_(jumps), jumpRates_(jumpRates), keyPrinter_(keyTitle,keyLabels) {}

private:
  size_t nAvr_v() const {return NOJ;}

  const Rates average_v(COND_ARG(ISTD) const LazyDensityOperator&) const;

  void actWithJ_v(COND_ARG_T(ISTD) StateVectorLow& psi, size_t jumpNo) const {jumps_(jumpNo)(BOOST_PP_EXPR_IIF(ISTD,t) BOOST_PP_COMMA_IF(ISTD) psi);}

  std::ostream& displayKey_v(std::ostream& os, size_t& i) const {return keyPrinter_.displayKey(os,i);}

  const JumpStrategies     jumps_    ;
  const JumpRateStrategies jumpRates_;

  const cpputils::KeyPrinter keyPrinter_;

};


template<int RANK>
class ElementLiouvillean<RANK,1,ISTD_true_false> : public Liouvillean<RANK,ISTD_true_false>
// This specialization can use the virtual-function technique of old
{
public:
  typedef Liouvillean<RANK,ISTD_true_false> Base;

  typedef typename Base::    StateVectorLow     StateVectorLow;

  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  typedef typename LiouvilleanAveragedCommon::DArray1D Rates;

  typedef cpputils::KeyPrinter::KeyLabels KeyLabels;

protected:
  ElementLiouvillean(const std::string& keyTitle, const std::string& keyLabel) 
    : keyPrinter_(keyTitle,KeyLabels(1,keyLabel)) {}

private:
  size_t nAvr_v() const {return 1;}

  const Rates average_v(COND_ARG_T(ISTD) const LazyDensityOperator& matrix) const {Rates rates(1); rates(0)=rate(BOOST_PP_EXPR_IIF(ISTD,t) BOOST_PP_COMMA_IF(ISTD) matrix); return rates;}

#ifndef   NDEBUG
  struct ElementLiouvilleanException : cpputils::Exception {};
#endif // NDEBUG

  void actWithJ_v(COND_ARG_T(ISTD) StateVectorLow& psi, size_t 
#ifndef   NDEBUG
                  jumpNo
#endif // NDEBUG
                  ) const {
#ifndef   NDEBUG
    if (jumpNo) throw ElementLiouvilleanException(); 
#endif // NDEBUG
    doActWithJ(BOOST_PP_EXPR_IIF(ISTD,t) BOOST_PP_COMMA_IF(ISTD) psi);
  }

  std::ostream& displayKey_v(std::ostream& os, size_t& i) const {return keyPrinter_.displayKey(os,i);}

  virtual void   doActWithJ(COND_ARG(ISTD)       StateVectorLow     &) const = 0;
  virtual double rate      (COND_ARG(ISTD) const LazyDensityOperator&) const = 0;

  const cpputils::KeyPrinter keyPrinter_;

};




template<int RANK, int NOJ>
const LiouvilleanAveragedCommon::DArray1D ElementLiouvillean<RANK,NOJ,ISTD_true_false>::average_v(COND_ARG_T(ISTD) const LazyDensityOperator& matrix) const
{
  Rates rates(NOJ);
  // Note that this cannot be anything like static because of the by-reference semantics of blitz::Array

  boost::transform(jumpRates_,rates.begin(),
                   bind(&JumpRateStrategy::operator(),_1,BOOST_PP_EXPR_IIF(ISTD,t) BOOST_PP_COMMA_IF(ISTD) boost::cref(matrix)));
  return rates;
}


#undef ISTD_true_false

#undef COND_ARG
#undef COND_ARG_T


#undef ISTD
#undef ELEMENTLIOUVILLEAN_ISTD_REENTERING

#endif  // ELEMENTLIOUVILLEAN_ISTD_REENTERING