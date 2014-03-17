// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
// -*- C++ -*-
/// \briefFile{Defines the structure::ElementLiouvillean & structure::ElementLiouvilleanStrategies}
#ifndef CPPQEDCORE_STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED
#define CPPQEDCORE_STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED

#include "ElementLiouvilleanFwd.h"

#include "Liouvillean.h"
#include "ElementLiouvilleanAveragedCommon.h"

#include <boost/function.hpp>

#include <boost/mpl/for_each.hpp>


namespace structure {


/// Tools for describing a single Lindblad operator
namespace lindblad {


/// A class to encapsulate Base::LindbladNo, so that the latter has information about the number of Lindblads
/** \tparam NLINDBLADS the number of Lindblads */
template<int NLINDBLADS>
class Base
{
protected:
  /// A tagging class for Lindblad
  template<int LINDBLAD_ORDINAL, typename OTHER=typename boost::enable_if_c< (LINDBLAD_ORDINAL<NLINDBLADS) >::type>
  class LindbladNo : mpl::int_<LINDBLAD_ORDINAL> {};

};


/// Class defining the virtual functions corresponding to a single Lindblad…
/**
 * … and all the Lindblads with lesser index, since the class uses base-class chaining. At the bottom of the chain sits Base`<NLINDBLADS>`
 * 
 * \tparamRANK
 * \tparam LINDBLAD_ORDINAL the index of the given Lindblad
 * \tparam IS_TIME_DEPENDENT governs time dependence
 * \tparam NLINDBLADS the total number of Lindblads
 * 
 * In the case of the _::typeErasedActWithJ and _::typeErasedRate functions, the simultaneous use of the type-erasure & non-virtual interface idioms is notable
 * 
 */
template<int RANK, int LINDBLAD_ORDINAL, bool IS_TIME_DEPENDENT, int NLINDBLADS=LINDBLAD_ORDINAL+1>
class _ : public mpl::if_c<LINDBLAD_ORDINAL,_<RANK,LINDBLAD_ORDINAL-1,IS_TIME_DEPENDENT,NLINDBLADS>,Base<NLINDBLADS> >::type
{
public:
  typedef typename time::DispatcherIsTimeDependent<IS_TIME_DEPENDENT>::type Time;
  
  /// calls the virtual doActWithJ for the given LINDBLAD_ORDINAL
  void typeErasedActWithJ(Time t, typename quantumdata::Types<RANK>::StateVectorLow& psi) const {doActWithJ(t,psi,typename Base<NLINDBLADS>::template LindbladNo<LINDBLAD_ORDINAL>());}
  
  /// calls the virtual rate for the given LINDBLAD_ORDINAL
  double typeErasedRate(Time t, const quantumdata::LazyDensityOperator<RANK>& matrix) const {return rate(t,matrix,typename Base<NLINDBLADS>::template LindbladNo<LINDBLAD_ORDINAL>());}
  
private:
  virtual void doActWithJ(Time, typename quantumdata::Types<RANK>::StateVectorLow&, typename Base<NLINDBLADS>::template LindbladNo<LINDBLAD_ORDINAL>) const = 0;
  
  virtual double rate(Time, const quantumdata::LazyDensityOperator<RANK>&, typename Base<NLINDBLADS>::template LindbladNo<LINDBLAD_ORDINAL>) const = 0;
  
};


} // lindblad


/// Thrown if the Lindblad index is not smaller than the total number of Lindblads
struct ElementLiouvilleanException : cpputils::TaggedException
{
  ElementLiouvilleanException(const std::string& tag) : cpputils::TaggedException(tag) {}
};


/// An implementation of Liouvillean for the case when the number of Lindblads is known @ compile time (which is very often the case with elements)…
/**
 * … if this is not the case, Liouvillean has to be used.
 * 
 * The class relies on lindblad::_ to declare a pair of purely virtual functions `doActWithJ` and `rate` for each Lindblad between 0 and NLINDBLADS-1.
 * The virtuals `average_v` and `actWithJ_v` inherited from Liouvillean are defined in terms of these purely virtual functions.
 * 
 * \tparamRANK
 * \tparam NLINDBLADS the total number of Lindblads
 * \tparam IS_TIME_DEPENDENT governs time dependence
 * 
 * \see Sec. \ref hierarchicaloscillator of the structure-bundle guide for an example of usage 
 * 
 */
template<int RANK, int NLINDBLADS, bool IS_TIME_DEPENDENT>
class ElementLiouvillean : public ElementLiouvilleanAveragedCommon<LiouvilleanTimeDependenceDispatched<RANK,IS_TIME_DEPENDENT> >,
                           public lindblad::_<RANK,NLINDBLADS-1,false>
{
private:
  typedef ElementLiouvilleanAveragedCommon<LiouvilleanTimeDependenceDispatched<RANK,IS_TIME_DEPENDENT> > Base;
  
public:
  typedef typename Base::StateVectorLow StateVectorLow;

  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  typedef typename Base::Rates Rates;
  
  typedef typename Base::Time Time;
  
protected:
  template<typename... KeyLabelsPack>
  ElementLiouvillean(const std::string& keyTitle, KeyLabelsPack&&... keyLabelsPack) : Base(keyTitle,keyLabelsPack...) {}
  
  ElementLiouvillean(const std::string& keyTitle, typename Base::KeyLabelsInitializer il) : Base(keyTitle,il) {}

private:
  class Average
  {
  public:
    Average(const ElementLiouvillean* ptr, Rates& rates, Time t, const LazyDensityOperator& matrix) : ptr_(ptr), rates_(rates), t_(t), matrix_(matrix) {}

    template<typename T> void operator()(T) const {rates_(T::value)=static_cast<const lindblad::_<RANK,T::value,false,NLINDBLADS>*const>(ptr_)->typeErasedRate(t_,matrix_);}

  private:
    const ElementLiouvillean*const ptr_;
    Rates& rates_;
    const Time t_;
    const LazyDensityOperator& matrix_;    
  };

  virtual const Rates rates_v(Time t, const LazyDensityOperator& matrix) const final {Rates rates(NLINDBLADS); mpl::for_each<tmptools::Ordinals<NLINDBLADS> >(Average(this,rates,t,matrix)); return rates;}

  class ActWithJ
  {
  public:
    ActWithJ(const ElementLiouvillean* ptr, Time t, StateVectorLow& psi, size_t lindbladNo) : ptr_(ptr), t_(t), psi_(psi), lindbladNo_(lindbladNo) {}

    template<typename T>
    void operator()(T) const {if (T::value==lindbladNo_) static_cast<const lindblad::_<RANK,T::value,false,NLINDBLADS>*const>(ptr_)->typeErasedActWithJ(t_,psi_);}

  private:
    const ElementLiouvillean*const ptr_;
    const Time t_;
    StateVectorLow& psi_;
    const size_t lindbladNo_;
  };

  virtual void actWithJ_v(Time t, StateVectorLow& psi, size_t lindbladNo) const final {if (lindbladNo>=NLINDBLADS) throw ElementLiouvilleanException(Base::getTitle()); mpl::for_each<tmptools::Ordinals<NLINDBLADS> >(ActWithJ(this,t,psi,lindbladNo));}
  
};



/// A specialization of ElementLiouvillean for the case of a single Lindblad 
/**
 * Here, the tagging class lindblad::_::LindbladNo is not needed for the signature of the pure virtuals in the implementer’s interface.
 * 
 * \tparamRANK
 * \tparam IS_TIME_DEPENDENT governs time dependence
 * 
 */
template<int RANK, bool IS_TIME_DEPENDENT>
class ElementLiouvillean<RANK,1,IS_TIME_DEPENDENT> : public ElementLiouvilleanAveragedCommon<LiouvilleanTimeDependenceDispatched<RANK,IS_TIME_DEPENDENT> >
{
private:
  typedef ElementLiouvilleanAveragedCommon<LiouvilleanTimeDependenceDispatched<RANK,IS_TIME_DEPENDENT> > Base;
  
public:
  typedef typename Base::StateVectorLow StateVectorLow;

  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  typedef typename Base::Rates Rates;

  typedef typename Base::Time Time;
  
protected:
  ElementLiouvillean(const std::string& keyTitle, const std::string& keyLabel) : Base(keyTitle,1,keyLabel) {}
  
private:
  virtual const Rates rates_v(Time t, const LazyDensityOperator& matrix) const final {Rates rates(1); rates(0)=rate(t,matrix); return rates;}

  virtual void actWithJ_v(Time t, StateVectorLow& psi, size_t lindbladNo) const final {if (lindbladNo) throw ElementLiouvilleanException(Base::getTitle()); doActWithJ(t,psi);}

  virtual void doActWithJ(Time, StateVectorLow&) const = 0;
  
  virtual double rate(Time, const LazyDensityOperator&) const = 0;

};


/// Besides ElementLiouvillean, this is another solution based on the strategy idiom to control the number of Lindblads @ compile time
/**
 * \tparamRANK
 * \tparam NLINDBLADS the total number of Lindblads
 * \tparam IS_TIME_DEPENDENT governs time dependence
 * 
 * \see Sec. \ref basicoscillator of the structure-bundle guide for an example of usage 
 */
template<int RANK, int NLINDBLADS, bool IS_TIME_DEPENDENT>
class ElementLiouvilleanStrategies : public ElementLiouvilleanAveragedCommon<LiouvilleanTimeDependenceDispatched<RANK,IS_TIME_DEPENDENT> >
{
private:
  typedef ElementLiouvilleanAveragedCommon<LiouvilleanTimeDependenceDispatched<RANK,IS_TIME_DEPENDENT> > Base;
  
public:
  typedef typename Base::StateVectorLow StateVectorLow;

  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  typedef typename Base::Rates Rates;
  
  typedef typename Base::Time Time;

  /// Strategy functional for acting with a given Lindblad operator (= performing a given jump) on a state
  /** The actual signature of the functional is decided on the basis of IS_TIME_DEPENDENT using the compile-time *if*-construct \refBoostConstruct{if_c,mpl/doc/refmanual/if-c.html}. */
  typedef typename mpl::if_c<IS_TIME_DEPENDENT,boost::function<void  (double,       StateVectorLow&     )>,boost::function<void  (      StateVectorLow&     )> >::type JumpStrategy;
  /// Strategy functional for calculating from a state the jump rate corresponding to a given Lindblad
  typedef typename mpl::if_c<IS_TIME_DEPENDENT,boost::function<double(double, const LazyDensityOperator&)>,boost::function<double(const LazyDensityOperator&)> >::type JumpRateStrategy;

  /// Tiny vector of length NLINDBLADS containing the #JumpStrategy instances
  typedef blitz::TinyVector<JumpStrategy    ,NLINDBLADS> JumpStrategies;
  /// Tiny vector of length NLINDBLADS containing the #JumpRateStrategy instances
  typedef blitz::TinyVector<JumpRateStrategy,NLINDBLADS> JumpRateStrategies;

protected:
  template<typename... KeyLabelsPack>
  ElementLiouvilleanStrategies(const JumpStrategies& jumps, const JumpRateStrategies& jumpRates, const std::string& keyTitle, KeyLabelsPack&&... keyLabelsPack)
    : Base(keyTitle,keyLabelsPack...), jumps_(jumps), jumpRates_(jumpRates) {}
  
  ElementLiouvilleanStrategies(const JumpStrategies& jumps, const JumpRateStrategies& jumpRates, const std::string& keyTitle, typename Base::KeyLabelsInitializer il)
    : Base(keyTitle,il), jumps_(jumps), jumpRates_(jumpRates) {}

private:
  virtual const Rates rates_v(Time t, const LazyDensityOperator& matrix) const final;

  virtual void actWithJ_v(Time t, StateVectorLow& psi, size_t lindbladNo) const final;
  
  const JumpStrategies     jumps_    ;
  const JumpRateStrategies jumpRates_;

};


} // structure


#endif // CPPQEDCORE_STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED


