// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines the structure::ElementLiouvillean & structure::ElementLiouvilleanStrategies}
#ifndef CPPQEDCORE_STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED
#define CPPQEDCORE_STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED

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
class _ : public mpl::if_c< (LINDBLAD_ORDINAL>0) ,_<RANK,LINDBLAD_ORDINAL-1,IS_TIME_DEPENDENT,NLINDBLADS>,Base<NLINDBLADS> >::type
{
public:
  typedef typename time::DispatcherIsTimeDependent<IS_TIME_DEPENDENT>::type Time;
  
  /// calls the virtual doActWithJ for the given LINDBLAD_ORDINAL
  void typeErasedActWithJ(Time t, typename quantumdata::Types<RANK>::StateVectorLow& psi) const {doActWithJ(t,psi,typename Base<NLINDBLADS>::template LindbladNo<LINDBLAD_ORDINAL>());}
  
  /// calls the virtual rate for the given LINDBLAD_ORDINAL
  double typeErasedRate(Time t, const quantumdata::LazyDensityOperator<RANK>& matrix) const {return rate(t,matrix,typename Base<NLINDBLADS>::template LindbladNo<LINDBLAD_ORDINAL>());}

  void typeErasedActWithSuperoperator(Time t,
                                      const typename quantumdata::Types<RANK>::DensityOperatorLow& rho,
                                      typename quantumdata::Types<RANK>::DensityOperatorLow& drhodt
                                     ) const {return doActWithSuperoperator(t,rho,drhodt,typename Base<NLINDBLADS>::template LindbladNo<LINDBLAD_ORDINAL>());}

private:
  virtual void doActWithJ(Time, typename quantumdata::Types<RANK>::StateVectorLow&, typename Base<NLINDBLADS>::template LindbladNo<LINDBLAD_ORDINAL>) const = 0;
  
  virtual double rate(Time, const quantumdata::LazyDensityOperator<RANK>&, typename Base<NLINDBLADS>::template LindbladNo<LINDBLAD_ORDINAL>) const = 0;

  virtual void doActWithSuperoperator(Time,
                                      const typename quantumdata::Types<RANK>::DensityOperatorLow&,
                                      typename quantumdata::Types<RANK>::DensityOperatorLow&,
                                      typename Base<NLINDBLADS>::template LindbladNo<LINDBLAD_ORDINAL>
                                    ) const {throw SuperoperatorNotImplementedException(LINDBLAD_ORDINAL);}
  
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
template<int RANK, int NLINDBLADS=1, bool IS_TIME_DEPENDENT=false>
class ElementLiouvillean : public ElementLiouvilleanAveragedCommon<LiouvilleanTimeDependenceDispatched<RANK,IS_TIME_DEPENDENT> >,
                           public lindblad::_<RANK,NLINDBLADS-1,false>
{
private:
  typedef ElementLiouvilleanAveragedCommon<LiouvilleanTimeDependenceDispatched<RANK,IS_TIME_DEPENDENT> > Base;
  
public:
#define TYPE_DEFINITION_FORWARD typedef typename Base::StateVectorLow StateVectorLow; typedef typename Base::DensityOperatorLow DensityOperatorLow; typedef typename Base::LazyDensityOperator LazyDensityOperator; typedef typename Base::Rates Rates; typedef typename Base::Time Time;  
  TYPE_DEFINITION_FORWARD
  
protected:
  template<typename... KeyLabelsPack>
  ElementLiouvillean(const std::string& keyTitle, KeyLabelsPack&&... keyLabelsPack) : Base(keyTitle,keyLabelsPack...) {}
  
  ElementLiouvillean(const std::string& keyTitle, typename Base::KeyLabelsInitializer il) : Base(keyTitle,il) {}

private:
  struct Average
  {
    template<typename T>
    void operator()(T) const {rates_(T::value)=static_cast<const lindblad::_<RANK,T::value,false,NLINDBLADS>*const>(ptr_)->typeErasedRate(t_,matrix_);}

    const ElementLiouvillean*const ptr_;
    Rates& rates_;
    const Time t_;
    const LazyDensityOperator& matrix_;    
  };

  const Rates rates_v(Time t, const LazyDensityOperator& matrix) const final
  {
    Rates rates(NLINDBLADS);
    mpl::for_each<tmptools::Ordinals<NLINDBLADS> >(Average{this,rates,t,matrix});
    return rates;
  }

  struct ActWithJ
  {
    template<typename T>
    void operator()(T) const {if (T::value==lindbladNo_) static_cast<const lindblad::_<RANK,T::value,false,NLINDBLADS>*const>(ptr_)->typeErasedActWithJ(t_,psi_);}

    const ElementLiouvillean*const ptr_;
    const Time t_;
    StateVectorLow& psi_;
    const size_t lindbladNo_;
  };

  void actWithJ_v(Time t, StateVectorLow& psi, size_t lindbladNo) const final
  {
    if (lindbladNo>=NLINDBLADS) throw ElementLiouvilleanException(Base::getTitle());
    mpl::for_each<tmptools::Ordinals<NLINDBLADS> >(ActWithJ{this,t,psi,lindbladNo});
  }
  
  struct ActWithSuperoperator
  {
    template<typename T>
    void operator()(T) const {if (T::value==lindbladNo_) static_cast<const lindblad::_<RANK,T::value,false,NLINDBLADS>*const>(ptr_)->typeErasedActWithSuperoperator(t_,rho_,drhodt_);}

    const ElementLiouvillean*const ptr_;
    const Time t_;
    const DensityOperatorLow& rho_;
    DensityOperatorLow& drhodt_;
    const size_t lindbladNo_;
  };

  void actWithSuperoperator_v(Time t, const DensityOperatorLow& rho, DensityOperatorLow& drhodt, size_t lindbladNo) const final
  {
    if (lindbladNo>=NLINDBLADS) throw ElementLiouvilleanException(Base::getTitle());
    mpl::for_each<tmptools::Ordinals<NLINDBLADS> >(ActWithSuperoperator{this,t,rho,drhodt,lindbladNo});    
  }

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
  TYPE_DEFINITION_FORWARD
  
protected:
  ElementLiouvillean(const std::string& keyTitle, const std::string& keyLabel) : Base(keyTitle,1,keyLabel) {}
  
private:
  const Rates rates_v(Time t, const LazyDensityOperator& matrix) const final {Rates rates(1); rates(0)=rate(t,matrix); return rates;}

  void actWithJ_v(Time t, StateVectorLow& psi, size_t lindbladNo) const final
  {
    if (lindbladNo) throw ElementLiouvilleanException(Base::getTitle());
    doActWithJ(t,psi);
  }

  void actWithSuperoperator_v(Time t, const DensityOperatorLow& rho, DensityOperatorLow& drhodt, size_t lindbladNo) const final
  {
    if (lindbladNo) throw ElementLiouvilleanException(Base::getTitle());
    doActWithSuperoperator(t,rho,drhodt);
  }
  
  virtual void doActWithJ(Time, StateVectorLow&) const = 0;
  
  virtual double rate(Time, const LazyDensityOperator&) const = 0;

  virtual void doActWithSuperoperator(Time, const DensityOperatorLow&, DensityOperatorLow&) const {throw SuperoperatorNotImplementedException(0);}

};


/// Simply a less templated base to ElementLiouvilleanStrategies defininig the strategy funcional types
template<int RANK, bool IS_TIME_DEPENDENT>
class ElementLiouvilleanStrategiesBase : public ElementLiouvilleanAveragedCommon<LiouvilleanTimeDependenceDispatched<RANK,IS_TIME_DEPENDENT> >
{
private:
  typedef ElementLiouvilleanAveragedCommon<LiouvilleanTimeDependenceDispatched<RANK,IS_TIME_DEPENDENT> > Base;

protected:
  using Base::Base; // inherit constructor
  
public:
  TYPE_DEFINITION_FORWARD

  /// Strategy functional for acting with a given Lindblad operator (= performing a given jump) on a state
  /** The actual signature of the functional is decided on the basis of IS_TIME_DEPENDENT using the compile-time *if*-construct \refBoostConstruct{if_c,mpl/doc/refmanual/if-c.html}. */
  typedef typename mpl::if_c<IS_TIME_DEPENDENT,boost::function<void  (double,       StateVectorLow&     )>,boost::function<void  (      StateVectorLow&     )> >::type JumpStrategy;
  /// Strategy functional for calculating from a state the jump rate corresponding to a given Lindblad
  typedef typename mpl::if_c<IS_TIME_DEPENDENT,boost::function<double(double, const LazyDensityOperator&)>,boost::function<double(const LazyDensityOperator&)> >::type JumpRateStrategy;
  
  typedef typename mpl::if_c<IS_TIME_DEPENDENT,
                            boost::function<void(double, const DensityOperatorLow&, DensityOperatorLow&)>,
                            boost::function<void(        const DensityOperatorLow&, DensityOperatorLow&)> >::type SuperoperatorStrategy;
};

  
/// Besides ElementLiouvillean, this is another solution based on the strategy idiom to control the number of Lindblads @ compile time
/**
 * \tparamRANK
 * \tparam NLINDBLADS the total number of Lindblads
 * \tparam IS_TIME_DEPENDENT governs time dependence
 * 
 * \see Sec. \ref basicoscillator of the structure-bundle guide for an example of usage 
 */
template<int RANK, int NLINDBLADS, bool IS_TIME_DEPENDENT=false>
class ElementLiouvilleanStrategies : public ElementLiouvilleanStrategiesBase<RANK,IS_TIME_DEPENDENT>
{
private:
  typedef ElementLiouvilleanStrategiesBase<RANK,IS_TIME_DEPENDENT> Base;
  
public:
  TYPE_DEFINITION_FORWARD
#undef TYPE_DEFINITION_FORWARD

  typedef typename Base::JumpStrategy JumpStrategy;
  typedef typename Base::JumpRateStrategy JumpRateStrategy;
  typedef typename Base::SuperoperatorStrategy SuperoperatorStrategy;

  /// Tiny vector of length NLINDBLADS containing the #JumpStrategy instances
  typedef blitz::TinyVector<JumpStrategy    ,NLINDBLADS> JumpStrategies;
  /// Tiny vector of length NLINDBLADS containing the #JumpRateStrategy instances
  typedef blitz::TinyVector<JumpRateStrategy,NLINDBLADS> JumpRateStrategies;

  typedef blitz::TinyVector<SuperoperatorStrategy,NLINDBLADS> SuperoperatorStrategies;

protected:
  template<typename... KeyLabelsPack>
  ElementLiouvilleanStrategies(const JumpStrategies& jumps, const JumpRateStrategies& jumpRates, const std::string& keyTitle, KeyLabelsPack&&... keyLabelsPack)
    : Base(keyTitle,keyLabelsPack...), jumps_(jumps), jumpRates_(jumpRates) {}
  
  ElementLiouvilleanStrategies(const JumpStrategies& jumps, const JumpRateStrategies& jumpRates, const std::string& keyTitle, typename Base::KeyLabelsInitializer il)
    : Base(keyTitle,il), jumps_(jumps), jumpRates_(jumpRates) {}

private:
  const Rates rates_v(Time t, const LazyDensityOperator& matrix) const final;

  void actWithJ_v(Time t, StateVectorLow& psi, size_t lindbladNo) const final;

  void actWithSuperoperator_v(Time t, const DensityOperatorLow& rho, DensityOperatorLow& drhodt, size_t lindbladNo) const final;
  
  const JumpStrategies     jumps_    ;
  const JumpRateStrategies jumpRates_;
  
  SuperoperatorStrategies superoperatorStrategies_;

};


} // structure


#endif // CPPQEDCORE_STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED


