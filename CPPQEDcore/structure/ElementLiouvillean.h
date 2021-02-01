// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines the structure::ElementLiouvillean & structure::ElementLiouvilleanStrategies}
#ifndef CPPQEDCORE_STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED
#define CPPQEDCORE_STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED

#include "Liouvillean.h"
#include "ElementLiouvilleanAveragedCommon.h"

#include <boost/range/algorithm/transform.hpp>

#include <boost/mpl/for_each.hpp>

#include <functional>


namespace mpl=boost::mpl;

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
  template<int LINDBLAD_ORDINAL, typename=std::enable_if_t< (LINDBLAD_ORDINAL<NLINDBLADS) > >
  using LindbladNo = tmptools::integral_c<LINDBLAD_ORDINAL>;

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
 * In the case of the _::typeErased… functions, the simultaneous use of the type-erasure & non-virtual interface idioms is notable
 * (on the other hand, we do override inherited non-virtual functions…)
 * 
 */
template<int RANK, int LINDBLAD_ORDINAL, bool IS_TIME_DEPENDENT, int NLINDBLADS>
class _ : public std::conditional_t< (LINDBLAD_ORDINAL>0) ,_<RANK,LINDBLAD_ORDINAL-1,IS_TIME_DEPENDENT,NLINDBLADS>,Base<NLINDBLADS> >
{
protected:
  typedef time::DispatcherIsTimeDependent_t<IS_TIME_DEPENDENT> Time;
  
  /// calls the virtual doActWithJ for the given LINDBLAD_ORDINAL
  void typeErasedActWithJ(Time t,
                          typename quantumdata::Types<RANK>::StateVectorLow& psi
                         ) const {doActWithJ(t,psi,tmptools::integral_c<LINDBLAD_ORDINAL>());}

  /// calls the virtual rate for the given LINDBLAD_ORDINAL
  double typeErasedRate(Time t,
                        const quantumdata::LazyDensityOperator<RANK>& matrix
                       ) const {return rate(t,matrix,tmptools::integral_c<LINDBLAD_ORDINAL>());}

  void typeErasedActWithSuperoperator(Time t,
                                      const typename quantumdata::Types<RANK>::DensityOperatorLow& rho,
                                      typename quantumdata::Types<RANK>::DensityOperatorLow& drhodt
                                     ) const {return doActWithSuperoperator(t,rho,drhodt,tmptools::integral_c<LINDBLAD_ORDINAL>());}

private:
  virtual void doActWithJ(Time, typename quantumdata::Types<RANK>::StateVectorLow&, tmptools::integral_c<LINDBLAD_ORDINAL>) const = 0;
  
  virtual double rate(Time, const quantumdata::LazyDensityOperator<RANK>&, tmptools::integral_c<LINDBLAD_ORDINAL>) const = 0;

  virtual void doActWithSuperoperator(Time,
                                      const typename quantumdata::Types<RANK>::DensityOperatorLow&,
                                      typename quantumdata::Types<RANK>::DensityOperatorLow&,
                                      tmptools::integral_c<LINDBLAD_ORDINAL>
                                    ) const {throw SuperoperatorNotImplementedException(LINDBLAD_ORDINAL);}
  
};


} // lindblad


/// Thrown if the Lindblad index is not smaller than the total number of Lindblads
struct ElementLiouvilleanException : std::runtime_error {using std::runtime_error::runtime_error;};


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
template<int RANK, int NLINDBLADS, bool IS_TIME_DEPENDENT=false>
class ElementLiouvillean : public ElementLiouvilleanAveragedCommon<LiouvilleanTimeDependenceDispatched<RANK,IS_TIME_DEPENDENT> >,
                           public lindblad::_<RANK,NLINDBLADS-1,false,NLINDBLADS>
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

  const Rates rates_v(Time t, const LazyDensityOperator& matrix) const override
  {
    Rates rates(NLINDBLADS);
/*
 * if 'this' is not captured, it’s an error ("'this' cannot be implicitly captured in this context"),
 * while if it’s captured, it’s a warning ("lambda capture 'this' is not used"), so we disable the warning
 */
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-lambda-capture"
    mpl::for_each<tmptools::Ordinals<NLINDBLADS> >([this,&rates,t,&matrix](auto arg) {
      using T = decltype(arg);
      rates(T::value)=lindblad::_<RANK,T::value,IS_TIME_DEPENDENT,NLINDBLADS>::typeErasedRate(t,matrix);
    });
    return rates;
  }


  void actWithJ_v(Time t, StateVectorLow& psi, size_t lindbladNo) const override
  {
    if (lindbladNo>=NLINDBLADS) throw ElementLiouvilleanException(Base::getTitle());
    mpl::for_each<tmptools::Ordinals<NLINDBLADS> >([this,t,&psi,lindbladNo](auto arg) {
      using T = decltype(arg);
      if (T::value==lindbladNo) lindblad::_<RANK,T::value,IS_TIME_DEPENDENT,NLINDBLADS>::typeErasedActWithJ(t,psi);
    });
  }
  

  void actWithSuperoperator_v(Time t, const DensityOperatorLow& rho, DensityOperatorLow& drhodt, size_t lindbladNo) const override
  {
    if (lindbladNo>=NLINDBLADS) throw ElementLiouvilleanException(Base::getTitle());
    mpl::for_each<tmptools::Ordinals<NLINDBLADS> >([this,t,&rho,&drhodt,lindbladNo](auto arg) {
      using T = decltype(arg);
      if (T::value==lindbladNo) lindblad::_<RANK,T::value,IS_TIME_DEPENDENT,NLINDBLADS>::typeErasedActWithSuperoperator(t,rho,drhodt);
    });
#pragma clang diagnostic pop
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
  
  const Rates rates_v(Time t, const LazyDensityOperator& matrix) const override {Rates rates(1); rates(0)=rate(t,matrix); return rates;}

  void actWithJ_v(Time t, StateVectorLow& psi, size_t lindbladNo) const override
  {
    if (lindbladNo) throw ElementLiouvilleanException(Base::getTitle());
    doActWithJ(t,psi);
  }

  void actWithSuperoperator_v(Time t, const DensityOperatorLow& rho, DensityOperatorLow& drhodt, size_t lindbladNo) const override
  {
    if (lindbladNo) throw ElementLiouvilleanException(Base::getTitle());
    doActWithSuperoperator(t,rho,drhodt);
  }
  
private:
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
  /** The actual signature of the functional is decided on the basis of IS_TIME_DEPENDENT */
  typedef std::conditional_t<IS_TIME_DEPENDENT,std::function<void  (double,       StateVectorLow&     )>,std::function<void  (      StateVectorLow&     )> > JumpStrategy;
  /// Strategy functional for calculating from a state the jump rate corresponding to a given Lindblad
  typedef std::conditional_t<IS_TIME_DEPENDENT,std::function<double(double, const LazyDensityOperator&)>,std::function<double(const LazyDensityOperator&)> > JumpRateStrategy;
  
  typedef std::conditional_t<IS_TIME_DEPENDENT,
                             std::function<void(double, const DensityOperatorLow&, DensityOperatorLow&)>,
                             std::function<void(        const DensityOperatorLow&, DensityOperatorLow&)> > SuperoperatorStrategy;
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

  const Rates rates_v(Time t, const LazyDensityOperator& matrix) const override
  {
    Rates rates(NLINDBLADS); // Note that this cannot be anything like static because of the by-reference semantics of blitz::Array

    boost::transform(jumpRates_,rates.begin(),
                     [&](const auto& jrs) -> double {
                       if constexpr (IS_TIME_DEPENDENT) return jrs(t,matrix);
                       else return jrs(matrix);
                     });
    
    return rates;
  }

  void actWithJ_v(Time t, StateVectorLow& psi, size_t lindbladNo) const override
  {
    if constexpr (IS_TIME_DEPENDENT) jumps_(lindbladNo)(t,psi);
    else jumps_(lindbladNo)(psi);
  }

  void actWithSuperoperator_v(Time t, const DensityOperatorLow& rho, DensityOperatorLow& drhodt, size_t lindbladNo) const override
  {
    if (superoperatorStrategies_(lindbladNo)==nullptr) throw SuperoperatorNotImplementedException(lindbladNo);

    if constexpr (IS_TIME_DEPENDENT) superoperatorStrategies_(lindbladNo)(t,rho,drhodt);
    else superoperatorStrategies_(lindbladNo)(rho,drhodt);
  }
  
  const JumpStrategies     jumps_    ;
  const JumpRateStrategies jumpRates_;
  
  SuperoperatorStrategies superoperatorStrategies_;

};


struct ElementLiouvilleanDiffusiveDimensionalityMismatchException : std::runtime_error {using std::runtime_error::runtime_error;};


/// Adds as many jump operators describing phase diffusion as the number of dimensions of the element
/**
 *  \tparam NLINDBLADS is now the number of additional Lindblads that is known @ compile time
 * 
 * In this way, any number of Liouvillean classes can be chained to compose arbitrarily complex Liouvillean evolutions.
 * For this, the implemented virtual functions must be at least protected.
 * 
 * This has to have access to the keys in order to be able to append, otherwise the ctor couldn’t be written for arbitrary base.
 * 
 * TODO: the class should check whether its base indeed has the functions that it relies on, otherwise replace with default funcionality
 * 
 */
template<int RANK, typename BASE, bool IS_TIME_DEPENDENT=false>
class ElementLiouvilleanDiffusive : public BASE
{
private:
  typedef BASE Base;
  
public:
  using DiffusionCoeffs=std::vector<double>;
  
  TYPE_DEFINITION_FORWARD
#undef TYPE_DEFINITION_FORWARD

  template <typename... BaseCtorPack>
  ElementLiouvilleanDiffusive(size_t dim, const DiffusionCoeffs& diffusionCoeffs, BaseCtorPack&&... baseCtorPack)
    : Base(std::forward<BaseCtorPack>(baseCtorPack)...), dim_(dim), base_nAvr_(this->getKeyPrinter().length()), diffusionCoeffs_(diffusionCoeffs)
  {
    for (size_t i=0; i<dim_; ++i) this->getKeyPrinter().getLabels().push_back("Phase flip for level "+std::to_string(i));
  }

  template <typename... BaseCtorPack>
  ElementLiouvilleanDiffusive(size_t dim, double diffusionCoeff, BaseCtorPack&&... baseCtorPack)
    : ElementLiouvilleanDiffusive(dim,DiffusionCoeffs(dim,diffusionCoeff),std::forward<BaseCtorPack>(baseCtorPack)...)
  {}

protected:
  const Rates rates_v(Time t, const LazyDensityOperator& matrix) const override
  {
    Rates res(base_nAvr_+dim_); res=-1.;
    res(blitz::Range(0,base_nAvr_-1))=Base::rates_v(t,matrix);
    return res;
  }

  void actWithJ_v(Time t, StateVectorLow& psi, size_t lindbladNo) const override
  {
    if (psi.extent(0)!=dim_) throw ElementLiouvilleanDiffusiveDimensionalityMismatchException("In actWithJ_v");
    if (lindbladNo<base_nAvr_) Base::actWithJ_v(t,psi,lindbladNo);
    else {
      psi*=sqrt(2.*diffusionCoeffs_[lindbladNo-base_nAvr_]);
      psi(lindbladNo-base_nAvr_)*=-1.;
    }
  }

  void actWithSuperoperator_v(Time t, const DensityOperatorLow& rho, DensityOperatorLow& drhodt, size_t lindbladNo) const override
  {
    if (rho.extent(0)!=dim_ || rho.extent(1)!=dim_) throw ElementLiouvilleanDiffusiveDimensionalityMismatchException("In actWithSuperoperator_v");
    if (lindbladNo<base_nAvr_) Base::actWithSuperoperator_v(t,rho,drhodt,lindbladNo);
    else if(lindbladNo==base_nAvr_) { // with this single case we cover the action of all the diffusive Lindblads
      for (int m=0; m<dim_; ++m) {
        drhodt(m,m)+=2*diffusionCoeffs_[lindbladNo-base_nAvr_]*rho(m,m);
        for (int n=m+1; n<dim_; ++n) {
          drhodt(m,n)-=6*diffusionCoeffs_[lindbladNo-base_nAvr_]*rho(m,n);
          drhodt(n,m)-=6*diffusionCoeffs_[lindbladNo-base_nAvr_]*rho(n,m);
        }
      }
    }
  }

private:
  const size_t dim_, base_nAvr_;
  const DiffusionCoeffs diffusionCoeffs_;
  
};


} // structure


#endif // CPPQEDCORE_STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED
