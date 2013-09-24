// -*- C++ -*-
/// \briefFile{Defines the hierarchical partial specializations of structure::ElementLiouvillean}
#ifndef STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED
#define STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED

#include "ElementLiouvilleanFwd.h"

#include "Liouvillean.h"
#include "ElementLiouvilleanAveragedCommon.h"

#include <boost/mpl/for_each.hpp>


namespace structure {


namespace details {


template<int NOJ>
class LindbladBase
{
protected:
  /// A tagging class for Lindblad
  template<int JUMP_ORDINAL, typename OTHER=typename boost::enable_if_c< (JUMP_ORDINAL<NOJ) >::type>
  class JumpNo : mpl::int_<JUMP_ORDINAL> {};

};


template<int RANK, int JUMP_ORDINAL, bool IS_TIME_DEPENDENT, int NOJ=JUMP_ORDINAL+1>
class Lindblad : public mpl::if_c<JUMP_ORDINAL,Lindblad<RANK,JUMP_ORDINAL-1,IS_TIME_DEPENDENT,NOJ>,LindbladBase<NOJ> >::type
{
public:
  virtual ~Lindblad () {}

  typedef typename time::DispatcherIsTimeDependent<IS_TIME_DEPENDENT>::type Time;
  
  /// Type-erasure & non-virtual interface idiom in one:
  //@{
  void typeErasedActWithJ(Time t, typename quantumdata::Types<RANK>::StateVectorLow& psi) const {doActWithJ(t,psi,typename LindbladBase<NOJ>::template JumpNo<JUMP_ORDINAL>());}
  
  double typeErasedRate(Time t, const quantumdata::LazyDensityOperator<RANK>& matrix) const {return rate(t,matrix,typename LindbladBase<NOJ>::template JumpNo<JUMP_ORDINAL>());}
  //@}
  
private:
  virtual void doActWithJ(Time, typename quantumdata::Types<RANK>::StateVectorLow&, typename LindbladBase<NOJ>::template JumpNo<JUMP_ORDINAL>) const = 0;
  
  virtual double rate(Time, const quantumdata::LazyDensityOperator<RANK>&, typename LindbladBase<NOJ>::template JumpNo<JUMP_ORDINAL>) const = 0;
  
};


} // details


struct ElementLiouvilleanException : cpputils::TaggedException
{
  ElementLiouvilleanException(const std::string& tag) : cpputils::TaggedException(tag) {}
};

  
template<int RANK, int NOJ, bool IS_TIME_DEPENDENT>
class ElementLiouvillean : public ElementLiouvilleanAveragedCommon<LiouvilleanTimeDependenceDispatched<RANK,IS_TIME_DEPENDENT> >,
                           public details::Lindblad<RANK,NOJ-1,false>
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

    template<typename T> void operator()(T) const {rates_(T::value)=static_cast<const details::Lindblad<RANK,T::value,false,NOJ>*const>(ptr_)->typeErasedRate(t_,matrix_);}

  private:
    const ElementLiouvillean*const ptr_;
    Rates& rates_;
    const Time t_;
    const LazyDensityOperator& matrix_;    
  };

  const Rates average_v(Time t, const LazyDensityOperator& matrix) const {Rates rates(NOJ); mpl::for_each<tmptools::Ordinals<NOJ> >(Average(this,rates,t,matrix)); return rates;}

  class ActWithJ
  {
  public:
    ActWithJ(const ElementLiouvillean* ptr, Time t, StateVectorLow& psi, size_t jumpNo) : ptr_(ptr), t_(t), psi_(psi), jumpNo_(jumpNo) {}

    template<typename T>
    void operator()(T) const {if (T::value==jumpNo_) static_cast<const details::Lindblad<RANK,T::value,false,NOJ>*const>(ptr_)->typeErasedActWithJ(t_,psi_);}

  private:
    const ElementLiouvillean*const ptr_;
    const Time t_;
    StateVectorLow& psi_;
    const size_t jumpNo_;    
  };

  void actWithJ_v(Time t, StateVectorLow& psi, size_t jumpNo) const {if (jumpNo>=NOJ) throw ElementLiouvilleanException(Base::getTitle()); mpl::for_each<tmptools::Ordinals<NOJ> >(ActWithJ(this,t,psi,jumpNo));}
  
};



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
  const Rates average_v(Time t, const LazyDensityOperator& matrix) const {Rates rates(1); rates(0)=rate(t,matrix); return rates;}

  void actWithJ_v(Time t, StateVectorLow& psi, size_t jumpNo) const {if (jumpNo) throw ElementLiouvilleanException(Base::getTitle()); doActWithJ(t,psi);}

  virtual void doActWithJ(Time, StateVectorLow&) const = 0;
  
  virtual double rate(Time, const LazyDensityOperator&) const = 0;

};



} // structure


#endif // STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED


