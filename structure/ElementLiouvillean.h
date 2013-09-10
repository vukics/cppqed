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
class Lindblad;


template<int RANK, int JUMP_ORDINAL, int NOJ>
class Lindblad<RANK,JUMP_ORDINAL,false,NOJ> : public mpl::if_c<JUMP_ORDINAL,Lindblad<RANK,JUMP_ORDINAL-1,false,NOJ>,LindbladBase<NOJ> >::type
{
public:
  virtual ~Lindblad () {}
  
  /// Type-erasure & non-virtual interface idiom in one:
  //@{
  void typeErasedActWithJ(typename quantumdata::Types<RANK>::StateVectorLow& psi) const {doActWithJ(psi,typename LindbladBase<NOJ>::template JumpNo<JUMP_ORDINAL>());}
  
  double typeErasedRate(const quantumdata::LazyDensityOperator<RANK>& matrix) const {return rate(matrix,typename LindbladBase<NOJ>::template JumpNo<JUMP_ORDINAL>());}
  //@}
  
private:
  virtual void doActWithJ(typename quantumdata::Types<RANK>::StateVectorLow&, typename LindbladBase<NOJ>::template JumpNo<JUMP_ORDINAL>) const = 0;
  
  virtual double rate(const quantumdata::LazyDensityOperator<RANK>&, typename LindbladBase<NOJ>::template JumpNo<JUMP_ORDINAL>) const = 0;
  
};


} // details


struct ElementLiouvilleanException : cpputils::TaggedException
{
  ElementLiouvilleanException(const std::string& tag) : cpputils::TaggedException(tag) {}
};

  
template<int RANK, int NOJ>
class ElementLiouvillean<RANK,NOJ,false> : public ElementLiouvilleanAveragedCommon<Liouvillean<RANK,false> >,
                                           public details::Lindblad<RANK,NOJ-1,false>
{
private:
  typedef ElementLiouvilleanAveragedCommon<Liouvillean<RANK,false> > Base;
  
public:
  typedef typename Base::StateVectorLow StateVectorLow;

  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  typedef typename Base::Rates Rates;
  
protected:
  template<typename... KeyLabelsPack>
  ElementLiouvillean(const std::string& keyTitle, KeyLabelsPack&&... keyLabelsPack) : Base(keyTitle,keyLabelsPack...) {}
  
  ElementLiouvillean(const std::string& keyTitle, typename Base::KeyLabelsInitializer il) : Base(keyTitle,il) {}

private:
  class Average
  {
  public:
    Average(const ElementLiouvillean* ptr, Rates& rates, const LazyDensityOperator& matrix) : ptr_(ptr), rates_(rates), matrix_(matrix) {}

    template<typename T> void operator()(T) const {rates_(T::value)=static_cast<const details::Lindblad<RANK,T::value,false,NOJ>*const>(ptr_)->typeErasedRate(matrix_);}

  private:
    const ElementLiouvillean*const ptr_;
    Rates& rates_;
    const LazyDensityOperator& matrix_;    
  };

  const Rates average_v(const LazyDensityOperator& matrix) const {Rates rates(NOJ); mpl::for_each<tmptools::Ordinals<NOJ> >(Average(this,rates,matrix)); return rates;}

  class ActWithJ
  {
  public:
    ActWithJ(const ElementLiouvillean* ptr, StateVectorLow& psi, size_t jumpNo) : ptr_(ptr), psi_(psi), jumpNo_(jumpNo) {}

    template<typename T>
    void operator()(T) const {if (T::value==jumpNo_) static_cast<const details::Lindblad<RANK,T::value,false,NOJ>*const>(ptr_)->typeErasedActWithJ(psi_);}

  private:
    const ElementLiouvillean*const ptr_;
    StateVectorLow& psi_;
    const size_t jumpNo_;    
  };

  void actWithJ_v(StateVectorLow& psi, size_t jumpNo) const {if (jumpNo>=NOJ) throw ElementLiouvilleanException(Base::getTitle()); mpl::for_each<tmptools::Ordinals<NOJ> >(ActWithJ(this,psi,jumpNo));}
  
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
  
protected:
  ElementLiouvillean(const std::string& keyTitle, const std::string& keyLabel) : Base(keyTitle,1,keyLabel) {}
  
private:
  const Rates average_v(const LazyDensityOperator& matrix) const {Rates rates(1); rates(0)=rate(matrix); return rates;}

  void actWithJ_v(StateVectorLow& psi, size_t jumpNo) const {if (jumpNo) throw ElementLiouvilleanException(Base::getTitle()); doActWithJ(psi);}

  virtual void doActWithJ(StateVectorLow&) const = 0;
  
  virtual double rate(const LazyDensityOperator&) const = 0;

};



} // structure


#endif // STRUCTURE_ELEMENTLIOUVILLEAN_H_INCLUDED


