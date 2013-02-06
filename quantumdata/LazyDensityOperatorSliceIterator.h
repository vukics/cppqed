// -*- C++ -*-
#ifndef QUANTUMDATA_LAZYDENSITYOPERATORSLICEITERATOR_H_INCLUDED
#define QUANTUMDATA_LAZYDENSITYOPERATORSLICEITERATOR_H_INCLUDED

#include "LazyDensityOperatorFwd.h"

#include <boost/operators.hpp>
#include <boost/shared_ptr.hpp>

#include <boost/mpl/size.hpp>


namespace mpl=boost::mpl;


///////////////////////////////////
//
// LazyDensityOperatorSliceIterator
//
///////////////////////////////////

// it's inherently a const iterator since LazyDensityOperator is immutable

namespace quantumdata {

namespace ldo {

#define INPUT_IteratorHelper boost::input_iterator_helper<DiagonalIterator<RANK,V>,const LazyDensityOperator<mpl::size<V>::value> >

template<int RANK, typename V>
class DiagonalIterator 
  : public INPUT_IteratorHelper 
{
private:
  typedef INPUT_IteratorHelper Base;

#undef INPUT_IteratorHelper

  typedef typename Base::value_type LazyDensityOperatorRes;

public:
  template<bool IS_END>
  DiagonalIterator(const LazyDensityOperator<RANK>& ldo, mpl::bool_<IS_END>);
    
  DiagonalIterator& operator++() {impl_->doIncrement(); return *this;}

  const LazyDensityOperatorRes& operator*() const {return impl_->dereference();}
  
  bool isEqualTo(const DiagonalIterator& other) const {return impl_->isEqualTo(*other.impl_);}

  static const bool IS_SPECIAL=(RANK==mpl::size<V>::value);

  class DI_ImplSpecial;
  class DI_Impl       ;

  typedef boost::shared_ptr<typename mpl::eval_if_c<IS_SPECIAL,
                                                    mpl::identity<DI_ImplSpecial>,
                                                    mpl::identity<DI_Impl       > 
                                                    >::type
                           > Impl;
  // eval_if guarantees that DI_ImplSpecial gets instantiated only in the special case

private:
  Impl impl_;
  // We are using the pointer-to-implementation technique, together with some compile-time dispatching.
  
};


template<int RANK, typename V>
inline
bool operator==(const DiagonalIterator<RANK,V>& i1, const DiagonalIterator<RANK,V>& i2)
{
  return i1.isEqualTo(i2);
}


} // ldo

} // quantumdata


#endif // QUANTUMDATA_LAZYDENSITYOPERATORSLICEITERATOR_H_INCLUDED



