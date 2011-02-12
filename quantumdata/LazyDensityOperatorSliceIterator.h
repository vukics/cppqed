// -*- C++ -*-
#ifndef _LAZY_DENSITY_OPERATOR_SMART_ITERATOR_H
#define _LAZY_DENSITY_OPERATOR_SMART_ITERATOR_H

#include "LazyDensityOperatorFwd.h"

#include<boost/operators.hpp>
#include<boost/shared_ptr.hpp>

#include<boost/mpl/size.hpp>


namespace mpl=boost::mpl;


///////////////////////////////////
//
// LazyDensityOperatorSliceIterator
//
///////////////////////////////////

// it's inherently a const iterator since LazyDensityOperator is immutable

namespace quantumdata {

namespace ldo {

namespace details {





template<typename>
class DI_ImplSpecial;

template<int, typename>
class DI_Impl;

} // details


#define TTD_LAZY_DENSITY_OPERATOR_RES const LazyDensityOperator<mpl::size<V>::value> 

#define TTD_FORWARD_ITERATOR_HELPER boost::input_iterator_helper<DiagonalIterator<RANK,V>,TTD_LAZY_DENSITY_OPERATOR_RES>

template<int RANK, typename V>
class DiagonalIterator 
  : public TTD_FORWARD_ITERATOR_HELPER 
{
public:
  static const bool IS_SPECIAL=(RANK==mpl::size<V>::value);

  typedef boost::shared_ptr<typename mpl::eval_if_c<IS_SPECIAL,
						    mpl::identity<details::DI_ImplSpecial<V> >,
						    mpl::identity<details::DI_Impl<RANK,V> > >::type> Impl;

  typedef TTD_FORWARD_ITERATOR_HELPER Base;

  typedef typename Base::value_type LazyDensityOperatorRes;

  template<bool IS_END>
  DiagonalIterator(const LazyDensityOperator<RANK>&, mpl::bool_<IS_END>);

  DiagonalIterator& operator++();

  const LazyDensityOperatorRes& operator*() const;
  
  bool isEqualTo(const DiagonalIterator&) const;

private:
  Impl impl_;
  // We are using the pointer-to-implementation technique, together with some compile-time dispatching.

};

#undef TTD_FORWARD_ITERATOR_HELPER

template<int RANK, typename V>
inline
bool operator==(const DiagonalIterator<RANK,V>& i1, const DiagonalIterator<RANK,V>& i2)
{
  return i1.isEqualTo(i2);
}


} // ldo

} // quantumdata


#include "impl/LazyDensityOperatorSliceIterator.tcc"

#undef TTD_LAZY_DENSITY_OPERATOR_RES

#endif // _LAZY_DENSITY_OPERATOR_SMART_ITERATOR_H



