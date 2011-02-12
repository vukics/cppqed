// -*- C++ -*-
#ifndef _VIEWED_AS_MATRIX_H
#define _VIEWED_AS_MATRIX_H

#include "LazyDensityOperatorFwd.h"

#include "DimensionsBookkeeper.h"

#include "ComplexExtensions.h"

#include<boost/mpl/if.hpp>


namespace mpl=boost::mpl;


namespace quantumdata {



template<int RANK, typename F, typename V, typename T>
const T
partialTrace(const LazyDensityOperator<RANK>&, F, V, T);


/*
NOTE: The following definition does not work because implicit type-conversions (here, from a functor to a boost::function) are not considered for template parameter deduction. Nevertheless, Boost.ConceptCheck could help here.

template<int RANK, typename V, typename T>
const T
partialTrace(const LazyDensityOperator<RANK>&,
	     boost::function<const T(const typename ldo::DiagonalIterator<RANK,V>::value_type&)>, V, T);

*/


template<int RANK> 
class LazyDensityOperator 
  : public DimensionsBookkeeper<RANK,true>
{
public:
  typedef DimensionsBookkeeper<RANK,true> Base;

  typedef typename Base::Dimensions Dimensions;

  typedef typename mpl::if_c<(RANK==1),int,TTD_IDXTINY(RANK)>::type Idx;
  // Idx is just an int if RANK==1, otherwise a tiny of ints


  virtual ~LazyDensityOperator() {}


  virtual const dcomp operator()(const Idx& i, const Idx& j) const = 0;

  double operator()(const Idx& i) const {return real((*this)(i,i));}

  // iterator
  
  template<typename V>
  ldo::DiagonalIterator<RANK,V> begin(V) const;

  template<typename V>
  ldo::DiagonalIterator<RANK,V> end  (V) const;

protected:
  LazyDensityOperator(const Dimensions& dims) : Base(dims) {}

};


} // quantumdata

#endif // _VIEWED_AS_MATRIX_H
