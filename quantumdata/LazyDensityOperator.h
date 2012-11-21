// -*- C++ -*-
#ifndef QUANTUMDATA_LAZYDENSITYOPERATOR_H_INCLUDED
#define QUANTUMDATA_LAZYDENSITYOPERATOR_H_INCLUDED

#include "LazyDensityOperatorFwd.h"

#include "DimensionsBookkeeper.h"

#include "ComplexExtensions.h"

#include "FFTFwd.h"

#include<boost/mpl/if.hpp>


namespace mpl=boost::mpl;


namespace quantumdata {



template<typename V, typename T, int RANK, typename F>
const T
partialTrace(const LazyDensityOperator<RANK>&, F);


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
  typedef boost::shared_ptr<const LazyDensityOperator> Ptr;
  
  typedef DimensionsBookkeeper<RANK,true> Base;

  typedef typename Base::Dimensions Dimensions;

  typedef typename mpl::if_c<(RANK==1),int,TTD_IDXTINY(RANK)>::type Idx;
  // Idx is just an int if RANK==1, otherwise a tiny of ints

  virtual ~LazyDensityOperator() {}

  const dcomp operator()(const Idx& i, const Idx& j) const {return index(i,j);}
  
  const Ptr ffTransform(fft::Direction dir) const {return ffTransform_v(dir);}

  double operator()(const Idx& i) const {return real((*this)(i,i));}

  // iterator
  
  template<typename V>
  const ldo::DiagonalIterator<RANK,V> begin() const;

  template<typename V>
  const ldo::DiagonalIterator<RANK,V> end  () const;

protected:
  LazyDensityOperator(const Dimensions& dims) : Base(dims) {}

private:
  virtual const dcomp index(const Idx& i, const Idx& j) const = 0;

  virtual const Ptr ffTransform_v(fft::Direction) const = 0;

};


} // quantumdata

#endif // QUANTUMDATA_LAZYDENSITYOPERATOR_H_INCLUDED
