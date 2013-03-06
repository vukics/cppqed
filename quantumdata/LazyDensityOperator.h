// -*- C++ -*-
#ifndef QUANTUMDATA_LAZYDENSITYOPERATOR_H_INCLUDED
#define QUANTUMDATA_LAZYDENSITYOPERATOR_H_INCLUDED

#include "LazyDensityOperatorFwd.h"

#include "DimensionsBookkeeper.h"

#include "BlitzArray.h"
#include "ComplexExtensions.h"

#include <boost/shared_ptr.hpp>

#include <boost/mpl/if.hpp>


namespace mpl=boost::mpl;


namespace quantumdata {



template<typename V, typename T, int RANK, typename F>
const T
partialTrace(const LazyDensityOperator<RANK>&, F);
// F must be an object callable with the signature:
// const T(const typename ldo::DiagonalIterator<RANK,V>::value_type&)



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

};


template<int RANK>
inline const typename LazyDensityOperator<RANK>::Idx dispatchLDO_index(const TTD_IDXTINY(RANK)& idx) {return idx   ;}

inline const          LazyDensityOperator<1   >::Idx dispatchLDO_index(const TTD_IDXTINY(1   )& idx) {return idx[0];}


template<int RANK>
const TTD_DARRAY(1) deflate(const LazyDensityOperator<RANK>&, bool offDiagonals);


} // quantumdata

#endif // QUANTUMDATA_LAZYDENSITYOPERATOR_H_INCLUDED
