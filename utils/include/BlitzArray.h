// -*- C++ -*-
/// \briefFile{Defines template aliases for real and complex arrays}
#ifndef   UTILS_INCLUDE_BLITZARRAY_H_INCLUDED
#define   UTILS_INCLUDE_BLITZARRAY_H_INCLUDED

#include "ComplexExtensions.h"

#include <blitz/array.h>

#include <boost/utility.hpp>

#include <boost/mpl/int.hpp>


namespace blitzplusplus {


struct ShallowCopy {}; // For referencing constructors
struct DeepCopy    {}; // For copying     constructors


// A wrapper for blitz::Array::rank_; which did change from the CVS (where it used to be called _bz_rank) to the Mercurial version.
template<typename A>
struct ArrayRankTraits;

/** \cond */

template<typename T, int RANK>
struct ArrayRankTraits<blitz::Array<T,RANK> > : boost::mpl::int_<RANK> {};



// The semantics of the "DeepCopy" of Array is different from blitz::Array in that the storage order is NOT copied but rather a user-specified (by default the C) storage order is applied in the copy. I believe it is more in line with the fact that usually the way we think about the storage in Array is the C way.

template<typename T, int RANK>
class Array : public blitz::Array<T,RANK>, private boost::noncopyable
{
public:
  typedef blitz::Array<T,RANK> Base;
  typedef blitz::GeneralArrayStorage<RANK> GAS;

  Array(const Base& base, ShallowCopy) : Base(base) {}

  template<typename T_OTHER>
  Array(const blitz::Array<T_OTHER,RANK>& base, DeepCopy, GAS storage=GAS()) 
    : Base(base.shape(),storage) {*this=base;}

  using Base::shape; using Base::stride;

private:
  using Base::storage_;

public:
  const Base clone(T* restrict data) const {return Base(data,shape(),stride(),blitz::neverDeleteData,storage_);}

  const Base& getBase() const {return *this;}
  Base& getBase() {const_cast<Base&>(static_cast<const Array*>(this)->base());}

};

/** \endcond */

} // blitzplusplus


#define TTD_DARRAY(r) blitz::Array<double,r>
#define TTD_CARRAY(r) blitz::Array<dcomp ,r>
// TTD stands for "template typedef"

/// An array of doubles of arbitrary arity
template <int RANK> using DArray=blitz::Array<double,RANK>;

/// A complex array of arbitrary arity
template <int RANK> using CArray=blitz::Array<dcomp ,RANK>;

#endif // UTILS_INCLUDE_BLITZARRAY_H_INCLUDED
