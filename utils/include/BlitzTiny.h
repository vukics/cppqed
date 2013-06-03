// -*- C++ -*-
/// Defines template aliases for `blitz::TinyVector`s used for characterising the size of multi-arrays and indexing them
#ifndef   UTILS_INCLUDE_BLITZTINY_H_INCLUDED
#define   UTILS_INCLUDE_BLITZTINY_H_INCLUDED

#include <blitz/array.h>

#include <boost/mpl/int.hpp>


namespace blitzplusplus {


// An indirection for accessing TinyVector's length at compile time.
template<typename V>
struct TinyVectorLengthTraits;

/** \cond */
template<typename T, int LENGTH>
struct TinyVectorLengthTraits<blitz::TinyVector<T,LENGTH> > : boost::mpl::int_<LENGTH> {};
/** \endcond */


} // blitzplusplus


#define TTD_EXTTINY(r) blitz::TinyVector<   size_t,r>
#define TTD_IDXTINY(r) blitz::TinyVector<ptrdiff_t,r>

/// A tiny vector describing extensions of objects of arbitrary arity
template <int RANK> using ExtTiny=blitz::TinyVector<   size_t,RANK>;

/// A tiny vector used for indexing of objects of arbitrary arity
template <int RANK> using IdxTiny=blitz::TinyVector<ptrdiff_t,RANK>;

#endif // UTILS_INCLUDE_BLITZTINY_H_INCLUDED
