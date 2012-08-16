// -*- C++ -*-
#ifndef   UTILS_INCLUDE_BLITZTINY_H_INCLUDED
#define   UTILS_INCLUDE_BLITZTINY_H_INCLUDED

#include <blitz/array.h>

#include <boost/mpl/int.hpp>


namespace blitzplusplus {


// An indirection for accessing TinyVector's length at compile time.
template<typename V>
struct TinyVectorLengthTraits;

template<typename T, int LENGTH>
struct TinyVectorLengthTraits<blitz::TinyVector<T,LENGTH> > : boost::mpl::int_<LENGTH> {};


} // blitzplusplus


#define TTD_EXTTINY(r) blitz::TinyVector<   size_t,r>

#define TTD_IDXTINY(r) blitz::TinyVector<ptrdiff_t,r>


#endif // UTILS_INCLUDE_BLITZTINY_H_INCLUDED
