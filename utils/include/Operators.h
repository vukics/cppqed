// -*- C++ -*-
#ifndef UTILS_INCLUDE_OPERATORS_H_INCLUDED
#define UTILS_INCLUDE_OPERATORS_H_INCLUDED

#include "OperatorsFwd.h"

#include "ComplexExtensions.h"

#include<boost/operators.hpp>


namespace linalg {

namespace details {

struct EmptyBase {};

} // details

template<typename T, typename B>
struct VectorSpace
  : boost::additive1      <T,        // Abel group
    boost::multiplicative2<T,double, // Vector space
    boost::multiplicative2<T,dcomp,  // "
    /*  boost::multipliable1  <T,        // Direct product */
    boost::equality_comparable1<T> > > > {};


} // linalg

#endif // UTILS_INCLUDE_OPERATORS_H_INCLUDED
