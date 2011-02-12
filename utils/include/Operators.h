// -*- C++ -*-
#ifndef _CPP_UTILS_OPERATORS_H
#define _CPP_UTILS_OPERATORS_H

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

#endif // _CPP_UTILS_OPERATORS_H
