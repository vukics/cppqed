// -*- C++ -*-
#ifndef UTILS_INCLUDE_OPERATORS_H_INCLUDED
#define UTILS_INCLUDE_OPERATORS_H_INCLUDED

#include "ComplexExtensions.h"

#include <boost/operators.hpp>

#include <boost/mpl/empty_base.hpp>


namespace linalg {


template<typename T, typename B=boost::mpl::empty_base>
struct VectorSpace
  : boost::additive1     <T,        // Abel group
    boost::multiplicative<T,double, // Vector space
    boost::multiplicative<T,dcomp,  // "
    /*  boost::multipliable1  <T,        // Direct product */
    boost::equality_comparable<T> > > > {};


} // linalg

#endif // UTILS_INCLUDE_OPERATORS_H_INCLUDED
