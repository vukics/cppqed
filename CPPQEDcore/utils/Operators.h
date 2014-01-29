/// \briefFile{Extensions built on top of \refBoost{Boost.Operator,utility/operators.htm}.}
// -*- C++ -*-
#ifndef UTILS_OPERATORS_H_INCLUDED
#define UTILS_OPERATORS_H_INCLUDED

#include "ComplexExtensions.h"

#include <boost/operators.hpp>

#include <boost/mpl/empty_base.hpp>


namespace linalg {


/// Operator aggregate for a complex vector space built on top of \refBoost{Boost.Operator,utility/operators.htm}.
/**
 * Comprises an Abel group and a complex outer product + comparison for equality
 * 
 * \note the real outer product must also be specified as there is no implicit conversion from `double` to dcomp
 * 
 * \tparam T the vector type
 */
template<typename T, typename B=boost::mpl::empty_base>
struct VectorSpace
  : boost::additive1     <T,        // Abel group
    boost::multiplicative<T,double, // Vector space
    boost::multiplicative<T,dcomp,  // "
    /*  boost::multipliable1  <T,        // Direct product */
    boost::equality_comparable<T> > > > {};


} // linalg

#endif // UTILS_OPERATORS_H_INCLUDED
