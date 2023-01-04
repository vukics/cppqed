// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Extensions built on top of \refBoost{Boost.Operator,utility/operators.htm}.}
#pragma once

#include "ComplexExtensions.h"

#include <boost/operators.hpp>

#include <optional>


namespace linalg {


/// Operator aggregate for a complex vector space built on top of \refBoost{Boost.Operator,utility/operators.htm}.
/**
 * Comprises an Abel group and a complex outer product + comparison for equality
 * 
 * \note the real outer product must also be specified as there is no implicit conversion from `double` to dcomp
 * 
 * \tparam T the vector type
 */
template<typename T, typename Base = std::nullopt_t>
struct VectorSpace
  : boost::additive1     <T,        // Abel group
    boost::multiplicative<T,double, // Vector space
    boost::multiplicative<T,dcomp,  // "
    /*  boost::multipliable1  <T,        // Direct product */
    boost::equality_comparable<T> > > >, Base {};


} // linalg


