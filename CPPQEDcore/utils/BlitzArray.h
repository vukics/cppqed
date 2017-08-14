// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines template aliases for real and complex arrays}
#ifndef   CPPQEDCORE_UTILS_BLITZARRAY_H_INCLUDED
#define   CPPQEDCORE_UTILS_BLITZARRAY_H_INCLUDED

#include "ComplexExtensions.h"

#include <blitz/array.h>


/// An array of doubles of arbitrary arity
template <int RANK> using DArray=blitz::Array<double,RANK>;

/// A complex array of arbitrary arity
template <int RANK> using CArray=blitz::Array<dcomp ,RANK>;

#endif // CPPQEDCORE_UTILS_BLITZARRAY_H_INCLUDED
