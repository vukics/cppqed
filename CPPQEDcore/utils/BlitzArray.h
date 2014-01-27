// -*- C++ -*-
/// \briefFile{Defines template aliases for real and complex arrays}
#ifndef   UTILS_BLITZARRAY_H_INCLUDED
#define   UTILS_BLITZARRAY_H_INCLUDED

#include "ComplexExtensions.h"

#include <blitz/array.h>


/// An array of doubles of arbitrary arity
template <int RANK> using DArray=blitz::Array<double,RANK>;

/// A complex array of arbitrary arity
template <int RANK> using CArray=blitz::Array<dcomp ,RANK>;

#endif // UTILS_BLITZARRAY_H_INCLUDED
