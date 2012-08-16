// -*- C++ -*-
#ifndef UTILS_INCLUDE_HERMITE_H_INCLUDED
#define UTILS_INCLUDE_HERMITE_H_INCLUDED

#include "ComplexExtensions.h"

namespace cpputils {

struct HermiteOverflow {};

const double Hermite (size_t, double      );
const dcomp  HermiteC(size_t, const dcomp&);

}

#endif // UTILS_INCLUDE_HERMITE_H_INCLUDED
