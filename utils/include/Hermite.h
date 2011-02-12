// -*- C++ -*-
#ifndef _HERMITE_H
#define _HERMITE_H

#include "ComplexExtensions.h"

namespace cpputils {

struct HermiteOverflow {};

const double Hermite (size_t, double      );
const dcomp  HermiteC(size_t, const dcomp&);

}

#endif // _HERMITE_H
