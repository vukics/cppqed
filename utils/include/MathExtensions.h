// -*- C++ -*-
/*
  Wrapper functions for mathematical functions taken from libraries,
  and several other mathematical functions.

  The principal aim of Math.h and the corresponding Math.cc is to
  localize dependence on GSL.

*/

#ifndef UTILS_INCLUDE_MATHEXTENSIONS_H_INCLUDED
#define UTILS_INCLUDE_MATHEXTENSIONS_H_INCLUDED

#include "MathExtensionsFwd.h"

#include "ComplexExtensions.h"
#include "Exception.h"


namespace mathutils {

struct FactOverflow  : public cpputils::Exception {};
  
extern const double PI    ;
extern const double SQRTPI;
  
int sign(double);
int fcmp(double, double, double);
  
template<class T> const T sqr(T x) {return x*x;}

double sqr(double x);

double sqrAbs(const dcomp&);

double fact  (unsigned          ) throw(FactOverflow);
double choose(unsigned, unsigned);

bool parity(         long);
bool parity(unsigned long);
// even --- false, odd --- true;


inline int round(double r) {
  return (r>0) ? floor(r+0.5) : ceil(r-0.5);
}


inline char minusOneToThePowerOf(unsigned long n) {return parity(n) ? -1 : 1;}
inline char minusOneToThePowerOf(         long n) {return parity(n) ? -1 : 1;}

} // mathutils

#endif // UTILS_INCLUDE_MATHEXTENSIONS_H_INCLUDED
