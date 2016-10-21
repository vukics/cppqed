// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
// -*- C++ -*-
#ifndef   CPPQEDCORE_UTILS_BOOLEANNEGATEDPROXY_H_INCLUDED
#define   CPPQEDCORE_UTILS_BOOLEANNEGATEDPROXY_H_INCLUDED

#include "BooleanNegatedProxyFwd.h"

#include <iosfwd>

namespace cpputils {


/// Bound to a boolean lvalue, it behaves like a boolean always with opposite value
class BooleanNegatedProxy
{
public:
  BooleanNegatedProxy(bool& v ///< The boolean lvalue to which the proxy is bound
                     ) : v_(v) {}
  
  /// Implicit conversion to boolean
  operator bool() const {return !v_;}
  
  BooleanNegatedProxy& operator=(const BooleanNegatedProxy& other) {v_=other.v_; return *this;}
  
  BooleanNegatedProxy& operator=(bool v) {v_=!v; return *this;}
  
private:
  bool &v_;
  
};


// std::ostream& operator<<(std::ostream&, const BooleanNegatedProxy&);

/// Input operator (output is trivial because of implicit conversion to boolean) \related BooleanNegatedProxy
std::istream& operator>>(std::istream&, BooleanNegatedProxy&);


} // cpputils

#endif // CPPQEDCORE_UTILS_BOOLEANNEGATEDPROXY_H_INCLUDED
