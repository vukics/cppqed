// -*- C++ -*-
#ifndef   UTILS_INCLUDE_BOOLEANNEGATEDPROXY_H_INCLUDED
#define   UTILS_INCLUDE_BOOLEANNEGATEDPROXY_H_INCLUDED

#include "BooleanNegatedProxyFwd.h"

#include <iosfwd>

namespace cpputils {


class BooleanNegatedProxy
{
public:
  BooleanNegatedProxy(bool& v) : v_(v) {}
  
  operator bool() const {return !v_;}
  
  BooleanNegatedProxy& operator=(bool v) {v_=!v; return *this;}
  
private:
  bool &v_;
  
};


// std::ostream& operator<<(std::ostream&, const BooleanNegatedProxy&);

std::istream& operator>>(std::istream&,       BooleanNegatedProxy&);


} // cpputils

#endif // UTILS_INCLUDE_BOOLEANNEGATEDPROXY_H_INCLUDED
