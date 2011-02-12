// -*- C++ -*-
#ifndef   BOOLEAN_NEGATED_PROXY_INCLUDED
#define   BOOLEAN_NEGATED_PROXY_INCLUDED

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

#endif // BOOLEAN_NEGATED_PROXY_INCLUDED
