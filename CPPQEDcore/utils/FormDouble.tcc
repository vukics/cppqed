// Copyright András Vukics 2006–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDCORE_UTILS_FORMDOUBLE_TCC_INCLUDED
#define   CPPQEDCORE_UTILS_FORMDOUBLE_TCC_INCLUDED

#include "FormDouble.h"

#include <sstream>


template<typename T>
std::ostream& formdouble::operator<<(std::ostream& os, const formdouble::Bound<T>& bf)
{
  using namespace std;
  ostringstream s; // intermediate stringstream taken so that the manips don't affect the NEXT output
  s.precision(bf.f.getPrecision());
  s.width(bf.f.getWidth());
  s.setf(ios_base::left,ios_base::adjustfield);
  s<<bf.val;
  return os<<s.str();
}


template<typename T>
const formdouble::Bound<T> FormDouble::operator()(const T& v) const
{
  return formdouble::Bound<T>(*this,v);
}


#endif // CPPQEDCORE_UTILS_FORMDOUBLE_TCC_INCLUDED
