// -*- C++ -*-
#ifndef   FORMDOUBLE_IMPL_H_INCLUDED
#define   FORMDOUBLE_IMPL_H_INCLUDED

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


#endif // FORMDOUBLE_IMPL_H_INCLUDED
