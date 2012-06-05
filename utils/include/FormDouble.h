// -*- C++ -*-
/*
  Collecting to one place all the issues of the formatting of doubles, necessary for output in various places.
  The idea is from Stroustrup's The C++ Programming Language (special edition) 21.4.6.3.

  Named FormDouble here, since it is only issues pertaining to the display of doubles that is treated. formdouble::Bound should be able to deal with everything consisting only doubles (eg dcomp).
*/

#ifndef _FORM_DOUBLE_H
#define _FORM_DOUBLE_H

#include "FormDoubleFwd.h"

#include <algorithm>
#include <iosfwd>


namespace formdouble {

template<typename T>
std::ostream& operator<<(std::ostream&, const Bound<T>&);

}


class FormDouble 
{
public:
  static const int defaultPrecision;

  explicit FormDouble(int precision);
  FormDouble(int precision, int width) : precision_(precision), width_(width) {}

  template<typename T>
  const formdouble::Bound<T> operator()(const T&) const;

  int getPrecision() const {return precision_;}
  int getWidth    () const {return     width_;}

private:
  int precision_;
  int width_;

};



namespace formdouble {


inline int actualPrecision(int precision) {return std::max(FormDouble::defaultPrecision,precision);}


template<typename T>
class Bound 
{
public:
  Bound(const FormDouble& ff, const T& v) : f(ff), val(v) {}

  const FormDouble& f;
  const T& val;

};


int widthPositive(int precision);
int widthAny     (int precision);


// Generic values:
inline const FormDouble low () {return FormDouble(FormDouble::defaultPrecision/2);}
inline const FormDouble high() {return FormDouble(FormDouble::defaultPrecision  );}


// The following exhibit at least the defaultPrecision:
inline const FormDouble positive (int precision) // positive-only quantities, such as time
{return FormDouble(actualPrecision(precision),widthPositive(actualPrecision(precision)));}

inline const FormDouble zeroWidth(int precision) // without additional width
{return FormDouble(actualPrecision(precision),0);}


} // formdouble


#include "impl/FormDouble.tcc"


#endif // _FORM_DOUBLE_H
