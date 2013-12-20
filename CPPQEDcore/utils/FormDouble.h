// -*- C++ -*-
/*
  Collecting to one place all the issues of the formatting of doubles, necessary for output in various places.
  The idea is from Stroustrup's The C++ Programming Language (special edition) 21.4.6.3.

  Named FormDouble here, since it is only issues pertaining to the display of doubles that are treated. formdouble::Bound should be able to deal with everything consisting only doubles (eg dcomp).
*/

#ifndef UTILS_FORMDOUBLE_H_INCLUDED
#define UTILS_FORMDOUBLE_H_INCLUDED

#include "FormDoubleFwd.h"

#include "Pars.h"

#include <algorithm>
#include <iosfwd>


namespace formdouble {

template<typename T>
std::ostream& operator<<(std::ostream&, const Bound<T>&);

} // formdouble


class FormDouble 
{
public:
  static const int defaultPrecision;
  static       int overallPrecision;

  explicit FormDouble(int precision);
  // Calculates an appropriate width from the maximal possible width of a number of the given precision in the representation of the given machine.
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



class Zero : public FormDouble
{
public:
  Zero() : FormDouble(FormDouble::defaultPrecision) {}
  Zero(int precision) : FormDouble(precision,0) {}
  // Implicit conversion from int possible.

  operator double() const {return getPrecision();}

};



} // formdouble



namespace parameters {


template<>
class Parameter<formdouble::Zero> : public Parameter<int>
{
public:
  typedef Parameter<int> Base;

  Parameter(const std::string& s, const std::string& d, const formdouble::Zero& v) : Base(s,d,v), v_() {}

  void read(std::istream& is) {Base::read(is); FormDouble::overallPrecision=get();}

  const formdouble::Zero& get() const {return v_=formdouble::Zero(Base::get());}

  formdouble::Zero& get() {return const_cast<formdouble::Zero&>(static_cast<const Parameter*>(this)->get());}

private:
  mutable formdouble::Zero v_; // Just a buffer, the actual value is stored in Base

};


} // parameters


#endif // UTILS_FORMDOUBLE_H_INCLUDED
