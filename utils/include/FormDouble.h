// -*- C++ -*-
/*
  Collecting to one place all the issues of the formatting of doubles, necessary for output in various places.
  The idea is from Stroustrup's The C++ Programming Language (special edition) 21.4.6.3.
*/
#ifndef _FORM_DOUBLE_H
#define _FORM_DOUBLE_H

#include "FormDoubleFwd.h"

#include<iosfwd>


namespace formdouble {

std::ostream& operator<<(std::ostream&, const Bound&);

}


class FormDouble 
{
public:
  explicit FormDouble(int precision);
  FormDouble(int precision, int width) : precision_(precision), width_(width) {}

  const formdouble::Bound operator()(double) const;

private:
  friend std::ostream& formdouble::operator<<(std::ostream&, const Bound&);
  
  int precision_;
  int width_;

};



namespace formdouble {


class Bound 
{
public:
  const FormDouble& f;
  double val;
  Bound(const FormDouble& ff, double v) : f(ff), val(v) {}
};


const FormDouble positive(int precision); // This is for quantities that can only be positive, such as time.


// Generic values:
inline const FormDouble low    () {return FormDouble(3);}
inline const FormDouble high   () {return FormDouble(6);}


} // formdouble



#endif // _FORM_DOUBLE_H
