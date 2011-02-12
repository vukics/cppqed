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


class FormDouble 
{
public:
  explicit FormDouble(int precision);
  FormDouble(int precision, int width) : precision_(precision), width_(width) {}

  const BoundFormDouble operator()(double) const;

private:
  friend std::ostream& operator<<(std::ostream&, const BoundFormDouble&);
  
  int precision_;
  int width_;

};


class BoundFormDouble 
{
public:
  const FormDouble& f;
  double val;
  BoundFormDouble(const FormDouble& ff, double v) : f(ff), val(v) {}
};


std::ostream& operator<<(std::ostream&, const BoundFormDouble&);


inline FormDouble low    () {return FormDouble(3   );}
inline FormDouble high   () {return FormDouble(6   );}
inline FormDouble special() {return FormDouble(6,13);} // This is for things that can only be positive, such as time.


} // formdouble



#endif // _FORM_DOUBLE_H
