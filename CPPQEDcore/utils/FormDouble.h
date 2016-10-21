// Copyright András Vukics 2006–2016. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_UTILS_FORMDOUBLE_H_INCLUDED
#define CPPQEDCORE_UTILS_FORMDOUBLE_H_INCLUDED

#include "FormDoubleFwd.h"

#include "Pars.h"

#include <algorithm>
#include <iosfwd>


/// Comprises tools related to FormDouble
/**
 * The aim of this module is to collect to one place all the issues of the formatting of doubles, necessary for output in various places.
 * The idea is from \cite stroustrup 21.4.6.3.
 * 
 * Named FormDouble here, since it is only issues pertaining to the display of doubles that are treated.
 */
namespace formdouble {

/// Output of a Bound value according to the corresponding FormDouble \related Bound
template<typename T>
std::ostream& operator<<(std::ostream&, const Bound<T>&);

} // formdouble


/// Class representing formatting options for output (stearming) of doubles
/**
 * It handles two characteristics, the \refStdCpp{precision,iomanip/setprecision} (number of digits) and \refStdCpp{width,iomanip/setw} of the output,
 * where the latter can also be left calculated automatically by FormDouble::FormDouble(int)
 */
class FormDouble 
{
public:
  static const int defaultPrecision; ///< Ultimate default precision for streaming doubles for the whole framework
  
  /// Generic overall precision accessible throughout the framework
  /** 
   * Can be set by reading in *any* variable of class FormDouble from the command line through the \link parameters parameter-bundle\endlink. 
   * \see parameters::Parameter<formdouble::Zero>, and trajectory::ParsRun::precision, where the overallPrecision is set by the \link trajectory\endlink-bundle
   */
  static       int overallPrecision;

  /// \name Constructors
  //@{
  /// Generic constructor
  FormDouble(int precision, int width) : precision_(precision), width_(width) {}
  
  /// Constructor calculating an appropriate width from the maximal possible width of a number of the given precision in the representation of the given machine.
  explicit FormDouble(int precision);
  //@}
  
  /// Binds a value of a double-based type (e.g. double, dcomp, `std::vector<double>`, etc.) to the present instant
  template<typename T>
  const formdouble::Bound<T> operator()(const T&) const;

  /// \name Getters
  //@{
  int getPrecision() const {return precision_;}
  int getWidth    () const {return     width_;}
  //@}
  
private:
  int precision_;
  int width_;

};



namespace formdouble {

/// If `precision` is larger than FormDouble::defaultPrecision, returns `precision`, otherwise the default \related FormDouble
inline int actualPrecision(int precision) {return std::max(FormDouble::defaultPrecision,precision);}

/// Essentially a compound of a FormDouble and a value of some double-based type (e.g. double, dcomp, `std::vector<double>`, etc.)
/** formdouble::Bound should be able to deal with everything consisting only of doubles (eg dcomp). */
template<typename T>
class Bound 
{
public:
  Bound(const FormDouble& ff, const T& v) : f(ff), val(v) {}

  const FormDouble& f;
  const T& val;

};


int widthPositive(int precision); ///< The maximal width in characters of a *positive* number streamed with the given `precision` \related FormDouble
int widthAny     (int precision); ///< The maximal width in characters of any (positive or negative) number streamed with the given `precision` \related FormDouble


inline const FormDouble low () {return FormDouble(FormDouble::defaultPrecision/2);} ///< Generic “low” precision (`FormDouble::defaultPrecision/2`) \related FormDouble
inline const FormDouble high() {return FormDouble(FormDouble::defaultPrecision  );} ///< Generic “high” precision (`FormDouble::defaultPrecision`) \related FormDouble


inline const FormDouble positive      (int precision) ///< FormDouble with at least FormDouble::defaultPrecision for positive-only quantities, such as time \related FormDouble
{return FormDouble(actualPrecision(precision),widthPositive(actualPrecision(precision)));}

inline const FormDouble zeroAdditional(int precision) ///< FormDouble with at least FormDouble::defaultPrecision and without additional width (padding) \related FormDouble
{return FormDouble(actualPrecision(precision),0);}


/// A dummy FormDouble used mainly for setting/getting the overallPrecision \see trajectory::ParsRun::precision.
/** Implicitely corvertible to/from an int */
class Zero : public FormDouble
{
public:
  Zero() : FormDouble(FormDouble::defaultPrecision) {}
  
  /// Implicit conversion from int possible.
  Zero(int precision) : FormDouble(precision,0) {}
  
  /// returns the precision in the style of implicit conversion
  operator int() const {return getPrecision();}

};



} // formdouble



namespace parameters {


/// Specialization which enables the setting of formdouble::FormDouble::overallPrecision in a read operation
template<>
class Parameter<formdouble::Zero> : public Parameter<int>
{
public:
  typedef Parameter<int> Base;

  Parameter(const std::string& s, const std::string& d, const formdouble::Zero& v) : Base(s,d,v), v_() {}

  const formdouble::Zero& get() const {return v_=formdouble::Zero(Base::get());}

  formdouble::Zero& get() {return const_cast<formdouble::Zero&>(static_cast<const Parameter*>(this)->get());}

protected:
  void read_v(std::istream& is) {Base::read_v(is); FormDouble::overallPrecision=get();}

private:
  mutable formdouble::Zero v_; // Just a buffer, the actual value is stored in Base

};


} // parameters


#endif // CPPQEDCORE_UTILS_FORMDOUBLE_H_INCLUDED
