// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \brief{Gathering to one place I/O formatting defaults used throughout the framework}
#ifndef CPPQEDCORE_UTILS_IO_MANIP_H_INCLUDED
#define CPPQEDCORE_UTILS_IO_MANIP_H_INCLUDED

#include "ComplexExtensions.h"

#include <boost/fusion/include/io.hpp>

#include <blitz/tinyvec2.cc>


class IO_Manipulator
{
public:
  static const char tuple_open='[', tuple_delimiter=',', tuple_close=']';
  
  template<typename STREAM>
  static auto& _(STREAM& stream) /// apply the manipulators
  {
    return stream<<boost::fusion::tuple_open(tuple_open)
                 <<boost::fusion::tuple_delimiter(tuple_delimiter)
                 <<boost::fusion::tuple_close(tuple_close);
  }
  
  template <typename T, int NL>
  static std::ostream& outputTinyVector(std::ostream& os, const blitz::TinyVector<T,NL>& x)
  {
    os<<tuple_open<<x[0];
    for (int i=1; i < NL; ++i) os<<tuple_delimiter<<x[i];
    return os<<tuple_close;
  }

  template <typename T, int NL>
  static std::istream& inputTinyVector(std::istream& is,       blitz::TinyVector<T,NL>& x)
  {
    using namespace std;

    char sep;
              
    is >> sep;
    BZPRECHECK(sep == tuple_open, "Format error while scanning input TinyVector" << endl << string(" (expected '") + tuple_open + string("' opening TinyVector)") );

    is >> x(0);
    for (int i = 1; i < NL; ++i) {
      if (!isspace(tuple_delimiter)) {
        is >> sep;
        BZPRECHECK(sep == tuple_delimiter, "Format error while scanning input TinyVector" << endl << string(" (expected '") + tuple_delimiter + string("' between TinyVector components)") );
      }
      BZPRECHECK(!is.bad(), "Premature end of input while scanning TinyVector");
      is >> x(i);
    }
    is >> sep;
    BZPRECHECK(sep == tuple_close, "Format error while scanning input TinyVector" << endl << string(" (expected '") + tuple_close + string("' closing TinyVector)") );
      
    return is;
  }

};


// unfortunately, blitz doesn’t take into account io manipulators, so we have to specialize the i/o operators here to be able to format the i/o
namespace blitz {

template<int NL> std::ostream& operator<<(std::ostream& s, const blitz::TinyVector<double,NL>& x) {return IO_Manipulator::outputTinyVector(s,x);}
template<int NL> std::istream& operator>>(std::istream& s,       blitz::TinyVector<double,NL>& x) {return IO_Manipulator:: inputTinyVector(s,x);}

template<int NL> std::ostream& operator<<(std::ostream& s, const blitz::TinyVector<dcomp,NL>& x) {return IO_Manipulator::outputTinyVector(s,x);}
template<int NL> std::istream& operator>>(std::istream& s,       blitz::TinyVector<dcomp,NL>& x) {return IO_Manipulator:: inputTinyVector(s,x);}

// alternative: (friend) function generating class

}

#endif // CPPQEDCORE_UTILS_IO_MANIP_H_INCLUDED
