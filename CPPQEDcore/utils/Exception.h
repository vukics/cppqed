// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines tentative base classes for the exception classes of the framework \todo Rationalize & uniformize the exceptions throughout}
#ifndef CPPQEDCORE_UTILS_EXCEPTION_H_INCLUDED
#define CPPQEDCORE_UTILS_EXCEPTION_H_INCLUDED

#include "ExceptionFwd.h"

#include <string>

namespace cpputils {

/// The class that is (meant to be, at least) the base of all exceptions in the framework
/**
 * \note It is extremely inconvenient to derive from \refStdCppConstruct{std::exception,exception/exception}, as this class defines `virtual ~exception() throw();`
 * so that each and every derived class also has to define its destructor in such a way. (Of course, this is just the point of `std::exception`, but here 
 * we do not strive for exception safety anyway as an exception usually means the end of the application.
 */
struct Exception // : public std::exception
{
  virtual ~Exception() /* throw() */ {}
};


/// Class reporting also the \link cpputils::TaggedException::what “what-ness”\endlink of the exception
class TaggedException : public Exception
{
public:
  TaggedException(const std::string& tag) : tag_(tag) {}

  virtual ~TaggedException() /* throw() */ {}

  const std::string getTag() const {return tag_;}

  /// Returns the short description
  const char* what() const throw() {return tag_.c_str();}

private:
  std::string tag_;

};


}

#endif // CPPQEDCORE_UTILS_EXCEPTION_H_INCLUDED
