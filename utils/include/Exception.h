// -*- C++ -*-
#ifndef EXCEPTION_HEADER_INCLUDED
#define EXCEPTION_HEADER_INCLUDED

#include "ExceptionFwd.h"

#include <string>

namespace cpputils {

struct Exception // : public std::exception
// It is extremely inconvenient to derive from std::exception, as this class defines
// virtual ~exception() throw();
// so that each and every derived class also has to define its destructor in such a way.
// Here, we do not strive for exception safety anyway.
{
  virtual ~Exception() /* throw() */ {}
};


class TaggedException : public Exception
{
public:
  TaggedException(const std::string& tag) : tag_(tag) {}

  virtual ~TaggedException() /* throw() */ {}

  const std::string getTag() const {return tag_;}

  const char* what() const throw() {return tag_.c_str();}

private:
  std::string tag_;

};


}

#endif // EXCEPTION_HEADER_INCLUDED
