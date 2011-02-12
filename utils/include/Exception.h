// -*- C++ -*-
#ifndef _ERROR_H
#define _ERROR_H

#include<string>

namespace cpputils {

struct Exception {
  virtual ~Exception() {}
};


class TaggedException {
public:
  TaggedException(const std::string& tag) : tag_(tag) {}

  virtual ~TaggedException() {}

  const std::string getTag() const {return tag_;}

private:
  std::string tag_;

};


}

#endif // _ERROR_H
