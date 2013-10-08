// -*- C++ -*-
#ifndef   UTILS_INCLUDE_KEYPRINTER_H_INCLUDED
#define   UTILS_INCLUDE_KEYPRINTER_H_INCLUDED

#include "KeyPrinterFwd.h"

#include <initializer_list>
#include <list>
#include <string>
#include <utility>


namespace cpputils {


class KeyPrinter
{
public:
  typedef std::list<std::string> KeyLabels;
  typedef std::initializer_list<std::string> KeyLabelsInitializer;

  template<typename... KeyLabelsPack>
  KeyPrinter(const std::string& keyTitle, KeyLabelsPack&&... keyLabelsPack) : keyTitle_(keyTitle), keyLabels_(std::forward<KeyLabelsPack>(keyLabelsPack)...) {}

  KeyPrinter(const std::string& keyTitle, KeyLabelsInitializer il) : keyTitle_(keyTitle), keyLabels_(il) {}

  size_t        length    ()                       const {return keyLabels_.size();}
  std::ostream& displayKey(std::ostream&, size_t&) const;

  const std::string& getTitle () const {return keyTitle_ ;}
  const KeyLabels  & getLabels() const {return keyLabels_;}

private:
  const std::string keyTitle_ ;  
  const KeyLabels   keyLabels_;
};


} // cpputils

#endif // UTILS_INCLUDE_KEYPRINTER_H_INCLUDED
