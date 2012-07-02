// -*- C++ -*-
#ifndef   CPPUTILS_KEY_PRINTER_INCLUDED
#define   CPPUTILS_KEY_PRINTER_INCLUDED

#include "KeyPrinterFwd.h"

#include <list>
#include <string>


namespace cpputils {

////////////////////////
//
// ElementAveragedCommon
//
////////////////////////


class KeyPrinter
{
public:
  typedef std::list<std::string> KeyLabels;

  KeyPrinter(const std::string&, const KeyLabels&);

  size_t length    ()                                                    const {return keyLabels_.size();}
  void   displayKey(std::ostream&, size_t&)                              const;

  const std::string& getTitle () const {return keyTitle_ ;}
  const KeyLabels  & getLabels() const {return keyLabels_;}

private:
  const std::string keyTitle_ ;  
  const KeyLabels   keyLabels_;
};


} // cpputils

#endif // CPPUTILS_KEY_PRINTER_INCLUDED
