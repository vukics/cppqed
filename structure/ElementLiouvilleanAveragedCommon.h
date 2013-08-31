// -*- C++ -*-
/// \briefFileDefault
#ifndef STRUCTURE_ELEMENTLIOUVILLEANAVERAGEDCOMMON_H_INCLUDED
#define STRUCTURE_ELEMENTLIOUVILLEANAVERAGEDCOMMON_H_INCLUDED

#include "ElementLiouvilleanAveragedCommonFwd.h"

#include "KeyPrinter.h"


namespace structure {


template<typename BASE>
class ElementLiouvilleanAveragedCommon : public BASE
{
public:
  typedef cpputils::KeyPrinter::KeyLabels KeyLabels;
  
  const std::string& getTitle () const {return keyPrinter_.getTitle ();}
  const KeyLabels  & getLabels() const {return keyPrinter_.getLabels();}

protected:
  template<typename... KeyLabelsPack>
  ElementLiouvilleanAveragedCommon(const std::string& keyTitle, KeyLabelsPack&&... keyLabelsPack) : keyPrinter_(keyTitle,keyLabelsPack...) {}

  ElementLiouvilleanAveragedCommon(const std::string& keyTitle, std::initializer_list<std::string> il) : keyPrinter_(keyTitle,il) {}

private:
  size_t nAvr_v() const {return keyPrinter_.length();}

  std::ostream& displayKey_v(std::ostream& os, size_t& i) const {return keyPrinter_.displayKey(os,i);}

  const cpputils::KeyPrinter keyPrinter_;

};

} // structure


#endif // STRUCTURE_ELEMENTLIOUVILLEANAVERAGEDCOMMON_H_INCLUDED


