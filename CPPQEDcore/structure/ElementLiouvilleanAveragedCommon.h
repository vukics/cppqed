// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_STRUCTURE_ELEMENTLIOUVILLEANAVERAGEDCOMMON_H_INCLUDED
#define CPPQEDCORE_STRUCTURE_ELEMENTLIOUVILLEANAVERAGEDCOMMON_H_INCLUDED

#include "ElementLiouvilleanAveragedCommonFwd.h"

#include "KeyPrinter.h"


namespace structure {


/// Implements LiouvilleanAveragedCommon::displayKey and LiouvilleanAveragedCommon::nAvr with the help of a cpputils::KeyPrinter
/**
 * The number of averages is taken simply to be the \link cpputils::KeyPrinter::length length of the key\endlink
 */
template<typename BASE>
class ElementLiouvilleanAveragedCommon : public BASE
{
public:
  typedef cpputils::KeyPrinter::KeyLabels            KeyLabels           ;
  typedef cpputils::KeyPrinter::KeyLabelsInitializer KeyLabelsInitializer;
  
  const std::string& getTitle () const {return keyPrinter_.getTitle ();} ///< redirect to cpputils::KeyPrinter
  const KeyLabels  & getLabels() const {return keyPrinter_.getLabels();} ///< \copydoc getTitle

protected:
  /// \copydoc getTitle
  template<typename... KeyLabelsPack>
  ElementLiouvilleanAveragedCommon(const std::string& keyTitle, KeyLabelsPack&&... keyLabelsPack) : keyPrinter_(keyTitle,keyLabelsPack...) {}

  /// \copydoc getTitle
  ElementLiouvilleanAveragedCommon(const std::string& keyTitle, KeyLabelsInitializer il) : keyPrinter_(keyTitle,il) {}

private:
  size_t nAvr_v() const final {return keyPrinter_.length();}

  std::ostream& displayKey_v(std::ostream& os, size_t& i) const final {return keyPrinter_.displayKey(os,i);}

  const cpputils::KeyPrinter keyPrinter_;

};

} // structure


#endif // CPPQEDCORE_STRUCTURE_ELEMENTLIOUVILLEANAVERAGEDCOMMON_H_INCLUDED


