// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef CPPQEDCORE_STRUCTURE_ELEMENTLIOUVILLEANAVERAGEDCOMMON_H_INCLUDED
#define CPPQEDCORE_STRUCTURE_ELEMENTLIOUVILLEANAVERAGEDCOMMON_H_INCLUDED

#include "KeyPrinter.h"


namespace structure {


/// Implements LiouvilleanAveragedCommon::streamKey and LiouvilleanAveragedCommon::nAvr with the help of a cppqedutils::KeyPrinter
/**
 * The number of averages is taken simply to be the \link cppqedutils::KeyPrinter::length length of the key\endlink
 */
template<typename BASE>
class ElementLiouvilleanAveragedCommon : public BASE
{
public:
  typedef cppqedutils::KeyPrinter::KeyLabels            KeyLabels           ;
  typedef cppqedutils::KeyPrinter::KeyLabelsInitializer KeyLabelsInitializer;
  
  const std::string& getTitle () const {return keyPrinter_.getTitle ();} ///< redirect to cppqedutils::KeyPrinter
  const KeyLabels  & getLabels() const {return keyPrinter_.getLabels();} ///< \copydoc getTitle

protected:
  /// \copydoc getTitle
  template<typename... KeyLabelsPack>
  ElementLiouvilleanAveragedCommon(const std::string& keyTitle, KeyLabelsPack&&... keyLabelsPack) : keyPrinter_(keyTitle,keyLabelsPack...) {}

  /// \copydoc getTitle
  ElementLiouvilleanAveragedCommon(const std::string& keyTitle, KeyLabelsInitializer il) : keyPrinter_(keyTitle,il) {}

  const cppqedutils::KeyPrinter& getKeyPrinter() const {return keyPrinter_;}
        cppqedutils::KeyPrinter& getKeyPrinter()       {return const_cast<cppqedutils::KeyPrinter&>(static_cast<const ElementLiouvilleanAveragedCommon*>(this)->getKeyPrinter());}
  
private:
  size_t nAvr_v() const final {return keyPrinter_.size();}

  std::ostream& streamKey_v(std::ostream& os, size_t& i) const final {return keyPrinter_.stream(os,i);}

  cppqedutils::KeyPrinter keyPrinter_;

};

} // structure


#endif // CPPQEDCORE_STRUCTURE_ELEMENTLIOUVILLEANAVERAGEDCOMMON_H_INCLUDED


