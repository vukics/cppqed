// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFileDefault
#ifndef   CPPQEDCORE_UTILS_KEYPRINTER_H_INCLUDED
#define   CPPQEDCORE_UTILS_KEYPRINTER_H_INCLUDED

#include <initializer_list>
#include <list>
#include <string>
#include <utility>


namespace cppqedutils {


/// Stores and prints a “key” (a.k.a. legend) to data, that is, an explanation to each element of a certain range of data
/**
 * The key is composed of a title and a list of labels
 * 
 * \see structure::LiouvilleanAveragedCommon::streamKey, structure::ElementLiouvilleanAveragedCommon
 */
class KeyPrinter
{
public:
  typedef std::list<std::string> KeyLabels; ///< the list of labels
  typedef std::initializer_list<std::string> KeyLabelsInitializer; ///< convenience typedef for initialization with an \refStdCppConstruct{initializer_list,initializer_list/initializer_list}

  /// construction from an arbitrary set of arguments capable to construct an object of type #KeyLabels
  template<typename... KeyLabelsPack>
  KeyPrinter(const std::string& keyTitle, KeyLabelsPack&&... keyLabelsPack) : keyTitle_(keyTitle), keyLabels_(std::forward<KeyLabelsPack>(keyLabelsPack)...) {}

  /// construction from a KeyLabelsInitializer list
  KeyPrinter(const std::string& keyTitle, KeyLabelsInitializer il) : keyTitle_(keyTitle), keyLabels_(il) {}

  size_t size() const {return keyLabels_.size();} ///< number of elements in the key
  std::ostream& stream(std::ostream&, size_t&) const; ///< streams the stored key in a nicely tabulated format

  /// \name Getters/setters
  //@{
  const std::string& getTitle () const {return keyTitle_ ;}
  const KeyLabels& getLabels() const {return keyLabels_;}
  KeyLabels& getLabels() {return const_cast<KeyLabels&>(static_cast<const KeyPrinter*>(this)->getLabels());}
  
  //@}
  
private:
  const std::string keyTitle_ ;  
  KeyLabels keyLabels_;
};


} // cppqedutils

#endif // CPPQEDCORE_UTILS_KEYPRINTER_H_INCLUDED
