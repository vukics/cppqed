// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Defines structure::ElementAveraged & structure::ClonableElementAveraged}
#ifndef CPPQEDCORE_STRUCTURE_ELEMENTAVERAGED_H_INCLUDED
#define CPPQEDCORE_STRUCTURE_ELEMENTAVERAGED_H_INCLUDED

#include "Averaged.h"
#include "ElementLiouvilleanAveragedCommon.h"



namespace structure {


namespace details {

std::ostream& displayCommon(const AveragedCommon::Averages&, std::ostream&, int);

} // details

/// An implementation of Averaged for the case when the number of quantum arevages is known @ compile time (which is very often the case with elements)…
/**
 * … if this is not the case, Averaged has to be used instead (even for elements).
 *
 * \tparamRANK
 * \tparam IS_TIME_DEPENDENT governs time dependence
 *
 * Implements AveragedCommon::display in such a way that the averages stemming from the given element are displayed in a nicely tabulated format
 *
 * \see Sec. \ref basicoscillator of the structure-bundle guide for an example of usage
 *
 */
template<int RANK, bool IS_TIME_DEPENDENT=false>
class ElementAveraged : public ElementLiouvilleanAveragedCommon<AveragedTimeDependenceDispatched<RANK,IS_TIME_DEPENDENT> >
{
private:
  typedef ElementLiouvilleanAveragedCommon<AveragedTimeDependenceDispatched<RANK,IS_TIME_DEPENDENT> > Base;
  
public:
  typedef AveragedCommon::Averages Averages;

protected:
  /// \name Constructors
  //@{
  template<typename... KeyLabelsPack>
  ElementAveraged(const std::string& keyTitle, KeyLabelsPack&&... keyLabelsPack) : Base(keyTitle,keyLabelsPack...) {} ///< The number of KeyLabel arguments in the constructors determines the number of calculated averages.
  
  ElementAveraged(const std::string& keyTitle, typename Base::KeyLabelsInitializer il) : Base(keyTitle,il) {} ///< ”
  //@}

  const Averages initializedAverages() const {Averages res(this->nAvr()); res=0.; return res;}

private:
  std::ostream& display_v(const Averages& a, std::ostream& os, int precision) const final {return details::displayCommon(a,os,precision);}

};



/// Besides being an ElementAveraged, it models the \refBoost{Clonable concept,ptr_container/doc/reference.html#the-clonable-concept}
template<int RANK, bool IS_TIME_DEPENDENT=false>
class ClonableElementAveraged : public ElementAveraged<RANK,IS_TIME_DEPENDENT>
{
private:
  typedef ElementAveraged<RANK,IS_TIME_DEPENDENT> Base;

protected:
  template<typename... KeyLabelsPack>
  ClonableElementAveraged(const std::string& keyTitle, KeyLabelsPack&&... keyLabelsPack) : Base(keyTitle,keyLabelsPack...) {}

  ClonableElementAveraged(const std::string& keyTitle, typename Base::KeyLabelsInitializer il) : Base(keyTitle,il) {}

public:
  typedef ClonableElementAveraged* ClonedPtr;

  const ClonedPtr clone() const {return do_clone();}

private:
  virtual const ClonedPtr do_clone() const = 0;

};


template<int RANK, bool IS_TIME_DEPENDENT>
inline ClonableElementAveraged<RANK,IS_TIME_DEPENDENT>*const new_clone(const ClonableElementAveraged<RANK,IS_TIME_DEPENDENT>& cea)
{
  return cea.clone();
}


} // structure


#endif // CPPQEDCORE_STRUCTURE_ELEMENTAVERAGED_H_INCLUDED
