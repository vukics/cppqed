// -*- C++ -*-
#ifndef STRUCTURE_ELEMENTAVERAGED_H_INCLUDED
#define STRUCTURE_ELEMENTAVERAGED_H_INCLUDED

#include "ElementAveragedFwd.h"

#include "Averaged.h"
#include "ElementLiouvilleanAveragedCommon.h"



namespace structure {


std::ostream& displayCommon(const AveragedCommon::Averages&, std::ostream&, int);



template<int RANK, bool IS_TIME_DEPENDENT>
class ElementAveraged : public ElementLiouvilleanAveragedCommon<AveragedTimeDependenceDispatched<RANK,IS_TIME_DEPENDENT> >
{
private:
  typedef ElementLiouvilleanAveragedCommon<AveragedTimeDependenceDispatched<RANK,IS_TIME_DEPENDENT> > Base;
  
public:
  typedef AveragedCommon::Averages Averages;

protected:
  template<typename... KeyLabelsPack>
  ElementAveraged(const std::string& keyTitle, KeyLabelsPack&&... keyLabelsPack) : Base(keyTitle,keyLabelsPack...) {}
  
  ElementAveraged(const std::string& keyTitle, typename Base::KeyLabelsInitializer il) : Base(keyTitle,il) {}

private:
  std::ostream& display_v(const Averages& a, std::ostream& os, int precision) const {return displayCommon(a,os,precision);}

};



template<int RANK, bool IS_TIME_DEPENDENT>
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


#endif // STRUCTURE_ELEMENTAVERAGED_H_INCLUDED