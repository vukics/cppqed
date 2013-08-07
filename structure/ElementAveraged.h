// -*- C++ -*-
#ifndef STRUCTURE_ELEMENTAVERAGED_H_INCLUDED
#define STRUCTURE_ELEMENTAVERAGED_H_INCLUDED

#include "ElementAveragedFwd.h"

#include "Averaged.h"

#include "KeyPrinter.h"



namespace structure {


void displayCommon(const AveragedCommon::Averages&, std::ostream&, int);


//////////////////
//
// ElementAveraged
//
//////////////////


template<int RANK, bool IS_TIME_DEPENDENT>
class ElementAveraged : public Averaged<RANK,IS_TIME_DEPENDENT>
{
public:
  typedef AveragedCommon::Averages Averages;
  typedef cpputils::KeyPrinter::KeyLabels KeyLabels;

  ElementAveraged(const std::string& keyTitle, const KeyLabels& keyLabels) : keyPrinter_(keyTitle,keyLabels) {}

  const std::string& getTitle () const {return keyPrinter_.getTitle ();}
  const KeyLabels  & getLabels() const {return keyPrinter_.getLabels();}

private: 
  size_t       nAvr_v()                                                   const {return keyPrinter_.length()         ;}

  void             display_v(const Averages& a, std::ostream& os, int precision) const {       displayCommon(a,os,precision);}
  std::ostream& displayKey_v(std::ostream& os, size_t& i)                        const {return keyPrinter_.displayKey(os,i) ;}

  const cpputils::KeyPrinter keyPrinter_;

};



template<int RANK, bool IS_TIME_DEPENDENT>
class ClonableElementAveraged : public ElementAveraged<RANK,IS_TIME_DEPENDENT>
{
public:
  typedef ElementAveraged<RANK,IS_TIME_DEPENDENT> Base;
  typedef typename Base::KeyLabels KeyLabels;

  typedef ClonableElementAveraged* ClonedPtr;

  ClonableElementAveraged(const std::string& keyTitle, const KeyLabels& keyLabels) : Base(keyTitle,keyLabels) {}

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
