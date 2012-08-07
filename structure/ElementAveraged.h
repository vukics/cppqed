// -*- C++ -*-
#ifndef STRUCTURE_ELEMENT_AVERAGED_H_INCLUDED
#define STRUCTURE_ELEMENT_AVERAGED_H_INCLUDED

#include "ElementAveragedFwd.h"

#include "Averaged.h"

#include "KeyPrinter.h"

#include <boost/ptr_container/ptr_list.hpp>



namespace structure {


void displayCommon(const AveragedCommon::Averages&, std::ostream&, int);


//////////////////
//
// ElementAveraged
//
//////////////////


template<int RANK, bool IS_TD>
class ElementAveraged : public Averaged<RANK,IS_TD>
{
public:
  typedef AveragedCommon::Averages Averages;
  typedef cpputils::KeyPrinter::KeyLabels KeyLabels;

  ElementAveraged(const std::string& keyTitle, const KeyLabels& keyLabels) : keyPrinter_(keyTitle,keyLabels) {}

  const std::string& getTitle () const {return keyPrinter_.getTitle ();}
  const KeyLabels  & getLabels() const {return keyPrinter_.getLabels();}

  size_t nAvr()                                                      const {return keyPrinter_.length()         ;}

private: 
  void   display(const Averages& a, std::ostream& os, int precision) const {       displayCommon(a,os,precision);}
  void   displayKey(std::ostream& os, size_t& i)                     const {       keyPrinter_.displayKey(os,i) ;}

  const cpputils::KeyPrinter keyPrinter_;

};



template<int RANK, bool IS_TD>
class ClonableElementAveraged : public ElementAveraged<RANK,IS_TD>
{
public:
  typedef ElementAveraged<RANK,IS_TD> Base;
  typedef typename Base::KeyLabels KeyLabels;

  typedef ClonableElementAveraged* ClonedPtr;

  ClonableElementAveraged(const std::string& keyTitle, const KeyLabels& keyLabels) : Base(keyTitle,keyLabels) {}

  const ClonedPtr clone() const {return do_clone();}

private:
  virtual const ClonedPtr do_clone() const = 0;

};


template<int RANK, bool IS_TD>
inline ClonableElementAveraged<RANK,IS_TD>*const new_clone(const ClonableElementAveraged<RANK,IS_TD>& cea)
{
  return cea.clone();
}


namespace averaged {


// NEEDS_WORK generalize this for arbitrary RANK

class DiagonalDO : public ClonableElementAveraged<1>
{
public:
  typedef ClonableElementAveraged<1> Base;

  DiagonalDO(const std::string&, size_t);

private:
  Base*const do_clone() const {return new DiagonalDO(*this);}

  const Averages average(const LazyDensityOperator&) const;
  void           process(Averages&                 ) const {}

  const size_t dim_;

};



template<int RANK>
class Collecting : public ClonableElementAveraged<RANK>
{
public:
  typedef ClonableElementAveraged<RANK> Element;
  typedef boost::ptr_list<Element> Collection;

  typedef ClonableElementAveraged<RANK> Base;
  typedef typename Base::KeyLabels KeyLabels;

  typedef AveragedCommon::Averages Averages;
  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  Collecting(const Collection&);
  Collecting(const Collecting&);

private:
  Base*const do_clone() const {return new Collecting(*this);}

  using Base::nAvr; using Base::getLabels;

  const Averages average(const LazyDensityOperator&) const;
  void           process(Averages&                 ) const;

  const Collection collection_;

};



} // averaged


} // structure


#include "impl/ElementAveraged.tcc"


#endif // STRUCTURE_ELEMENT_AVERAGED_H_INCLUDED
