// -*- C++ -*-
#ifndef STRUCTURE_ELEMENTAVERAGED_H_INCLUDED
#define STRUCTURE_ELEMENTAVERAGED_H_INCLUDED

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

private: 
  size_t       nAvr_v()                                                   const {return keyPrinter_.length()         ;}

  void      display_v(const Averages& a, std::ostream& os, int precision) const {       displayCommon(a,os,precision);}
  void   displayKey_v(std::ostream& os, size_t& i)                        const {       keyPrinter_.displayKey(os,i) ;}

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

class ReducedDensityOperator : public ClonableElementAveraged<1>
{
public:
  typedef ClonableElementAveraged<1> Base;

  ReducedDensityOperator(const std::string&, size_t, bool offDiagonals=false);

private:
  Base*const do_clone() const {return new ReducedDensityOperator(*this);}

  const Averages average_v(const LazyDensityOperator&) const;
  void           process_v(Averages&                 ) const {}

  const size_t dim_;
  const bool offDiagonals_;

};



template<int RANK, bool IS_TD>
class Collecting : public ClonableElementAveraged<RANK,IS_TD>
{
public:
  typedef ClonableElementAveraged<RANK,IS_TD> Element;
  typedef boost::ptr_list<Element> Collection;

  typedef ClonableElementAveraged<RANK,IS_TD> Base;
  typedef typename Base::KeyLabels KeyLabels;

  typedef AveragedCommon::Averages Averages;
  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  Collecting(const Collection&);
  Collecting(const Collecting&);

private:
  Base*const do_clone() const {return new Collecting(*this);}

  using Base::nAvr; using Base::getLabels;

  const Averages average_v(double, const LazyDensityOperator&) const;

  const Averages average_v(const LazyDensityOperator& matrix) const {return average_v(0.,matrix);}

  void  process_v(Averages&) const;

  const Collection collection_;

};



} // averaged


} // structure


#endif // STRUCTURE_ELEMENTAVERAGED_H_INCLUDED
