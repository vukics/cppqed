// -*- C++ -*-
#ifndef _STRUCTURE_ELEMENT_AVERAGED_H
#define _STRUCTURE_ELEMENT_AVERAGED_H

#include "ElementAveragedFwd.h"

#include "Averaged.h"

#include<list>


namespace structure {


////////////////////////
//
// ElementAveragedCommon
//
////////////////////////


class ElementAveragedCommon
{
public:
  typedef std::list<std::string> KeyLabels;

  ElementAveragedCommon(const std::string&, const KeyLabels&);

  void   display   (const AveragedCommon::Averages&, std::ostream&, int) const;
  size_t nAvr      ()                                                    const {return keyLabels_.size();}
  void   displayKey(std::ostream&, size_t&)                              const;

private:
  const std::string keyTitle_ ;  
  const KeyLabels   keyLabels_;
};


//////////////////
//
// ElementAveraged
//
//////////////////


template<int RANK>
class ElementAveraged : public Averaged<RANK>, private ElementAveragedCommon
{
public:
  typedef AveragedCommon::Averages Averages;
  typedef ElementAveragedCommon::KeyLabels KeyLabels;

  ElementAveraged(const std::string& keyTitle, const KeyLabels& keyLabels) : ElementAveragedCommon(keyTitle,keyLabels) {}

private: 
  void   display(const Averages& a, std::ostream& os, int precision) const {       ElementAveragedCommon::display(a,os,precision);}
  size_t nAvr()                                                      const {return ElementAveragedCommon::nAvr()                 ;}
  void   displayKey(std::ostream& os, size_t& i)                     const {       ElementAveragedCommon::displayKey(os,i)       ;}

};



namespace averaged {


// NEEDS_WORK generalize this for arbitrary RANK

class DiagonalDO : public structure::ElementAveraged<1>
{
public:
  typedef structure::ElementAveraged<1> Base;

  DiagonalDO(size_t);

private:
  const Averages average(const LazyDensityOperator&) const;
  void           process(Averages&                 ) const {}

  const size_t dim_;

};


} // averaged


} // structure


#endif // _STRUCTURE_ELEMENT_AVERAGED_H
