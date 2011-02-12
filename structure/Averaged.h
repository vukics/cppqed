// -*- C++ -*-
#ifndef _STRUCTURE_AVERAGED_H
#define _STRUCTURE_AVERAGED_H

#include "AveragedFwd.h"

#include "Types.h"

#include "Exception.h"

#include<iosfwd>


namespace structure {


/////////////////
//
// AveragedCommon
//
/////////////////


class AveragedCommon : public virtual LiouvilleanAveragedCommon
{
public:
  typedef TTD_DARRAY(1) Averages;


  static size_t nAvr(const AveragedCommon* averagedCommon)
  {
    return averagedCommon ? averagedCommon->nAvr() : 0;
  }


  static void displayKey(std::ostream& o, size_t& i, const AveragedCommon* averagedCommon)
  {
    if (averagedCommon) averagedCommon->displayKey(o,i);
  }


  static void process(Averages& averages, const AveragedCommon* averagedCommon) 
  {
    if (averagedCommon) averagedCommon->process(averages);
  }


  virtual ~AveragedCommon() {}


protected:
  virtual void   process   (Averages&)                           const = 0;
  virtual void   display   (const Averages&, std::ostream&, int) const = 0;

private:
  virtual void   displayKey(std::ostream&, size_t&) const = 0;
  virtual size_t nAvr      ()                       const = 0;

};


///////////
//
// Averaged
//
///////////


#ifndef   NDEBUG
struct AveragesFishyException : cpputils::Exception {};
#endif // NDEBUG


template<int RANK> class Averaged 
  : public quantumdata::Types<RANK,AveragedCommon>
{
public:
  typedef quantumdata::Types<RANK,AveragedCommon> Base;

  typedef AveragedCommon::Averages Averages;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

  using Base::display;

  static const Averages average(const LazyDensityOperator& matrix, const Averaged* averaged, StaticTag=theStaticOne)
  {
    const Averages averages(averaged ? averaged->average(matrix) : Averages());
#ifndef   NDEBUG
    if (size_t(averages.size())!=nAvr(averaged))
      throw AveragesFishyException();
#endif // NDEBUG
    return averages;
  }


  static void display(const LazyDensityOperator& matrix, std::ostream& os, int precision, const Averaged* averaged, StaticTag=theStaticOne)
  {
    if (averaged) {
      Averages averages(averaged->average(matrix));
      averaged->process(averages);
      averaged->display(averages,os,precision);
    }
  }


  virtual ~Averaged() {}

private:
  virtual const Averages average(const LazyDensityOperator&)     const = 0;

};



} // structure

#endif // _STRUCTURE_AVERAGED_H
