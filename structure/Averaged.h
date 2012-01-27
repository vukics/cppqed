// -*- C++ -*-
#ifndef _STRUCTURE_AVERAGED_H
#define _STRUCTURE_AVERAGED_H

#include "AveragedFwd.h"

#include "Types.h"

#include "Exception.h"

#include "BlitzArrayExtensions.h"

#include <iosfwd>


namespace structure {


/////////////////
//
// AveragedCommon
//
/////////////////


class AveragedCommon : public virtual LiouvilleanAveragedCommon
{
public:
  typedef DArray1D Averages;


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


  virtual void   process   (Averages&)                           const = 0;
  virtual void   display   (const Averages&, std::ostream&, int) const = 0;

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

struct InfiniteDetectedException : cpputils::Exception {};


template<int RANK>
class Averaged<RANK,true>
  : public quantumdata::Types<RANK,AveragedCommon>
{
public:
  typedef quantumdata::Types<RANK,AveragedCommon> Base;

  typedef AveragedCommon::Averages Averages;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

  using Base::display;

  static const Averages average(double t, const LazyDensityOperator& matrix, const Averaged* averaged, StaticTag=theStaticOne)
  {
    return averaged ? averaged->average(t,matrix) : Averages();
  }


  static void display(double t, const LazyDensityOperator& matrix, std::ostream& os, int precision, const Averaged* averaged, StaticTag=theStaticOne)
  {
    if (averaged) {
      Averages averages(averaged->average(t,matrix));
#ifndef   NDEBUG
      if (size_t(averages.size())!=Base::nAvr(averaged)) throw AveragesFishyException();
#endif // NDEBUG
      if (!all(blitzplusplus::isfinite(averages))) throw InfiniteDetectedException();
      averaged->process(averages);
      averaged->display(averages,os,precision);
    }
  }


  virtual ~Averaged() {}

  virtual const Averages average(double, const LazyDensityOperator&) const = 0;

};


template<int RANK>
class Averaged<RANK,false>
  : public Averaged<RANK,true>
{
public:
  typedef typename Averaged<RANK,true>::Averages            Averages           ;
  typedef typename Averaged<RANK,true>::LazyDensityOperator LazyDensityOperator;

  virtual const Averages average(const LazyDensityOperator&) const = 0;

  const Averages average(double, const LazyDensityOperator& matrix) const {return average(matrix);}

};



} // structure

#endif // _STRUCTURE_AVERAGED_H
