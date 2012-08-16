// -*- C++ -*-
#ifndef STRUCTURE_AVERAGED_H_INCLUDED
#define STRUCTURE_AVERAGED_H_INCLUDED

#include "AveragedFwd.h"

#include "LiouvilleanAveragedCommon.h"
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


class AveragedCommon : public LiouvilleanAveragedCommon
{
public:
  typedef DArray1D Averages;


  static size_t nAvr(const AveragedCommon* averagedCommon)
  {
    return averagedCommon ? averagedCommon->nAvr() : 0;
  }


  static void process(Averages& averages, const AveragedCommon* averagedCommon) 
  {
    if (averagedCommon) averagedCommon->process(averages);
  }


  virtual ~AveragedCommon() {}


  virtual void   process   (Averages&)                           const = 0;
  virtual void   display   (const Averages&, std::ostream&, int) const = 0;

  virtual size_t nAvr      ()                                    const = 0;

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


namespace averaged {


template<int RANKFROM, int RANKTO, bool IS_TD>
class Transferring : public Averaged<RANKFROM,IS_TD>
// Transfers the calculation of averages to another Averaged class, possibly with different RANK.
// The LazyDensityOperator for the other class should reference the same data.
// For a usage example cf. scripts/QbitMode_Matrix.cc
{
public:
  typedef typename Averaged<RANKFROM,IS_TD>::Averages            Averages           ;
  typedef typename Averaged<RANKFROM,IS_TD>::LazyDensityOperator LazyDensityOperator;

  typedef Averaged<RANKTO,IS_TD> AveragedTo;
  typedef quantumdata::LazyDensityOperator<RANKTO> LazyDensityOperatorTo;

  Transferring(const AveragedTo& averaged, const LazyDensityOperatorTo& ldo)
    : averaged_(averaged), ldo_(ldo) {}

  void process(Averages& averages) const {averaged_.process(averages);}

  void display(const Averages& averages, std::ostream& os, int n) const {averaged_.display(averages,os,n);}

  void displayKey(std::ostream& os, size_t& n) const {averaged_.displayKey(os,n);}

  size_t nAvr() const {return averaged_.nAvr();}

  const Averages average(double time, const LazyDensityOperator&) const {return averaged_.average(time,ldo_);}

  const Averages average(const LazyDensityOperator&) const {return averaged_.average(ldo_);}

private:
  const AveragedTo& averaged_;
  const LazyDensityOperatorTo& ldo_;

};


} // averaged



} // structure

#endif // STRUCTURE_AVERAGED_H_INCLUDED
