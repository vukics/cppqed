// -*- C++ -*-
#ifndef STRUCTURE_AVERAGED_H_INCLUDED
#define STRUCTURE_AVERAGED_H_INCLUDED

#include "AveragedFwd.h"

#include "LiouvilleanAveragedCommon.h"

#include "Exception.h"

#include <iosfwd>


namespace structure {


/////////////////
//
// AveragedCommon
//
/////////////////


class AveragedCommon
{
public:
  typedef boost::shared_ptr<const AveragedCommon> Ptr;

  typedef LiouvilleanAveragedCommon::DArray1D Averages;

  virtual ~AveragedCommon() {}

  void process(Averages& averages) const {process_v(averages);}

  void display(const Averages& averages, std::ostream& os, int i) const {display_v(averages,os,i);}

private:
  virtual void process_v(      Averages&                    ) const = 0;
  virtual void display_v(const Averages&, std::ostream&, int) const = 0;

};


///////////
//
// Averaged
//
///////////



template<int RANK>
class Averaged<RANK,true>
  : public quantumdata::Types<RANK,LiouvilleanAveragedCommonRanked<RANK> >,public AveragedCommon
{
public:
  static const int N_RANK=RANK;

  typedef boost::shared_ptr<const Averaged> Ptr;

  typedef quantumdata::Types<RANK,LiouvilleanAveragedCommonRanked<RANK> > Base;

  typedef typename AveragedCommon::Averages Averages;

  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;

};


template<int RANK>
class Averaged<RANK,false>
  : public Averaged<RANK,true>
{
public:
  typedef typename Averaged<RANK,true>::Averages            Averages           ;
  typedef typename Averaged<RANK,true>::LazyDensityOperator LazyDensityOperator;

private:
  const Averages average_v(double, const LazyDensityOperator& matrix) const {return average_v(matrix);}

  virtual const Averages average_v(const LazyDensityOperator&) const = 0;

};


} // structure

#endif // STRUCTURE_AVERAGED_H_INCLUDED
