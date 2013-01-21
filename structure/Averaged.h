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

  std::ostream& display(const Averages& averages, std::ostream& os, int i) const {display_v(averages,os,i); return os;}

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

  typedef typename Averaged<RANKTO,IS_TD>::Ptr AveragedToPtr;
  typedef quantumdata::LazyDensityOperator<RANKTO> LazyDensityOperatorTo;

  Transferring(AveragedToPtr averaged, const LazyDensityOperatorTo& ldo)
    : averaged_(averaged), ldo_(ldo) {}

private:
  void displayKey_v(std::ostream& os, size_t& n) const {averaged_->displayKey(os,n);}

  void process_v(Averages& averages) const {averaged_->process(averages);}

  void display_v(const Averages& averages, std::ostream& os, int n) const {averaged_->display(averages,os,n);}

  size_t nAvr_v() const {return averaged_->nAvr();}

  const Averages average_v(double time, const LazyDensityOperator&) const {return averaged_->average(time,ldo_);}
  const Averages average_v(             const LazyDensityOperator&) const {return averaged_->average(     ldo_);}

  const AveragedToPtr averaged_;
  const LazyDensityOperatorTo& ldo_;

};


} // averaged



} // structure

#endif // STRUCTURE_AVERAGED_H_INCLUDED
