// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#ifndef   CPPQEDELEMENTS_UTILS_AVERAGINGUTILS_H_INCLUDED
#define   CPPQEDELEMENTS_UTILS_AVERAGINGUTILS_H_INCLUDED

#include "ElementAveraged.h"

#include "DimensionsBookkeeper.h"

#include <boost/assign/list_of.hpp>

#include <boost/ptr_container/ptr_list.hpp>



template<int RANK>
class ReducedDensityOperator : private DimensionsBookkeeper<RANK>, public structure::ClonableElementAveraged<RANK>
{
private:
  typedef structure::ClonableElementAveraged<RANK> Base;

public:
  typedef typename Base::Averages Averages;
  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;
  
  typedef typename Base::KeyLabels KeyLabels;
  
  typedef typename DimensionsBookkeeper<RANK>::Dimensions Dimensions;
  
  using DimensionsBookkeeper<RANK>::getDimensions; using DimensionsBookkeeper<RANK>::getTotalDimension;

  ReducedDensityOperator(const std::string&, const Dimensions&, bool offDiagonals=false, const KeyLabels& =KeyLabels());

private:
  static const KeyLabels helper(const Dimensions&, bool offDiagonals, const KeyLabels&);

  Base*const do_clone() const {return new ReducedDensityOperator(*this);}

protected:
  const Averages average_v(structure::NoTime, const LazyDensityOperator&) const;

private:
  const bool offDiagonals_;

};


template<int RANK, typename V>
class ReducedDensityOperatorNegativity : public ReducedDensityOperator<RANK>
{
public:
  typedef ReducedDensityOperator<RANK> Base;
  typedef typename Base::Dimensions Dimensions;
  typedef typename Base::Averages Averages;
  typedef typename Base::LazyDensityOperator LazyDensityOperator;
  
  using Base::getDimensions; using Base::getTotalDimension;
  
  ReducedDensityOperatorNegativity(const std::string& label, const Dimensions& dim)
    : Base(label,dim,true,boost::assign::list_of("negativity")) {}

private:
  const Averages average_v(structure::NoTime, const LazyDensityOperator&) const;
  void           process_v(                   Averages&                 ) const;
  
};


namespace averagingUtils {


template<int RANK, bool IS_TIME_DEPENDENT=false>
class Collecting : public structure::ClonableElementAveraged<RANK,IS_TIME_DEPENDENT>
{
private:
  typedef structure::ClonableElementAveraged<RANK,IS_TIME_DEPENDENT> Base;

public:
  typedef Base Element;
  typedef boost::ptr_list<Element> Collection;

  typedef typename Base::KeyLabels KeyLabels;

  typedef structure::AveragedCommon::Averages Averages;
  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  typedef typename Base::Time Time;
  
  Collecting(const Collection&);
  Collecting(const Collecting&);

private:
  Base*const do_clone() const {return new Collecting(*this);}

  const Averages average_v(Time, const LazyDensityOperator&) const;

  void  process_v(Averages&) const;

  const Collection collection_;

};



template<int RANKFROM, int RANKTO, bool IS_TIME_DEPENDENT=false>
class Transferring : public structure::AveragedTimeDependenceDispatched<RANKFROM,IS_TIME_DEPENDENT>
// Transfers the calculation of averages to another Averaged class, possibly with different RANK.
// The LazyDensityOperator for the other class should reference the same data.
// For a usage example cf. scripts/QbitMode_Matrix.cc
{
private:
  typedef structure::AveragedTimeDependenceDispatched<RANKFROM,IS_TIME_DEPENDENT> Base;
  
public:
  typedef typename Base::Averages            Averages           ;
  typedef typename Base::LazyDensityOperator LazyDensityOperator;

  typedef typename Base::Time Time;
  
  typedef typename structure::AveragedTimeDependenceDispatched<RANKTO,IS_TIME_DEPENDENT>::Ptr AveragedToPtr;
  typedef quantumdata::LazyDensityOperator<RANKTO> LazyDensityOperatorTo;

  Transferring(AveragedToPtr averaged, const LazyDensityOperatorTo& ldo)
    : averaged_(averaged), ldo_(ldo) {}

private:
  std::ostream& streamKey_v(std::ostream& os, size_t& n) const {return averaged_->streamKey(os,n);}

  void process_v(Averages& averages) const {averaged_->process(averages);}

  std::ostream& stream_v(const Averages& averages, std::ostream& os, int n) const {return averaged_->stream(averages,os,n);}

  size_t nAvr_v() const {return averaged_->nAvr();}

  const Averages average_v(Time t, const LazyDensityOperator&) const;

  const AveragedToPtr averaged_;
  const LazyDensityOperatorTo& ldo_;

};


} // averagingUtils


#endif // CPPQEDELEMENTS_UTILS_AVERAGINGUTILS_H_INCLUDED
