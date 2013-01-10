// -*- C++ -*-
#ifndef   ELEMENTS_UTILS_AVERAGINGUTILS_H_INCLUDED
#define   ELEMENTS_UTILS_AVERAGINGUTILS_H_INCLUDED

#include "AveragingUtilsFwd.h"

#include "ElementAveraged.h"

#include "DimensionsBookkeeper.h"

#include <boost/ptr_container/ptr_list.hpp>



template<int RANK>
class ReducedDensityOperator : private DimensionsBookkeeper<RANK>, public structure::ClonableElementAveraged<RANK>
{
public:
  typedef structure::ClonableElementAveraged<RANK> Base;
  typedef typename Base::Averages Averages;
  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;
  
  typedef typename DimensionsBookkeeper<RANK>::Dimensions Dimensions;
  
  ReducedDensityOperator(const std::string&, const Dimensions&, bool offDiagonals=false);

private:
  struct Helper;

  using DimensionsBookkeeper<RANK>::getDimensions; using DimensionsBookkeeper<RANK>::getTotalDimension;
  
  Base*const do_clone() const {return new ReducedDensityOperator(*this);}

  const Averages average_v(const LazyDensityOperator&) const;
  void           process_v(Averages&                 ) const {}

  const bool offDiagonals_;

};


namespace averagingUtils {
  

template<int RANK, bool IS_TD>
class Collecting : public structure::ClonableElementAveraged<RANK,IS_TD>
{
public:
  typedef structure::ClonableElementAveraged<RANK,IS_TD> Element;
  typedef boost::ptr_list<Element> Collection;

  typedef structure::ClonableElementAveraged<RANK,IS_TD> Base;
  typedef typename Base::KeyLabels KeyLabels;

  typedef structure::AveragedCommon::Averages Averages;
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



template<int RANKFROM, int RANKTO, bool IS_TD>
class Transferring : public structure::Averaged<RANKFROM,IS_TD>
// Transfers the calculation of averages to another Averaged class, possibly with different RANK.
// The LazyDensityOperator for the other class should reference the same data.
// For a usage example cf. scripts/QbitMode_Matrix.cc
{
public:
  typedef typename structure::Averaged<RANKFROM,IS_TD>::Averages            Averages           ;
  typedef typename structure::Averaged<RANKFROM,IS_TD>::LazyDensityOperator LazyDensityOperator;

  typedef typename structure::Averaged<RANKTO,IS_TD>::Ptr AveragedToPtr;
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


} // averagingUtils


#endif // ELEMENTS_UTILS_AVERAGINGUTILS_H_INCLUDED
