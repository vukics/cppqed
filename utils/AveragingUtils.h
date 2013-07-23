// -*- C++ -*-
#ifndef   ELEMENTS_UTILS_AVERAGINGUTILS_H_INCLUDED
#define   ELEMENTS_UTILS_AVERAGINGUTILS_H_INCLUDED

#include "AveragingUtilsFwd.h"

#include "ElementAveraged.h"

#include "DimensionsBookkeeper.h"

#include <boost/assign/list_of.hpp>

#include <boost/ptr_container/ptr_list.hpp>



template<int RANK>
class ReducedDensityOperator : private DimensionsBookkeeper<RANK>, public structure::ClonableElementAveraged<RANK>
{
public:
  typedef structure::ClonableElementAveraged<RANK> Base;
  typedef typename Base::Averages Averages;
  typedef quantumdata::LazyDensityOperator<RANK> LazyDensityOperator;
  
  typedef typename Base::KeyLabels KeyLabels;
  
  typedef typename DimensionsBookkeeper<RANK>::Dimensions Dimensions;
  
  ReducedDensityOperator(const std::string&, const Dimensions&, bool offDiagonals=false, const KeyLabels& =KeyLabels());

  using DimensionsBookkeeper<RANK>::getDimensions; using DimensionsBookkeeper<RANK>::getTotalDimension; using Base::nAvr;
  
private:
  struct Helper;

  Base*const do_clone() const {return new ReducedDensityOperator(*this);}

protected:
  const Averages average_v(const LazyDensityOperator&) const;

private:
  void           process_v(Averages&                 ) const {}

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
  const Averages average_v(const LazyDensityOperator&) const;
  void           process_v(Averages&                 ) const;
  
};


namespace averagingUtils {
  

template<int RANK, bool IS_TIME_DEPENDENT>
class Collecting : public structure::ClonableElementAveraged<RANK,IS_TIME_DEPENDENT>
{
public:
  typedef structure::ClonableElementAveraged<RANK,IS_TIME_DEPENDENT> Element;
  typedef boost::ptr_list<Element> Collection;

  typedef structure::ClonableElementAveraged<RANK,IS_TIME_DEPENDENT> Base;
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



template<int RANKFROM, int RANKTO, bool IS_TIME_DEPENDENT>
class Transferring : public structure::Averaged<RANKFROM,IS_TIME_DEPENDENT>
// Transfers the calculation of averages to another Averaged class, possibly with different RANK.
// The LazyDensityOperator for the other class should reference the same data.
// For a usage example cf. scripts/QbitMode_Matrix.cc
{
public:
  typedef typename structure::Averaged<RANKFROM,IS_TIME_DEPENDENT>::Averages            Averages           ;
  typedef typename structure::Averaged<RANKFROM,IS_TIME_DEPENDENT>::LazyDensityOperator LazyDensityOperator;

  typedef typename structure::Averaged<RANKTO,IS_TIME_DEPENDENT>::Ptr AveragedToPtr;
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
