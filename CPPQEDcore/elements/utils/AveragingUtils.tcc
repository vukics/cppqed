#ifndef   ELEMENTS_UTILS_AVERAGINGUTILS_TCC_INCLUDED
#define   ELEMENTS_UTILS_AVERAGINGUTILS_TCC_INCLUDED

#include "AveragingUtils.h"

#include "DensityOperator.tcc"
#include "LazyDensityOperator.tcc"
#include "NegPT.tcc"

#include "Algorithm.h"
#include "MathExtensions.h"
#include "MultiIndexIterator.tcc"
#include "Range.h"

#include <boost/bind.hpp>
#include <boost/assign/std.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <algorithm>


template<int RANK>
const typename ReducedDensityOperator<RANK>::KeyLabels
ReducedDensityOperator<RANK>::helper(const Dimensions& dim, bool offDiagonals, const KeyLabels& subsequent)
{
  typedef cpputils::MultiIndexIterator<RANK> Iterator;
  using std::stringstream;

  KeyLabels res;
  Iterator i(Dimensions(size_t(0)),dim-1,cpputils::mii::begin);

  // Diagonals
  for (; i!=i.getEnd(); ++i) {
    stringstream ss(stringstream::out);
    ss<<"rho_"<<*i<<';'<<*i;
    res.push_back(ss.str());
  }

  if (!offDiagonals) return res;

  // Offdiagonals
  for (i.setToBegin(); i!=i.getEnd(); ++i) 
    for (Iterator j(std::next(i)); j!=j.getEnd(); ++j) {
      {
        stringstream ss(stringstream::out);
        ss<<"real["<<"rho_"<<*i<<','<<*j<<']';
        res.push_back(ss.str());
      }
      {
        stringstream ss(stringstream::out);
        ss<<"imag["<<"rho_"<<*i<<','<<*j<<']';
        res.push_back(ss.str());
      }
  }

  return res;

}


template<int RANK>
ReducedDensityOperator<RANK>::ReducedDensityOperator(const std::string& label, const Dimensions& dim, bool offDiagonals, const KeyLabels& subsequent) :
  DimensionsBookkeeper<RANK>(dim),
  Base(label,helper(getDimensions(),offDiagonals,subsequent)),
  offDiagonals_(offDiagonals)
{
}


template<int RANK>
const typename ReducedDensityOperator<RANK>::Averages 
ReducedDensityOperator<RANK>::average_v(structure::NoTime, const LazyDensityOperator& matrix) const
{
  return deflate(matrix,offDiagonals_);
}


template<int RANK, typename V>
const typename ReducedDensityOperatorNegativity<RANK,V>::Averages 
ReducedDensityOperatorNegativity<RANK,V>::average_v(structure::NoTime, const LazyDensityOperator& matrix) const
{
  Averages res(Base::average_v(matrix));
  res.resize(res.size()+1);
  return res;
}


template<int RANK, typename V>
void ReducedDensityOperatorNegativity<RANK,V>::process_v(Averages& averages) const
{
  quantumdata::DensityOperator<RANK> rho(getDimensions());
  linalg::CMatrix matrix(rho.matrixView()); // references the same data
  const size_t dim=getTotalDimension();
  
  {
    int idx=0;  
    for (int i=0; i<dim; ++i, ++idx) matrix(i,i)=averages(idx);
    for (int i=0; i<dim; ++i) for (int j=i+1; j<dim; ++j, idx+=2) matrix(j,i)=conj(matrix(i,j)=dcomp(averages(idx),averages(idx+1)));
  }
  
  averages(averages.size()-1)=quantumdata::negPT(rho,V());
  
}


#define TRANSFORMED_iterator(beginend) boost::make_transform_iterator(collection.beginend(),boost::bind(&Element::getLabels,_1))

template<int RANK, bool IS_TIME_DEPENDENT>
averagingUtils::Collecting<RANK,IS_TIME_DEPENDENT>::Collecting(const Collection& collection)
  : Base(collection.front().getTitle(),
         cpputils::concatenateGrow(make_iterator_range(TRANSFORMED_iterator(begin),TRANSFORMED_iterator(end)),KeyLabels())),
    collection_(collection.clone())
{}


template<int RANK, bool IS_TIME_DEPENDENT>
averagingUtils::Collecting<RANK,IS_TIME_DEPENDENT>::Collecting(const Collecting& collecting)
  : Base(collecting),
    collection_(collecting.collection_.clone())
{}

#undef TRANSFORMED_iterator


namespace averagingUtils { namespace details {

double convert(structure::OneTime t) {return t ;}
double convert(structure:: NoTime  ) {return 0.;}

} } // averagingUtils::details

#define TRANSFORMED_iterator(beginend) boost::make_transform_iterator(collection_.beginend(),boost::bind(&Element::average,_1,details::convert(t),boost::cref(matrix)))

template<int RANK, bool IS_TIME_DEPENDENT>
auto
averagingUtils::Collecting<RANK,IS_TIME_DEPENDENT>::average_v(Time t, const LazyDensityOperator& matrix) const -> const Averages
{
  Averages res(nAvr()); res=0;
  return cpputils::concatenate(make_iterator_range(TRANSFORMED_iterator(begin),TRANSFORMED_iterator(end)),res);

}

#undef TRANSFORMED_iterator



template<int RANK, bool IS_TIME_DEPENDENT>
void
averagingUtils::Collecting<RANK,IS_TIME_DEPENDENT>::process_v(Averages& avr) const
{
  struct Helper
  {
    static void doIt(const Element& eav, Averages& avr, ptrdiff_t& l, ptrdiff_t& u)
    {
      using blitz::Range;
      if ((u=l+eav.nAvr())>l) {
        Averages temp(avr(Range(l+1,u)));
        eav.process(temp);
      }
      std::swap(l,u);
    }
  };
 
  ptrdiff_t l=-1, u=0;
  for_each(collection_,boost::bind(Helper::doIt,_1,avr,l,u));
}


template<int RANKFROM, int RANKTO, bool IS_TIME_DEPENDENT>
auto
averagingUtils::Transferring<RANKFROM,RANKTO,IS_TIME_DEPENDENT>::average_v(Time t, const LazyDensityOperator&) const -> const Averages
{
  return averaged_->average(details::convert(t),ldo_);
}


/* The earlier solution fails with the C++11 standard due to a "bug" in Boost.Assign, cf. https://svn.boost.org/trac/boost/ticket/7364

template<int RANK>
struct ReducedDensityOperator<RANK>::Helper
{
  typedef cpputils::MultiIndexIterator<RANK> Iterator;
  
  Helper(const Dimensions& dim) : i_(Dimensions(size_t(0)),dim-1,cpputils::mii::begin), j_(i_), real_(true), offDiagonals_(false) {++i_; ++j_;}
  
  Helper() : i_(Dimensions(size_t(0)),Dimensions(size_t(0)),cpputils::mii::begin), j_(i_), real_(true), offDiagonals_(false) {}

  const std::string operator()()
  {
    using namespace std;
    stringstream ss(stringstream::out);

    // Diagonals
    if (!offDiagonals_ && i_!=i_.getEnd()) {
      ss<<"rho_"<<*i_<<';'<<*i_;
      ++i_;
    }
    else if (!offDiagonals_ && i_==i_.getEnd()) {
      offDiagonals_=true;
      i_.setToBegin();
    }
    
    // Offdiagonals
    if (offDiagonals_) {
      ss<<(real_ ? "real[" : "imag[")<<"rho_"<<*i_<<','<<*j_<<']';
      if (real_)
        real_=false;
      else {
        real_=true;
        ++j_;
        if (j_==j_.getEnd())
          ++(j_=++i_);
      }
    }
    
    return ss.str();
  }

private:
  Iterator i_, j_;
  bool real_, offDiagonals_;
  
};



template<int RANK>
ReducedDensityOperator<RANK>::ReducedDensityOperator(const std::string& label, const Dimensions& dim, bool offDiagonals, const KeyLabels& subsequent) :
  DimensionsBookkeeper<RANK>(dim),
  Base(label,boost::assign::repeat_fun((offDiagonals ? mathutils::sqr(getTotalDimension()) : getTotalDimension())-1,
                                       Helper(getDimensions())).range(subsequent)),
  offDiagonals_(offDiagonals)
{
}

*/

#endif // ELEMENTS_UTILS_AVERAGINGUTILS_TCC_INCLUDED