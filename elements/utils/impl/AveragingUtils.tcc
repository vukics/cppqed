#ifndef   ELEMENTS_UTILS_IMPL_AVERAGINGUTILS_TCC_INCLUDED
#define   ELEMENTS_UTILS_IMPL_AVERAGINGUTILS_TCC_INCLUDED

#include "AveragingUtils.h"

#include "impl/DensityOperator.tcc"
#include "impl/LazyDensityOperator.tcc"
#include "impl/NegPT.tcc"

#include "Algorithm.h"
#include "MathExtensions.h"
#include "impl/MultiIndexIterator.tcc"
#include "Range.h"

#include <boost/bind.hpp>

#include <boost/iterator/transform_iterator.hpp>

#include <algorithm>


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
  Base(label,boost::assign::list_of(Helper()()).repeat_fun((offDiagonals ? mathutils::sqr(getTotalDimension()) : getTotalDimension())-1,Helper(getDimensions())).range(subsequent)),
  offDiagonals_(offDiagonals)
{
}


template<int RANK>
const typename ReducedDensityOperator<RANK>::Averages 
ReducedDensityOperator<RANK>::average_v(const LazyDensityOperator& matrix) const
{
  return deflate(matrix,offDiagonals_);
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

template<int RANK, bool IS_TD>
averagingUtils::Collecting<RANK,IS_TD>::Collecting(const Collection& collection)
  : Base(collection.begin()->getTitle(),
         cpputils::concatenateGrow(make_iterator_range(TRANSFORMED_iterator(begin),TRANSFORMED_iterator(end)),KeyLabels())),
    collection_(collection.clone())
{}


template<int RANK, bool IS_TD>
averagingUtils::Collecting<RANK,IS_TD>::Collecting(const Collecting& collecting)
  : Base(collecting),
    collection_(collecting.collection_.clone())
{}

#undef TRANSFORMED_iterator



#define TRANSFORMED_iterator(beginend) boost::make_transform_iterator(collection_.beginend(),boost::bind(&Element::average,_1,t,boost::cref(matrix)))

template<int RANK, bool IS_TD>
const typename averagingUtils::Collecting<RANK,IS_TD>::Averages
averagingUtils::Collecting<RANK,IS_TD>::average_v(double t, const LazyDensityOperator& matrix) const
{
  Averages res(nAvr()); res=0;
  return cpputils::concatenate(make_iterator_range(TRANSFORMED_iterator(begin),TRANSFORMED_iterator(end)),res);

}

#undef TRANSFORMED_iterator



template<int RANK, bool IS_TD>
void
averagingUtils::Collecting<RANK,IS_TD>::process_v(Averages& avr) const
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


#endif // ELEMENTS_UTILS_IMPL_AVERAGINGUTILS_TCC_INCLUDED