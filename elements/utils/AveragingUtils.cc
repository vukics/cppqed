#include "AveragingUtils.h"

#include "LazyDensityOperator.h"

#include "Algorithm.h"
#include "MathExtensions.h"
#include "MultiIndexIterator.h"

#include <boost/bind.hpp>

#include <boost/iterator/transform_iterator.hpp>

#include <boost/assign/list_of.hpp>

#include <algorithm>


using namespace std;
using mathutils::sqr;


template<int RANK>
struct ReducedDensityOperator<RANK>::Helper
{
  typedef cpputils::MultiIndexIterator<RANK> Iterator;
  
  Helper(const Dimensions& dim) : i_(Dimensions(0ul),dim-1,Iterator::begin), j_(i_), real_(true), offDiagonals_(false) {++i_; ++j_;}
  
  Helper() : i_(Dimensions(0ul),Dimensions(0ul),Iterator::begin), j_(i_), real_(true), offDiagonals_(false) {}

  const string operator()()
  {
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
ReducedDensityOperator<RANK>::ReducedDensityOperator(const string& label, const Dimensions& dim, bool offDiagonals) :
  DimensionsBookkeeper<RANK>(dim),
  Base(label,boost::assign::list_of(Helper()()).repeat_fun((offDiagonals ? sqr(getTotalDimension()) : getTotalDimension())-1,Helper(getDimensions()))),
  offDiagonals_(offDiagonals)
{
/*  using namespace std; cerr<<dim<<endl;
  typedef typename Base::KeyLabels KeyLabels;
  const KeyLabels labels=Base::getLabels();
  for (typename KeyLabels::const_iterator i=labels.begin(); i!=labels.end(); ++i) cerr<<*i<<endl;*/
}


template<int RANK>
const typename ReducedDensityOperator<RANK>::Averages 
ReducedDensityOperator<RANK>::average_v(const LazyDensityOperator& matrix) const
{
  const size_t dim=getTotalDimension();
  Averages averages(offDiagonals_ ? sqr(dim) : dim);
  
  typedef cpputils::MultiIndexIterator<RANK> Iterator;
  const Iterator etalon(Dimensions(0ul),getDimensions()-1,Iterator::begin);
  
  size_t idx=0;

  for (Iterator i(etalon); idx<dim; ++i)
    averages(idx++)=matrix(quantumdata::dispatchLDO_index(*i));
  
  if (offDiagonals_)
    for (Iterator i=etalon.getBegin(); idx<sqr(dim); ++i)
      for (Iterator j=++Iterator(i); j!=etalon.getEnd(); ++j) {
        averages(idx++)=real(matrix(quantumdata::dispatchLDO_index(*i),quantumdata::dispatchLDO_index(*j)));
        averages(idx++)=imag(matrix(quantumdata::dispatchLDO_index(*i),quantumdata::dispatchLDO_index(*j)));
      }

  return averages;
  
}


template class ReducedDensityOperator<1>;
template class ReducedDensityOperator<2>;


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
      swap(l,u);
    }
  };
 
  ptrdiff_t l=-1, u=0;
  for_each(collection_,boost::bind(Helper::doIt,_1,avr,l,u));
}

