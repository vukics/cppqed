#include "AveragingUtils.h"

#include "LazyDensityOperator.h"

#include "Algorithm.h"
#include "MathExtensions.h"

#include <boost/bind.hpp>

#include <boost/iterator/transform_iterator.hpp>

#include <boost/assign/list_of.hpp>

#include <algorithm>


using namespace std;

namespace {

struct Helper
{
  Helper(size_t dim) : dim_(dim), i_(1), j_(1), real_(true) {}

  const string operator()()
  {
    stringstream ss(stringstream::out);
    if (i_<dim_) {
      ss<<"rho_"<<i_<<','<<i_;
      ++i_;
    }
    else {
      size_t ii=i_-dim_;
      ss<<(real_ ? "real(" : "imag(")<<"rho_"<<ii<<','<<j_<<')';
      if (real_)
        real_=false;
      else {
        real_=true;
        if (j_<dim_-1)
          ++j_;
        else
          j_=++i_-dim_+1;
      }
    }
    return ss.str();
  }

private:
  const size_t dim_;
  size_t i_, j_;
  bool real_;
  
};

}


template<int RANK>
ReducedDensityOperator<RANK>::ReducedDensityOperator(const std::string& label, size_t dim, bool offDiagonals) : 
  Base(label,boost::assign::list_of(string("rho_0,0")).repeat_fun((offDiagonals ? mathutils::sqr(dim) : dim)-1,Helper(dim))),
  dim_(dim), offDiagonals_(offDiagonals)
{}


template<int RANK>
const typename ReducedDensityOperator<RANK>::Averages 
ReducedDensityOperator<RANK>::average_v(const LazyDensityOperator& matrix) const
{
  Averages averages(offDiagonals_ ? mathutils::sqr(dim_): dim_);
  for (size_t i=0; i<dim_; ++i)
    averages(i)=matrix(i);
  if (offDiagonals_)
    for (size_t i=0, idx=dim_; i<dim_; ++i)
      for (size_t j=i+1; j<dim_; ++j) {
        averages(idx++)=real(matrix(i,j));
        averages(idx++)=imag(matrix(i,j));
      }
  return averages;
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

