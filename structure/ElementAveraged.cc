#include "ElementAveraged.h"

#include "LazyDensityOperator.h"

#include "impl/FormDouble.tcc"
#include "MathExtensions.h"
#include "Range.h"

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include <boost/assign/list_of.hpp>

#include <algorithm>


using namespace std;



void structure::displayCommon(const AveragedCommon::Averages& averages, std::ostream& os, int precision)
{
  using namespace formdouble;

  namespace bll=boost::lambda;

  os<<'\t';
  {
    const FormDouble fd(precision);
    boost::for_each(averages,os<<bll::bind(&FormDouble::operator()<double>,&fd,bll::_1));
  }

}


namespace structure { namespace averaged {

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

ReducedDensityOperator::ReducedDensityOperator(const std::string& label, size_t dim, bool offDiagonals) : 
  Base(label,boost::assign::list_of(string("rho_0,0")).repeat_fun((offDiagonals ? mathutils::sqr(dim) : dim)-1,Helper(dim))),
  dim_(dim), offDiagonals_(offDiagonals)
{}


const ReducedDensityOperator::Averages 
ReducedDensityOperator::average_v(const LazyDensityOperator& matrix) const
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

} // averaged
} // structure
