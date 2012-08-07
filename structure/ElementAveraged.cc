#include "ElementAveraged.h"

#include "LazyDensityOperator.h"

#include "FormDouble.h"
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
  Helper() : i_(1) {}

  const string operator()()
  {
    stringstream ss(stringstream::out);
    ss<<"rho_"<<i_<<','<<i_;
    i_++;
    return ss.str();
  }

private:
  size_t i_;
};

}

DiagonalDO::DiagonalDO(const std::string& label, size_t dim) : 
  Base(label,boost::assign::list_of(string("rho_0,0")).repeat_fun(dim-1,Helper())),
  dim_(dim)
{}


const DiagonalDO::Averages 
DiagonalDO::average(const LazyDensityOperator& matrix) const
{
  Averages averages(dim_);
  for (size_t i=0; i<dim_; i++)
    averages(i)=matrix(i);
  return averages;
}

} // averaged
} // structure
