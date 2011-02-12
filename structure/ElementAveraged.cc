#include "ElementAveraged.h"

#include "LazyDensityOperator.h"

#include "FormDouble.h"

#include "Range.h"

#include<boost/lambda/lambda.hpp>
#include<boost/lambda/bind.hpp>

#include<boost/assign/list_of.hpp>

#include<algorithm>


using namespace boost::assign;
using namespace std;
namespace bll=boost::lambda;

using boost::for_each;


structure::ElementAveragedCommon::ElementAveragedCommon(const std::string& keyTitle, const KeyLabels& keyLabels) 
  : keyTitle_(keyTitle), keyLabels_(keyLabels) 
{
}



void structure::ElementAveragedCommon::displayKey(std::ostream& os, size_t& i) const
{
  os<<"# "<<keyTitle_;
  for_each(keyLabels_,os<<bll::constant("\n# ")<<bll::constant(setw(2))<<bll::var(i)++<<". "<<bll::_1);
  os<<endl;
}



void structure::ElementAveragedCommon::display(const AveragedCommon::Averages& averages, std::ostream& os, int precision) const
{
  using namespace cpputils;
  using namespace formdouble;

  namespace bll=boost::lambda;

  os<<'\t';
  {
    const FormDouble fd(precision);
    for_each(averages,os<<bll::bind(&FormDouble::operator(),&fd,bll::_1));
  }

}


namespace structure { namespace averaged {

namespace details {

struct Helper
{
  Helper() : i_(1) {}

  const string operator()()
  {
    stringstream ss(stringstream::out);
    ss<<"rho"<<i_<<i_;
    i_++;
    return ss.str();
  }

private:
  size_t i_;
};

} // details

DiagonalDO::DiagonalDO(size_t dim) : 
  Base("Multilevel",list_of(string("rho00")).repeat_fun(dim-1,details::Helper())),
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

} // structure
} // averaged
