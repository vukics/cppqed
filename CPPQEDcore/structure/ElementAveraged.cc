#include "ElementAveraged.h"

#include "FormDouble.tcc"
#include "Range.h"

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>


using namespace std;



std::ostream& structure::displayCommon(const AveragedCommon::Averages& averages, std::ostream& os, int precision)
{
  using namespace formdouble;

  namespace bll=boost::lambda;

  os<<'\t';
  {
    const FormDouble fd(precision);
    boost::for_each(averages,os<<bll::bind(&FormDouble::operator()<double>,&fd,bll::_1));
  }
  return os;
}

