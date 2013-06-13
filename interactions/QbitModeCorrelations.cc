#include "QbitModeCorrelations.h"

#include "impl/LazyDensityOperator.tcc"

#include <boost/assign/list_of.hpp>


QbitModeCorrelations::QbitModeCorrelations()
  : EA_Base(
	    "QbitModeCorrelations",
	    boost::assign::list_of("real(<sigma*adaggerr>")("imag(\")")
	    )
{
}


const QbitModeCorrelations::Averages
QbitModeCorrelations::average_v(const LazyDensityOperator& matrix) const
{
  typedef LazyDensityOperator::Idx Idx;

  Averages averages(2);
  averages=0;

  for (int n=0; n<int(matrix.getDimension(1))-1; n++) {
    dcomp temp=sqrt(n+1)*matrix(Idx(1,n),Idx(0,n+1));
    averages(0)+=real(temp);
    averages(1)+=imag(temp);
  }

  return averages;

}

