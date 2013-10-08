#include "QbitModeCorrelations.h"

#include "impl/LazyDensityOperator.tcc"


QbitModeCorrelations::QbitModeCorrelations()
  : EA_Base("QbitModeCorrelations",{"real(<sigma*adagger>","imag(\")"})
{
}


const QbitModeCorrelations::Averages
QbitModeCorrelations::average_v(NoTime, const LazyDensityOperator& matrix) const
{
  typedef LazyDensityOperator::Idx Idx;

  Averages averages(6);
  averages=0;

  for (int n=0; n<int(matrix.getDimension(1))-1; n++) {
    dcomp temp=sqrt(n+1)*matrix(Idx(1,n),Idx(0,n+1));
    averages(0)+=real(temp);
    averages(1)+=imag(temp);
  }
  for (int n=0; n<int(matrix.getDimensions()[1])-1; n++) {
    dcomp temp=sqrt(n+1)*matrix(Idx(1,n+1),Idx(0,n));
    averages(2)+=real(temp);
    averages(3)+=imag(temp);
  }
  for (int n=0; n<int(matrix.getDimensions()[1])-1; n++) {
    dcomp temp=.5*sqrt(n+1)*(matrix(Idx(1,n+1),Idx(1,n))-matrix(Idx(0,n+1),Idx(0,n)));
    averages(4)+=real(temp);
    averages(5)+=imag(temp);
  }

  return averages;

}

