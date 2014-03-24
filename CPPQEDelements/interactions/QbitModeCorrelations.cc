#include "QbitModeCorrelations.h"

#include "LazyDensityOperator.tcc"


QbitModeCorrelations::QbitModeCorrelations()
  : EA_Base("QbitModeCorrelations",{"real(<sigma*adagger>","imag(\")","real(<sigma*a>)","imag(\")","real(<sigmaZ*a>)","imag(\")"})
{
}


const QbitModeCorrelations::Averages
QbitModeCorrelations::average_v(NoTime, const LazyDensityOperator& matrix) const
{
  Averages averages(6);
  averages=0;

  for (int n=0; n<int(matrix.getDimension(1))-1; n++) {
    dcomp temp=sqrt(n+1)*matrix(1,n)(0,n+1);
    averages(0)+=real(temp);
    averages(1)+=imag(temp);
  }
  for (int n=0; n<int(matrix.getDimension(1))-1; n++) {
    dcomp temp=sqrt(n+1)*matrix(1,n+1)(0,n);
    averages(2)+=real(temp);
    averages(3)+=imag(temp);
  }
  for (int n=0; n<int(matrix.getDimension(1))-1; n++) {
    dcomp temp=.5*sqrt(n+1)*(matrix(1,n+1)(1,n)-matrix(0,n+1)(0,n));
    averages(4)+=real(temp);
    averages(5)+=imag(temp);
  }

  return averages;

}

