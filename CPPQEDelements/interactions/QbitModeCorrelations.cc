// Copyright András Vukics 2006–2023. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "QbitModeCorrelations.h"

#include "LazyDensityOperator.h"


QbitModeCorrelations::QbitModeCorrelations()
  : EA_Base("QbitModeCorrelations",{"real(<sigma*adagger>","imag(\")","real(<sigma*a>)","imag(\")","real(<sigmaZ*a>)","imag(\")"})
{
}


const structure::Averages QbitModeCorrelations::average_v(NoTime, const quantumdata::LazyDensityOperator<2>& matrix) const
{
  auto averages(initializedAverages());

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

