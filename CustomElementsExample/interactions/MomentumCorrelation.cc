// Copyright Raimar Sandner 2012–2017. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "MomentumCorrelation.h"
#include <boost/assign.hpp>

using namespace boost::assign;

MomentumCorrelation::MomentumCorrelation(particle::Ptr p0, particle::Ptr p1) :
  structure::Interaction<2>(Frees(p0,p1),RealFreqs(),ComplexFreqs()),
  EA_Base("MomentumCorrelation",{"<P1P2>"})
{}

const MomentumCorrelation::Averages MomentumCorrelation::average_v(Time, const LazyDensityOperator& matrix) const
{
  auto averages(initializedAverages());

  size_t dim1 = matrix.getDimensions()[0];
  size_t dim2 = matrix.getDimensions()[1];

  for (int n=0; n < dim1; n++)
    for (int m=0; m < dim2;m++) {
      averages(0)+=(-(dim1/2.)+n)*(-(dim2/2.)+m)*matrix(n,m);
    }
  return averages;
}
