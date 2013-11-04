#include "MomentumCorrelation.h"
#include <boost/assign.hpp>

using namespace boost::assign;

MomentumCorrelation::MomentumCorrelation() :
   EA_Base("MomentumCorrelation",list_of("<P1P2>"))
{}

const MomentumCorrelation::Averages MomentumCorrelation::average_v(const LazyDensityOperator& matrix) const
{
  typedef LazyDensityOperator::Idx Idx;
  
  Averages averages(1);
  averages=0;
  double diag,temp;
  size_t dim1 = matrix.getDimensions()[0];
  size_t dim2 = matrix.getDimensions()[1];
  for (int n=0; n < dim1; n++)
    for (int m=0; m < dim2;m++) {
      averages(0)+=(-(dim1/2.)+n)*(-(dim2/2.)+m)*matrix(Idx(n,m));
    }
  return averages;
}

namespace momentumcorrelationinteraction {

Base::Base(particle::Ptr p1, particle::Ptr p2): 
   IA_Base(Frees(p1,p2),RealFreqs(),ComplexFreqs())
{

}

} // momentumcorrelationinteraction
