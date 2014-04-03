// Copyright András Vukics 2006–2014. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "ExampleInteraction.h"

basic::InteractionX_X::InteractionX_X(const PumpedLossyMode& m0, const PumpedLossyMode& m1, double g)
  : Interaction<2>(m0,m1,RF{"g",g,sqrt(m0.getDimension()*m1.getDimension())}),
    TridiagonalHamiltonian<2,false>(g*
                                    (aop(m0)+aop(m0).dagger())*
                                    (aop(m1)+aop(m1).dagger())
                                    )
{}


hierarchical::InteractionX_X::InteractionX_X(const ModeBase& m0, const ModeBase& m1, double g)
  : Interaction<2>(m0,m1,RF{"g",g,sqrt(m0.getDimension()*m1.getDimension())}),
    TridiagonalHamiltonian<2,true>(g*
                                   (aop(m0)+aop(m0).dagger())*
                                   (aop(m1)+aop(m1).dagger())
                                   )
{}


auto hierarchical::InteractionX_X_Correlations::average_v(NoTime, const LazyDensityOperator& matrix) const -> const Averages
{
  auto averages(initializedAverages());

  for (int n=0; n<int(matrix.getDimension(0)); n++) for (int m=1; m<int(matrix.getDimension(1)); m++) {
      if(n<int(matrix.getDimension(0))-1) {
        dcomp temp=sqrt(m*(n+1))*matrix(n,m)(n+1,m-1);
        averages(0)+=real(temp);
        averages(1)+=imag(temp);
      }
      if(n>0) {
        dcomp temp=sqrt(m*n)*matrix(n,m)(n-1,m-1);
        averages(2)+=real(temp);
        averages(3)+=imag(temp);
      }
    }

  double
    xq= averages(0)+averages(2),
    xp= averages(1)+averages(3),
    yq=-averages(1)+averages(3),
    yp= averages(0)-averages(2);

  averages=xq,xp,yq,yp;

  return averages;

}
