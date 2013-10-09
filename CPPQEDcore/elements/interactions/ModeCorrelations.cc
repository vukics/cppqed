#include "ModeCorrelations.h"

#include "LazyDensityOperator.tcc"
#include "Structure.h"

#include "SmartPtr.h"

using namespace boost;
using namespace assign;
using cpputils::sharedPointerize;


ModeCorrelations::ModeCorrelations()
  : EA_Base("ModeCorrelations",{"<a number operator>","VAR(a number operator)","real(<a>)","imag(\")","<X^2>-<X>^2","<Y^2>-<Y>^2","<(XY+YX)/2>-<X><Y>",
                                "<b number operator>","VAR(b number operator)","real(<b>)","imag(\")","<Q^2>-<Q>^2","<P^2>-<P>^2","<(QP+PQ)/2>-<Q><P>",
                                "<XQ>-<X><Q>","<XP>-<X><P>","<YQ>-<Y><Q>","<YP>-<Y><P>"}
           ),
    averagedMode_()
{
}


namespace {

#include "details/BinaryHelper.h"

}


const ModeCorrelations::Averages
ModeCorrelations::average_v(structure::NoTime, const LazyDensityOperator& matrix) const
{
  using quantumdata::partialTrace;
  typedef LazyDensityOperator::Idx Idx;

  Averages averages(18);
  averages=0;

  {
    const Averaged1::Ptr p(sharedPointerize(averagedMode_));

    const Averages 
      a0(partialTrace<V0,Averages>(matrix,bind(structure::average<1>,p,0,_1))),
      a1(partialTrace<V1,Averages>(matrix,bind(structure::average<1>,p,0,_1)));

    copy(a1,copy(a0,averages.begin()));
  }

  for (int n=0; n<int(matrix.getDimension(0)); n++) for (int m=1; m<int(matrix.getDimension(1)); m++) {
      if(n<int(matrix.getDimension(0))-1) {
        dcomp temp=sqrt(m*(n+1))*matrix(Idx(n,m),Idx(n+1,m-1));
        averages(14)+=real(temp);
        averages(15)+=imag(temp);
      }
      if(n>0) {
        dcomp temp=sqrt(m*n)*matrix(Idx(n,m),Idx(n-1,m-1));
        averages(16)+=real(temp);
        averages(17)+=imag(temp);
      }
    }

  return averages;

}


void
ModeCorrelations::process_v(Averages& averages) const
{
  {
    Averages ranged(averages(blitz::Range(0, 6)));
    averagedMode_.process(ranged);
  }
  {
    Averages ranged(averages(blitz::Range(7,13)));
    averagedMode_.process(ranged);
  }

  double 
    xAvr=sqrt(2)*averages( 2),
    yAvr=sqrt(2)*averages( 3),
    qAvr=sqrt(2)*averages( 9),
    pAvr=sqrt(2)*averages(10);

  double
    xq= averages(14)+averages(16)-xAvr*qAvr,
    xp= averages(15)+averages(17)-xAvr*pAvr,
    yq=-averages(15)+averages(17)-yAvr*qAvr,
    yp= averages(14)-averages(16)-yAvr*pAvr;

  averages(14)=xq;
  averages(15)=xp;
  averages(16)=yq;
  averages(17)=yp;

}
