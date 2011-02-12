#include "Transformation.h"


namespace quantumdata {

namespace transformation {


namespace specializations {

using namespace linalg;
using namespace blitz;
using namespace tensor;

void performTransformation(const CMatrix& trafo, const CVector& in, CVector& out)
{
  out=sum(trafo(i,j)*in(j),j);
}


void performTransformation(const TTD_CARRAY(4)& trafo, const TTD_CARRAY(2)& in, TTD_CARRAY(2)& out)
{
  out=sum(sum(trafo(i,j,k,l)*in(k,l),l),k);
}

} // specializations

} // transformation

} // quantumdata
