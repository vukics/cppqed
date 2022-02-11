// Copyright András Vukics 2006–2022. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
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


void performTransformation(const CArray<4>& trafo, const CArray<2>& in, CArray<2>& out)
{
  out=sum(sum(trafo(i,j,k,l)*in(k,l),l),k);
}

} // specializations

} // transformation

} // quantumdata
