#include "DensityOperator.h"

namespace quantumdata {

template<> 
inline
const dcomp 
DensityOperator<1>::index(const Idx& i, const Idx& j) const
{
  return getArray()(i,j);
}

} // quantumdata


