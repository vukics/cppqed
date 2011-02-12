#include "DensityOperator.h"

namespace quantumdata {

template<> 
inline
const dcomp 
DensityOperator<1>::operator()(const Idx& i, const Idx& j) const
{
  return operator()()(i,j);
}

} // quantumdata


