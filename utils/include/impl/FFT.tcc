// -*- C++ -*-
#ifndef   UTILS_INCLUDE_IMPL_FFT_TCC_INCLUDED
#define   UTILS_INCLUDE_IMPL_FFT_TCC_INCLUDED

#include "FFT.h"

#include "ArrayTraitsFwd.h"

#include<cstddef>

namespace fft {

template<typename A>
void transform(A& a, Direction dir) throw(FFT_Exception)
{
  typedef cpputils::ArrayTraversalTraits<A> Traits;
  details::transform(Traits::data(a),Traits::stride(a),Traits::ssLimit(a),dir);
}


} // fft

#endif // UTILS_INCLUDE_IMPL_FFT_TCC_INCLUDED
