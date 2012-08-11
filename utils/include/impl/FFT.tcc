// -*- C++ -*-
#ifndef   _FAST_FOURIER_TRANSFORM_IMPL_H
#define   _FAST_FOURIER_TRANSFORM_IMPL_H

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

#endif // _FAST_FOURIER_TRANSFORM_IMPL_H
