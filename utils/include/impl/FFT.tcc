// -*- C++ -*-
#ifndef   _FAST_FOURIER_TRANSFORM_IMPL_H
#define   _FAST_FOURIER_TRANSFORM_IMPL_H

#include "ArrayTraitsFwd.h"

#include<cstddef>

namespace fft {


namespace details {

void transform(double*, size_t, size_t, Direction) throw(FFT_Exception);

} // details


template<typename A>
void transform(A& a, Direction dir) throw(FFT_Exception)
{
  typedef cpputils::ArrayTraversalTraits<A> Traits;
  details::transform(Traits::data(a),Traits::stride(a),Traits::ssLimit(a),dir);
}


} // fft

#endif // _FAST_FOURIER_TRANSFORM_IMPL_H
