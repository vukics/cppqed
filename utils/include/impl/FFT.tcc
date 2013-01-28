// -*- C++ -*-
#ifndef   UTILS_INCLUDE_IMPL_FFT_TCC_INCLUDED
#define   UTILS_INCLUDE_IMPL_FFT_TCC_INCLUDED

#include "FFT.h"


namespace fft {

template<typename A>
void transform(A& a, Direction dir) throw(FFT_Exception)
{
  details::transform(cpputils::data(a),cpputils::stride(a),cpputils::subscriptLimit(a),dir);
}


} // fft

#endif // UTILS_INCLUDE_IMPL_FFT_TCC_INCLUDED
