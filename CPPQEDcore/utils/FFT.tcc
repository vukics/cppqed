// -*- C++ -*-
#ifndef   CPPQEDCORE_UTILS_FFT_TCC_INCLUDED
#define   CPPQEDCORE_UTILS_FFT_TCC_INCLUDED

#include "FFT.h"
// #include "ArrayTraits.h"
// The same note applies as with EvolvedGSL.tcc

namespace fft {

template<typename A>
void transform(A& a, Direction dir)
{
  details::transform(cpputils::data(a),cpputils::stride(a),cpputils::subscriptLimit(a),dir);
}


} // fft

#endif // CPPQEDCORE_UTILS_FFT_TCC_INCLUDED
