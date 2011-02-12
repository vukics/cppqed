// -*- C++ -*-
#ifndef   _FAST_FOURIER_TRANSFORM_H
#define   _FAST_FOURIER_TRANSFORM_H

#include "FFTFwd.h"

#include "Exception.h"


namespace fft {

struct FFT_Exception : public cpputils::Exception {};

template<typename A>
void transform(A&, Direction) throw(FFT_Exception);

} // fft

#include "impl/FFT.tcc"

#endif // _FAST_FOURIER_TRANSFORM_H
