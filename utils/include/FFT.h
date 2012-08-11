// -*- C++ -*-
#ifndef   _FAST_FOURIER_TRANSFORM_H
#define   _FAST_FOURIER_TRANSFORM_H

#include "FFTFwd.h"

#include "Exception.h"


namespace fft {

struct FFT_Exception : public cpputils::Exception {};

namespace details {

void transform(double*, size_t, size_t, Direction) throw(FFT_Exception);

} // details


template<typename A>
void transform(A&, Direction) throw(FFT_Exception);

} // fft

#endif // _FAST_FOURIER_TRANSFORM_H
