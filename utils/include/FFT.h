// -*- C++ -*-
#ifndef   UTILS_INCLUDE_FFT_H_INCLUDED
#define   UTILS_INCLUDE_FFT_H_INCLUDED

#include "FFTFwd.h"

#include "Exception.h"

#include <cstddef>

namespace fft {

struct FFT_Exception : public cpputils::Exception {};

namespace details {

void transform(double*, size_t, size_t, Direction) throw(FFT_Exception);

} // details


template<typename A>
void transform(A&, Direction) throw(FFT_Exception);


inline const Direction reverse(Direction dir) {return dir==DIR_KX ? DIR_XK : DIR_KX;}

} // fft

#endif // UTILS_INCLUDE_FFT_H_INCLUDED
