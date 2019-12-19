// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
/// \briefFile{Fast Fourier transformation}
#ifndef   CPPQEDCORE_UTILS_FFT_H_INCLUDED
#define   CPPQEDCORE_UTILS_FFT_H_INCLUDED

#include "FFTFwd.h"

#include "Exception.h"

#include <cstddef>

/// Fast Fourier transformation
namespace fft {

/// direction of FFT
enum Direction {
  DIR_XK, ///< “forward” from x- to k-space (by convention)
  DIR_KX  ///< “backward” from k- to x-space
};


struct FFT_Exception : public cpputils::Exception {};

namespace details {

void transform(double*, size_t, size_t, Direction);

} // details


/// Radix-2 FFT transformation for complex data
/**
 * Implemented via [GSL](http://www.gnu.org/software/gsl/manual/html_node/Radix_002d2-FFT-routines-for-complex-data.html).
 * 
 * In the spirit of evolved::Evolved it can take any array type that is transformable to a 
 * [complex packed array](https://www.gnu.org/software/gsl/manual/html_node/Overview-of-complex-data-FFTs.html)
 * via the \link ArrayTraits.h array traits\endlink functions.
 * 
 * \tparam A the complex array type
 * 
 * \throw FFT_Exception if the size of the array is not a power of two.
 * 
 */
template<typename A>
void transform(A&, Direction);


/// reverses the direction
inline const Direction reverse(Direction dir) {return dir==DIR_KX ? DIR_XK : DIR_KX;}

} // fft

#endif // CPPQEDCORE_UTILS_FFT_H_INCLUDED
