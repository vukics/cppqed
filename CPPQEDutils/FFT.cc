// Copyright András Vukics 2006–2020. Distributed under the Boost Software License, Version 1.0. (See accompanying file LICENSE.txt)
#include "FFT.h"

#include<gsl/gsl_errno.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf.h>
#include<gsl/gsl_fft_complex.h>


namespace fft {


void details::transform(double* data, size_t stride, size_t n, Direction sign)
{
  if (gsl_fft_complex_radix2_transform(data,stride,n,(sign==DIR_XK ? gsl_fft_forward : gsl_fft_backward))!=GSL_SUCCESS) 
    throw FFT_Exception();  
}


} // fft
