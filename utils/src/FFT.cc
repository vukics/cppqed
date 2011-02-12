#include "FFT.h"

#include<gsl/gsl_errno.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf.h>
#include<gsl/gsl_fft_complex.h>


namespace fft {


namespace details {

void transform(double* data, size_t stride, size_t n, Direction sign) throw(FFT_Exception)
{
  if (gsl_fft_complex_radix2_transform(data,stride,n,(sign==FFTDIR_XK ? forward : backward))!=GSL_SUCCESS) 
    throw FFT_Exception();  
}


} // details


} // fft
